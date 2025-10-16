from flask import Blueprint, request, jsonify
from rdkit import Chem
from rdktools import smarts
from rdkit.Chem import FilterCatalog
from utils.request_processing import process_smiles_input, process_smarts_input
import tempfile

smarts_filter = Blueprint("smarts_filter", __name__, url_prefix="/smarts_filter")


class MoleculeCollector:
    def __init__(self):
        self.accepted = []
        self.accepted_names = []

    def SetProps(self, props):
        pass

    def write(self, mol):
        canon = Chem.MolToSmiles(mol, canonical=True)
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        self.accepted.append(canon)
        self.accepted_names.append([name, canon])  # Use list, not tuple


class MatchCountCollector:
    def __init__(self):
        self.results = []

    def SetProps(self, props):
        pass

    def write(self, mol):
        self.results.append({
            "smiles": Chem.MolToSmiles(mol, canonical=True),
            "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
            "n_matches": mol.GetIntProp("n_matches")
        })


class MatchMultiCountCollector:
    def __init__(self, named_smarts):
        self.named_smarts = named_smarts
        # Compile RDKit query molecules once for efficiency
        self.smarts_queries = []
        for pattern, name in named_smarts:
            qmol = Chem.MolFromSmarts(pattern)
            if not qmol:
                raise ValueError(f"Invalid SMARTS pattern: {pattern}")
            self.smarts_queries.append((qmol, pattern, name))

        self.results = {}

    def SetProps(self, props):
        # Not used here, but kept for compatibility
        pass

    def write(self, mol):
        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        if mol_name not in self.results:
            self.results[mol_name] = {}

        for qmol, pattern, name in self.smarts_queries:
            matches = mol.GetSubstructMatches(qmol)
            count = len(matches)
            highlight_atoms = [list(match) for match in matches] if matches else []
            self.results[mol_name][name] = {
                "count": count,
                "highlight_atoms": highlight_atoms,
            }


class MatchFilter:
    def __init__(self):
        self.accepted = []

    def SetProps(self, props):
        pass

    def write(self, mol):
        self.accepted.append({
            "SMILES": Chem.MolToSmiles(mol, canonical=True),
            "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
        })


# Initialize the PAINS filter catalog once (reuse on every request)
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
catalog = FilterCatalog.FilterCatalog(params)


def process_request_data(request, is_get_request=True):
    """Helper function to process request data from both GET and POST requests"""
    if is_get_request:
        # GET request - use query parameters
        smiles_list = process_smiles_input(request, "SMILES", 50)  # Limited to 50 for GET
        names_list = process_smiles_input(request, "Smile_Names", 50)
        return smiles_list, names_list
    else:
        # POST request - use JSON body
        json_data = request.get_json()
        if not json_data:
            return None, None
        
        smiles_list = json_data.get("SMILES", [])
        names_list = json_data.get("Smile_Names", [])
        
        if not isinstance(smiles_list, list):
            return None, None
            
        return smiles_list, names_list


def get_request_params(request, is_get_request=True):
    """Helper function to get parameters from both GET and POST requests"""
    if is_get_request:
        return request.args
    else:
        json_data = request.get_json()
        return json_data if json_data else {}


def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if not value:
        return False
    return value.lower() in ['true', '1', 'yes', 'y']


@smarts_filter.route('/get_filterpains', methods=['GET', 'POST'])
def get_filterpains():
    is_get_request = request.method == 'GET'
    smiles_list, names_list = process_request_data(request, is_get_request)
    
    if smiles_list is None:
        return jsonify({"error": "Invalid request format"}), 400

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    # Limit check for GET requests
    if is_get_request and len(smiles_list) > 50:
        return jsonify({"error": "GET requests are limited to 50 SMILES for testing purposes"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({
            "error": f"'SMILES' and 'Names' must be the same length: got {len(smiles_list)} and {len(names_list)}"
        }), 400

    if not names_list:
        names_list = smiles_list

    parsed = []
    invalid = []

    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append((name, smiles, mol))
        else:
            invalid.append({"name": name, "smiles": smiles})

    if not parsed:
        return jsonify({"error": "No valid molecules found"}), 400

    # Get all PAINS filters
    all_pains_filters = []
    num_entries = catalog.GetNumEntries()
    for i in range(num_entries):
        entry = catalog.GetEntry(i)
        all_pains_filters.append(entry.GetDescription())

    results = []

    for name, smiles, mol in parsed:
        reasons = []
        highlight_atom_sets = []

        try:
            matched_entries = catalog.GetMatches(mol)
            for entry in matched_entries:
                reasons.append(entry.GetDescription())

                try:
                    matches = entry.GetFilterMatches(mol)
                    for match in matches:
                        atom_indices = [atom_idx for _, atom_idx in match.atomPairs]
                        if atom_indices:
                            highlight_atom_sets.append(atom_indices)
                except Exception as match_err:
                    print(f"[Match Error] {entry.GetDescription()} on '{name}': {match_err}")
        except Exception as e:
            print(f"[Catalog Match Error] Failed on molecule '{name}': {e}")

        results.append({
            "name": name,
            "smiles": smiles,
            "failed": bool(reasons),
            "reasons": reasons,
            "highlight_atoms": highlight_atom_sets
        })

    return jsonify({
        "results": results,
        "all_pains_filters": all_pains_filters,
        "invalid": invalid
    }), 200


@smarts_filter.route('/get_matchcounts', methods=['GET', 'POST'])
def get_matchcounts():
    is_get_request = request.method == 'GET'
    smiles_list, names_list = process_request_data(request, is_get_request)
    params = get_request_params(request, is_get_request)
    
    if smiles_list is None:
        return jsonify({"error": "Invalid request format"}), 400

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    # Limit check for GET requests
    if is_get_request and len(smiles_list) > 50:
        return jsonify({"error": "GET requests are limited to 50 SMILES for testing purposes"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({"error": f"'SMILES' and 'Names' must be the same length: got {len(smiles_list)} and {len(names_list)}"}), 400

    if not names_list:
        names_list = smiles_list

    if is_get_request:
        smart = process_smarts_input(request)[0]
    else:
        smart = params.get("smarts", "")
        if isinstance(smart, list) and smart:
            smart = smart[0]
    
    if not smart:
        return jsonify({"error": "Smarts pattern is required"}), 400

    exclude_mol_props = params.get("ExcludeMolProp", False)
    usa = params.get("usa", False)
    nonzero_rows = params.get("nonzero_rows", False)
    
    if is_get_request:
        exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)
        usa = request.args.get("usa", type=bool, default=False)
        nonzero_rows = request.args.get("nonzero_rows", type=bool, default=False)
    
    parsed = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append(mol)

    if not parsed:
        return jsonify({"error": "No valid molecules parsed from input"}), 400

    collector = MatchCountCollector()

    smarts.MatchCounts(
        smarts=smart,
        usa=usa,
        molReader=parsed,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
        nonzero_rows=nonzero_rows,
    )

    return jsonify(collector.results), 200


@smarts_filter.route('/get_matchfilter', methods=['GET', 'POST'])
def get_matchfilter():
    is_get_request = request.method == 'GET'
    smiles_list, names_list = process_request_data(request, is_get_request)
    params = get_request_params(request, is_get_request)
    
    if smiles_list is None:
        return jsonify({"error": "Invalid request format"}), 400

    # Limit check for GET requests
    if is_get_request and len(smiles_list) > 50:
        return jsonify({"error": "GET requests are limited to 50 SMILES for testing purposes"}), 400

    if is_get_request:
        smart = process_smarts_input(request)[0]
    else:
        smart = params.get("smarts", "")
        if isinstance(smart, list) and smart:
            smart = smart[0]
    
    smart_name = params.get("Smart_Names", smart)
    exclude_mol_props = params.get("exclude_mol_props", False)
    
    if is_get_request:
        exclude_mol_props = request.args.get("exclude_mol_props", type=bool, default=False)

    if len(smiles_list) != len(names_list):
        return jsonify({"error": "SMILES and Names must be same length"}), 400

    parsed = []
    invalid = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append((name, smiles, mol))
        else:
            invalid.append({"name": name, "smiles": ""})

    collector = MatchFilter()
    smarts.MatchFilter(
        smarts=smart,
        molReader=[mol for _, _, mol in parsed],
        molWriter=collector,
        exclude_mol_props = exclude_mol_props
    )

    accepted_smiles = set(item['SMILES'] for item in collector.accepted)
    failed = []
    passed = []
    qmol = Chem.MolFromSmarts(smart)

    for name, smiles, mol in parsed:
        canon = Chem.MolToSmiles(mol, canonical=True)
        matches = mol.GetSubstructMatches(qmol)
        if canon in accepted_smiles:
            failed.append({
                "name": name,
                "smiles": smiles,
                "failed": True,
                "reason": smart_name,
                "highlight_atoms": [list(m) for m in matches]
            })
        else:
            passed.append({
                "name": name,
                "smiles": smiles,
                "failed": False
            })

    failed.extend(invalid)
    return jsonify({"passed": passed, "failed": failed}), 200


@smarts_filter.route('/get_multi_matchcounts', methods=['GET', 'POST'])
def get_multi_matchcounts():
    is_get_request = request.method == 'GET'
    smiles_list, names_list = process_request_data(request, is_get_request)
    params = get_request_params(request, is_get_request)
    
    if smiles_list is None:
        return jsonify({"error": "Invalid request format"}), 400

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    # Limit check for GET requests
    if is_get_request and len(smiles_list) > 50:
        return jsonify({"error": "GET requests are limited to 50 SMILES for testing purposes"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({
            "error": f"'SMILES' and 'Smile_Names' must be the same length: got {len(smiles_list)} and {len(names_list)}"
        }), 400

    if not names_list:
        names_list = smiles_list

    try:
        if is_get_request:
            smarts_list = request.args.getlist("smarts")
            smart_names = request.args.getlist("Smart_Names")
            isomeric_smiles = str_to_bool(request.args.get("isomericSmiles", 'false'))
            kekule_smiles = str_to_bool(request.args.get("kekuleSmiles", 'false'))
        else:
            smarts_list = params.get("smarts", [])
            smart_names = params.get("Smart_Names", [])
            isomeric_smiles = str_to_bool(params.get("isomericSmiles", 'false'))
            kekule_smiles = str_to_bool(params.get("kekuleSmiles", 'false'))
    except Exception as e:
        return jsonify({"error": str(e)}), 400

    if not smarts_list or not smart_names:
        return jsonify({"error": "Missing SMARTS patterns or Smart_Names"}), 400

    if len(smarts_list) != len(smart_names):
        return jsonify({
            "error": f"'Smart_Names' and 'smarts' must be the same length: got {len(smart_names)} and {len(smarts_list)}"
        }), 400

    named_smarts = list(zip(smarts_list, smart_names))

    # Write SMARTS and names to temporary file
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_smarts_file:
        for pattern, name in named_smarts:
            temp_smarts_file.write(f"{pattern} {name}\n")
        temp_smarts_file_path = temp_smarts_file.name

    exclude_mol_props = str_to_bool(params.get("ExcludeMolProp", 'false'))
    usa = str_to_bool(params.get("unique_set", 'false'))
    raiseError = str_to_bool(params.get("strict_error", 'false'))
    nonzero_rows = str_to_bool(params.get('only_rows', 'false'))

    parsed = []
    name_to_mol = {}
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append(mol)
            name_to_mol[name] = mol

    if not parsed:
        return jsonify({"error": "No valid molecules parsed from input"}), 400

    def mol_to_custom_smiles(mol):
        if kekule_smiles:
            return Chem.MolToSmiles(mol, kekuleSmiles=True)
        if isomeric_smiles:
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        return Chem.MolToSmiles(mol)    

    # Run the match counts
    collector = MatchMultiCountCollector(named_smarts)
    smarts.MatchCountsMulti(
        smarts_file_path=temp_smarts_file_path,
        strict_smarts=raiseError,
        usa=usa,
        molReader=parsed,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
        nonzero_rows=nonzero_rows
    )

    results = []
    for mol_name in collector.results:
        mol = name_to_mol.get(mol_name)
        if mol is None:
            continue  # safety

        smiles_out = mol_to_custom_smiles(mol)
        result = {
            "name": mol.GetProp("_Name"),
            "smiles": smiles_out,
            "matches": []
        }

        for pattern, name in named_smarts:
            match_info = collector.results.get(mol_name, {}).get(name, {"count": 0, "highlight_atoms": []})
            result["matches"].append({
                "smarts": pattern,
                "name": name,
                "count": match_info["count"],
                "highlight_atoms": match_info["highlight_atoms"]
            })

        results.append(result)

    return jsonify(results), 200


@smarts_filter.route('/get_multi_matchfilter', methods=['GET', 'POST'])
def get_multi_matchfilter():
    is_get_request = request.method == 'GET'
    smiles_list, names_list = process_request_data(request, is_get_request)
    params = get_request_params(request, is_get_request)
    
    if smiles_list is None:
        return jsonify({"error": "Invalid request format"}), 400

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    # Limit check for GET requests
    if is_get_request and len(smiles_list) > 50:
        return jsonify({"error": "GET requests are limited to 50 SMILES for testing purposes"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({"error": "'SMILES' and 'Names' must be the same length"}), 400

    if not names_list:
        names_list = smiles_list

    if is_get_request:
        smarts_list = request.args.getlist("smarts")
        smart_names = request.args.getlist("Smart_Names")
    else:
        smarts_list = params.get("smarts", [])
        smart_names = params.get("Smart_Names", [])

    if not smarts_list or not smart_names:
        return jsonify({"error": "Missing SMARTS or Smart_Names"}), 400

    if len(smarts_list) != len(smart_names):
        return jsonify({"error": "'Smart_Names' and 'smarts' must be the same length"}), 400

    exclude_mol_props = params.get("exclude_mol_props", False)
    strict = params.get("strict", False)
    
    if is_get_request:
        exclude_mol_props = request.args.get("exclude_mol_props", type=bool, default=False)
        strict = request.args.get("strict", type=bool, default=False)

    # Write SMARTS file for MatchFilterMulti
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_smarts_file:
        for smarts_pattern in smarts_list:
            temp_smarts_file.write(smarts_pattern.strip() + "\n")
        temp_smarts_file_path = temp_smarts_file.name

    parsed = []
    invalid = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append((name, smiles, mol))
        else:
            invalid.append({
                "name": name,
                "smiles": "",
                "failed": True,
                "reason": "Invalid SMILES",
                "highlight_atoms": []
            })

    if not parsed and not invalid:
        return jsonify({"error": "No valid molecules parsed"}), 400

    collector = MatchFilter()
    smarts.MatchFilterMulti(
        smarts_file_path=temp_smarts_file_path,
        strict_smarts=strict,
        molReader=[mol for _, _, mol in parsed],
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    accepted_smiles = set(item['SMILES'] for item in collector.accepted)

    failed = []
    passed = []

    # For each molecule, test which SMARTS it failed
    for name, smiles, mol in parsed:
        canon = Chem.MolToSmiles(mol, canonical=True)

        if canon in accepted_smiles:
            # Failed: find which SMARTS matched (may be more than one)
            for smarts_pattern, smart_name in zip(smarts_list, smart_names):
                qmol = Chem.MolFromSmarts(smarts_pattern)
                if not qmol:
                    continue
                matches = mol.GetSubstructMatches(qmol)
                if matches:
                    failed.append({
                        "name": name,
                        "smiles": smiles,
                        "failed": True,
                        "reason": smart_name,
                        "highlight_atoms": [list(m) for m in matches]
                    })
                    break  # stop after first match for same mol
        else:
            passed.append({
                "name": name,
                "smiles": smiles,
                "failed": False
            })

    failed.extend(invalid)

    return jsonify({"passed": passed, "failed": failed}), 200
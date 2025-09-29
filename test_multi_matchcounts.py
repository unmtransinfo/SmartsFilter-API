import requests
import json
import time
import random

def generate_test_smiles(count=100):
    """Generate a list of diverse SMILES for testing"""
    base_smiles = [
        "CCO",  # ethanol
        "CCC",  # propane
        "c1ccccc1",  # benzene
        "CCN",  # ethylamine
        "CC(C)C",  # isobutane
        "CCCCC",  # pentane
        "c1ccc2ccccc2c1",  # naphthalene
        "c1ccncc1",  # pyridine
        "CC(=O)O",  # acetic acid
        "CCO[CH2]",  # ether
        "c1ccc(cc1)O",  # phenol
        "CC(C)(C)O",  # tert-butanol
        "c1ccc(cc1)N",  # aniline
        "CC(=O)N",  # acetamide
        "c1ccc2c(c1)cccc2",  # naphthalene alt
        "CC(C)CO",  # isobutanol
        "c1ccc(cc1)C",  # toluene
        "CCc1ccccc1",  # ethylbenzene
        "CC(C)CCO",  # isopentanol
        "c1ccc(cc1)CC",  # ethyltoluene
        "CCCCCO",  # pentanol
        "c1ccc2c(c1)ccc2",  # naphthalene variation
        "CC(C)C(C)C",  # 2,3-dimethylbutane
        "c1ccc(cc1)CCC",  # propylbenzene
        "CCCCCCCO",  # heptanol
        "c1ccc(cc1)CCCC",  # butylbenzene
        "CC(C)C(C)(C)C",  # 2,2,3-trimethylbutane
        "c1ccc2c(c1)cccc2C",  # methylnaphthalene
        "CCCCCCCCO",  # octanol
        "c1ccc(cc1)CCCCC",  # pentylbenzene
        "O=C1NC(=O)C(=Cc2ccc(Cl)cc2)C(=O)N1",  # PAINS compound
        "O=C1NC(=S)NC(=O)C1=Cc1ccc(Cl)cc1",  # Another PAINS compound
        "c1ccc2c(c1)c(c(n2)C)C(=O)O",  # indole derivative
        "CC(C)(C)c1ccc(cc1)O",  # BHT-like structure
        "c1ccc(cc1)S(=O)(=O)N",  # sulfonamide
        "CC(=O)Nc1ccc(cc1)O",  # acetaminophen-like
        "c1ccc(cc1)C(=O)O",  # benzoic acid
        "CC(C)Nc1ccccc1",  # N-isopropylaniline
        "c1ccc2c(c1)nc(n2)N",  # aminobenzimidazole
        "CC(=O)c1ccc(cc1)N",  # aminoacetophenone
        "c1ccc(cc1)C(=O)N",  # benzamide
        "CC(C)c1ccc(cc1)N",  # isopropylaniline
        "c1ccc2c(c1)c(cn2)C",  # methylindole
        "CC(=O)Nc1ccccc1",  # acetanilide
        "c1ccc(cc1)C(=O)C",  # acetophenone
        "CC(C)Oc1ccccc1",  # isopropyl phenyl ether
        "c1ccc2c(c1)cc(n2)C",  # methylindole variation
        "CC(=O)c1ccc(cc1)O",  # hydroxyacetophenone
        "c1ccc(cc1)C(=O)CC",  # propiophenone
        "CC(C)c1ccc(cc1)C",  # cymene
    ]
    
    # Extend the list by creating variations
    extended_smiles = []
    for i in range(count):
        if i < len(base_smiles):
            extended_smiles.append(base_smiles[i])
        else:
            # Create variations by modifying existing SMILES
            base = base_smiles[i % len(base_smiles)]
            # Simple modifications - add methyl groups, change chain lengths, etc.
            if "CC" in base and random.random() > 0.5:
                extended_smiles.append(base.replace("CC", "CCC", 1))
            elif "c1ccccc1" in base and random.random() > 0.5:
                extended_smiles.append(base.replace("c1ccccc1", "c1ccc(cc1)C"))
            else:
                extended_smiles.append(base)
    
    return extended_smiles[:count]

def generate_test_smarts(count=200):
    """Generate a list of diverse SMARTS patterns for testing"""
    base_smarts = [
        "[OH]",  # hydroxyl
        "[NH2]",  # amino
        "[CH3]",  # methyl
        "c1ccccc1",  # benzene ring
        "[C]=O",  # carbonyl
        "[S]",  # sulfur
        "[N]",  # nitrogen
        "[O]",  # oxygen
        "[Cl]",  # chlorine
        "[Br]",  # bromine
        "[F]",  # fluorine
        "[I]",  # iodine
        "[P]",  # phosphorus
        "C=C",  # alkene
        "C#C",  # alkyne
        "C-C",  # alkane
        "[nH]",  # NH in aromatic
        "c1ccc2ccccc2c1",  # naphthalene
        "[NH1]",  # secondary amine
        "[OH1]",  # hydroxyl explicit
        "C(=O)O",  # carboxylic acid
        "C(=O)N",  # amide
        "S(=O)(=O)",  # sulfone
        "S(=O)",  # sulfoxide
        "C(=O)C",  # ketone
        "[CH2]",  # methylene
        "[CH]",  # methine
        "[C+]",  # positive carbon
        "[C-]",  # negative carbon
        "[N+]",  # positive nitrogen
        "[N-]",  # negative nitrogen
        "[O-]",  # negative oxygen
        "c1ccncc1",  # pyridine
        "c1cccnc1",  # pyridine alt
        "c1ccco1",  # furan
        "c1cccs1",  # thiophene
        "c1cc[nH]c1",  # pyrrole
        "c1ccc2[nH]c3ccccc3c2c1",  # carbazole
        "c1ccc2c(c1)cccc2",  # naphthalene variation
        "C1CCC2CCCCC2C1",  # decalin
        "C1CCCCC1",  # cyclohexane
        "C1CCCC1",  # cyclopentane
        "C1CCC1",  # cyclobutane
        "C1CC1",  # cyclopropane
        "[R]",  # ring atom
        "[!R]",  # non-ring atom
        "[R1]",  # size 1 ring
        "[R2]",  # size 2 ring
        "[r3]",  # 3-membered ring
        "[r4]",  # 4-membered ring
        "[r5]",  # 5-membered ring
        "[r6]",  # 6-membered ring
        "[r7]",  # 7-membered ring
        "[r8]",  # 8-membered ring
        "[D1]",  # degree 1
        "[D2]",  # degree 2
        "[D3]",  # degree 3
        "[D4]",  # degree 4
        "[X1]",  # connectivity 1
        "[X2]",  # connectivity 2
        "[X3]",  # connectivity 3
        "[X4]",  # connectivity 4
        "[H0]",  # no hydrogens
        "[H1]",  # 1 hydrogen
        "[H2]",  # 2 hydrogens
        "[H3]",  # 3 hydrogens
        "[v1]",  # valence 1
        "[v2]",  # valence 2
        "[v3]",  # valence 3
        "[v4]",  # valence 4
        "[#6]",  # carbon by atomic number
        "[#7]",  # nitrogen by atomic number
        "[#8]",  # oxygen by atomic number
        "[#16]",  # sulfur by atomic number
        "[#15]",  # phosphorus by atomic number
        "[#9]",  # fluorine by atomic number
        "[#17]",  # chlorine by atomic number
        "[#35]",  # bromine by atomic number
        "[#53]",  # iodine by atomic number
        "[$([OH])]",  # recursive SMARTS for OH
        "[$([NH2])]",  # recursive SMARTS for NH2
        "[$(C=O)]",  # recursive SMARTS for C=O
        "[$(C(=O)O)]",  # recursive SMARTS for carboxylic acid
        "[$(C(=O)N)]",  # recursive SMARTS for amide
        "[$(S(=O)(=O))]",  # recursive SMARTS for sulfone
        "C=C-C",  # specific alkene pattern
        "C=C-C=C",  # conjugated system
        "C-C-C",  # propyl chain
        "C-C-C-C",  # butyl chain
        "C-C-C-C-C",  # pentyl chain
        "N-C-C-N",  # diamine pattern
        "O-C-C-O",  # diol pattern
        "c-C-c",  # aromatic-aliphatic-aromatic
        "c-C-C-c",  # longer aromatic-aliphatic bridge
        "c-N-c",  # aromatic amine
        "c-O-c",  # aromatic ether
        "c-S-c",  # aromatic thioether
        "C(=O)-C(=O)",  # dicarbonyl
        "N=N",  # azo group
        "N-N",  # hydrazine
        "O-O",  # peroxide
        "S-S",  # disulfide
        "[C]1[C][C][C][C][C]1",  # cyclohexane explicit
        "[c]1[c][c][c][c][c]1",  # benzene explicit
        "*-*",  # any bond
        "[*]",  # any atom
        "[!C]",  # not carbon
        "[!H]",  # not hydrogen
        "[C,N]",  # carbon or nitrogen
        "[C,N,O]",  # carbon, nitrogen, or oxygen
        "[C&H1]",  # carbon with 1 hydrogen
        "[C&H2]",  # carbon with 2 hydrogens
        "[C&H3]",  # carbon with 3 hydrogens
        "[C&R]",  # carbon in ring
        "[C&!R]",  # carbon not in ring
        "[N&H0]",  # nitrogen with no hydrogens
        "[N&H1]",  # nitrogen with 1 hydrogen
        "[N&H2]",  # nitrogen with 2 hydrogens
        "[O&H0]",  # oxygen with no hydrogens
        "[O&H1]",  # oxygen with 1 hydrogen
        "C~C",  # any bond between carbons
        "C:C",  # aromatic bond
        "C-C-C-C-C-C",  # hexyl chain
        "C1CC2CCC1C2",  # bridged ring system
        "C1CCC(CC1)C",  # substituted cyclohexane
        "c1ccc(cc1)c2ccccc2",  # biphenyl
    ]
    
    # Extend the list with variations
    extended_smarts = []
    for i in range(count):
        if i < len(base_smarts):
            extended_smarts.append(base_smarts[i])
        else:
            # Create variations
            base = base_smarts[i % len(base_smarts)]
            extended_smarts.append(base)
    
    return extended_smarts[:count]

def test_multi_matchcounts_large_scale():
    """Test the get_multi_matchcounts endpoint with 100+ molecules and 200+ SMARTS"""
    
    url = "http://localhost:11412/api/v1/smarts_filter/get_multi_matchcounts"
    
    # Generate test data
    print("Generating test data...")
    smiles_list = generate_test_smiles(150)  # 150 molecules
    names_list = [f"mol_{i+1}" for i in range(len(smiles_list))]
    smarts_list = generate_test_smarts(250)  # 250 SMARTS patterns
    smart_names = [f"smarts_{i+1}" for i in range(len(smarts_list))]
    
    print(f"Generated {len(smiles_list)} SMILES and {len(smarts_list)} SMARTS patterns")
    
    # Prepare request data with ALL parameters
    data = {
        "SMILES": smiles_list,
        "Smile_Names": names_list,
        "smarts": smarts_list,
        "Smart_Names": smart_names,
        "isomericSmiles": True,  # Test isomeric SMILES
        "kekuleSmiles": True,   # Test kekule SMILES (will override isomeric)
        "unique_set": True,     # Test unique set-of-atoms match counts
        "strict_error": False,  # Don't raise errors on invalid SMARTS
        "only_rows": True,      # Only include rows with at least one match
        "ExcludeMolProp": True  # Exclude molecular properties
    }
    
    print(f"Request payload size: ~{len(json.dumps(data))} characters")
    
    # Send request and measure time
    print("Sending POST request...")
    start_time = time.time()
    
    try:
        response = requests.post(url, json=data, timeout=300)  # 5 minute timeout
        end_time = time.time()
        
        print(f"Response received in {end_time - start_time:.2f} seconds")
        print(f"Status Code: {response.status_code}")
        
        if response.status_code == 200:
            result = response.json()
            print(f"Success! Processed {len(result)} molecules")
            
            # Print parameter configuration used
            print(f"Parameters used:")
            print(f"  - isomericSmiles: True")
            print(f"  - kekuleSmiles: True (takes precedence over isomeric)")
            print(f"  - unique_set: True")
            print(f"  - strict_error: False")
            print(f"  - only_rows: True (nonzero rows only)")
            print(f"  - ExcludeMolProp: True")
            
            # Analyze results
            total_matches = 0
            molecules_with_matches = 0
            
            for mol_result in result:
                mol_matches = sum(match['count'] for match in mol_result['matches'])
                total_matches += mol_matches
                if mol_matches > 0:
                    molecules_with_matches += 1
            
            print(f"Total matches found: {total_matches}")
            print(f"Molecules with at least one match: {molecules_with_matches}/{len(result)}")
            print(f"Average matches per molecule: {total_matches/len(result):.2f}")
            
            # Show sample results
            print("\nSample results (first 3 molecules):")
            for i, mol_result in enumerate(result[:3]):
                print(f"  Molecule {i+1}: {mol_result['name']}")
                print(f"    SMILES: {mol_result['smiles']}")
                matches_count = sum(match['count'] for match in mol_result['matches'])
                print(f"    Total matches: {matches_count}")
                
                # Show first few matches
                positive_matches = [m for m in mol_result['matches'] if m['count'] > 0]
                if positive_matches:
                    print(f"    First few matching patterns:")
                    for match in positive_matches[:3]:
                        print(f"      - {match['name']}: {match['count']} matches")
                print()
            
            return True
            
        else:
            print(f"Error: {response.status_code}")
            try:
                error_data = response.json()
                print(f"Error details: {error_data}")
            except:
                print(f"Response text: {response.text}")
            return False
            
    except requests.exceptions.Timeout:
        print("Request timed out after 5 minutes")
        return False
    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error: {e}")
        return False

def test_performance_comparison():
    """Compare performance between different data sizes"""
    sizes = [
        (10, 20),   # Small
        (50, 100),  # Medium  
        (100, 200), # Large
        (150, 250)  # Extra Large
    ]
    
    url = "http://localhost:11412/api/v1/smarts_filter/get_multi_matchcounts"
    
    for mol_count, smarts_count in sizes:
        print(f"\nTesting with {mol_count} molecules and {smarts_count} SMARTS...")
        
        smiles_list = generate_test_smiles(mol_count)
        names_list = [f"mol_{i+1}" for i in range(len(smiles_list))]
        smarts_list = generate_test_smarts(smarts_count)
        smart_names = [f"smarts_{i+1}" for i in range(len(smarts_list))]
        
        data = {
            "SMILES": smiles_list,
            "Smile_Names": names_list,
            "smarts": smarts_list,
            "Smart_Names": smart_names,
            "ExcludeMolProp": True
        }
        
        start_time = time.time()
        try:
            response = requests.post(url, json=data, timeout=300)
            end_time = time.time()
            
            if response.status_code == 200:
                result = response.json()
                total_comparisons = mol_count * smarts_count
                processing_time = end_time - start_time
                comparisons_per_second = total_comparisons / processing_time
                
                print(f"  ✓ Success: {processing_time:.2f}s")
                print(f"  Total comparisons: {total_comparisons:,}")
                print(f"  Rate: {comparisons_per_second:.0f} comparisons/second")
            else:
                print(f"  ✗ Failed: {response.status_code}")
                
        except Exception as e:
            print(f"  ✗ Error: {e}")

if __name__ == "__main__":
    print("=== Large Scale Test for get_multi_matchcounts ===")
    print("Testing POST endpoint with 150 molecules and 250 SMARTS patterns...")
    
    success = test_multi_matchcounts_large_scale()
    
    if success:
        print("\n=== Performance Comparison ===")
        test_performance_comparison()
    else:
        print("Large scale test failed. Check your endpoint and try again.")

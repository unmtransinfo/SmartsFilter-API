"""
SMARTS Pattern Validator
A Python script to validate SMARTS patterns using RDKit

Installation:
    pip install rdkit pandas

Usage:
    python smarts_validator.py
"""

from rdkit import Chem
import pandas as pd
from typing import Dict, List
import sys


class SMARTSValidator:
    """Class to validate SMARTS patterns"""
    
    def __init__(self):
        self.results = []
    
    def validate_smarts(self, smarts: str, name: str = "") -> Dict:
        """
        Validate a single SMARTS pattern
        
        Args:
            smarts: SMARTS pattern string
            name: Optional name/identifier for the pattern
            
        Returns:
            Dictionary with validation results
        """
        result = {
            'name': name,
            'smarts': smarts,
            'valid': False,
            'error': None,
            'pattern': None
        }
        
        try:
            pattern = Chem.MolFromSmarts(smarts)
            
            if pattern is None:
                result['error'] = "Invalid SMARTS syntax - pattern could not be parsed"
            else:
                result['valid'] = True
                result['pattern'] = pattern
                result['num_atoms'] = pattern.GetNumAtoms()
                result['num_bonds'] = pattern.GetNumBonds()
                
        except Exception as e:
            result['error'] = f"Exception: {str(e)}"
        
        self.results.append(result)
        return result
    
    def validate_from_file(self, filepath: str, delimiter: str = None) -> pd.DataFrame:
        """
        Validate SMARTS patterns from a file
        
        Args:
            filepath: Path to file with SMARTS patterns
            delimiter: Delimiter between SMARTS and name (auto-detected if None)
            
        Returns:
            DataFrame with validation results
        """
        self.results = []
        
        try:
            with open(filepath, 'r') as f:
                content = f.read()
        except FileNotFoundError:
            print(f"Error: File '{filepath}' not found")
            return pd.DataFrame()
        
        return self.validate_from_text(content, delimiter)
    
    def validate_from_text(self, content: str, delimiter: str = None) -> pd.DataFrame:
        """
        Validate SMARTS patterns from text content
        
        Args:
            content: Text content with SMARTS patterns
            delimiter: Delimiter between SMARTS and name (auto-detected if None)
            
        Returns:
            DataFrame with validation results
        """
        lines = content.strip().split('\n')
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Auto-detect delimiter
            if delimiter is None:
                for delim in ['\t', '  ', ' ', ',']:
                    if delim in line:
                        delimiter = delim
                        break
            
            # Split line into SMARTS and name
            if delimiter and delimiter in line:
                parts = line.split(delimiter, 1)
                smarts = parts[0].strip()
                name = parts[1].strip() if len(parts) > 1 else ""
            else:
                smarts = line
                name = ""
            
            self.validate_smarts(smarts, name)
        
        return self.get_results_df()
    
    def test_smarts_on_molecule(self, smarts: str, smiles: str) -> Dict:
        """
        Test if a SMARTS pattern matches a molecule
        
        Args:
            smarts: SMARTS pattern
            smiles: SMILES string of molecule to test
            
        Returns:
            Dictionary with match results
        """
        result = {
            'smarts': smarts,
            'smiles': smiles,
            'matches': False,
            'num_matches': 0,
            'error': None
        }
        
        try:
            pattern = Chem.MolFromSmarts(smarts)
            mol = Chem.MolFromSmiles(smiles)
            
            if pattern is None:
                result['error'] = "Invalid SMARTS pattern"
            elif mol is None:
                result['error'] = "Invalid SMILES string"
            else:
                matches = mol.GetSubstructMatches(pattern)
                result['matches'] = len(matches) > 0
                result['num_matches'] = len(matches)
                result['match_atoms'] = matches
                
        except Exception as e:
            result['error'] = f"Exception: {str(e)}"
        
        return result
    
    def get_results_df(self) -> pd.DataFrame:
        """Convert results to pandas DataFrame"""
        df_data = []
        for r in self.results:
            df_data.append({
                'Name': r['name'],
                'SMARTS': r['smarts'][:60] + '...' if len(r['smarts']) > 60 else r['smarts'],
                'Valid': '✓' if r['valid'] else '✗',
                'Error': r['error'][:50] + '...' if r['error'] and len(r['error']) > 50 else (r['error'] if r['error'] else ''),
                'Atoms': r.get('num_atoms', ''),
                'Bonds': r.get('num_bonds', '')
            })
        return pd.DataFrame(df_data)
    
    def print_summary(self):
        """Print validation summary"""
        total = len(self.results)
        valid = sum(1 for r in self.results if r['valid'])
        invalid = total - valid
        
        print(f"\n{'='*70}")
        print(f"SMARTS VALIDATION SUMMARY")
        print(f"{'='*70}")
        print(f"Total patterns: {total}")
        print(f"Valid: {valid} ({valid/total*100:.1f}%)" if total > 0 else "Valid: 0")
        print(f"Invalid: {invalid} ({invalid/total*100:.1f}%)" if total > 0 else "Invalid: 0")
        print(f"{'='*70}\n")
    
    def save_results(self, filepath: str):
        """Save validation results to CSV file"""
        df = self.get_results_df()
        df.to_csv(filepath, index=False)
        print(f"Results saved to: {filepath}")


def main():
    """Main function with example usage"""
    
    print("SMARTS Pattern Validator")
    print("=" * 70)
    
    validator = SMARTSValidator()
    
    # Example 1: Validate individual patterns
    print("\n1. Testing individual SMARTS patterns:")
    print("-" * 70)
    
    test_patterns = [
        ("[C,c]", "Carbon (aliphatic or aromatic)"),
        ("[OH]", "Hydroxyl group"),
        ("c1ccccc1", "Benzene ring"),
        ("[N+](=O)[O-]", "Nitro group"),
        ("OO", "Peroxide"),
    ]
    
    for smarts, name in test_patterns:
        result = validator.validate_smarts(smarts, name)
        status = "✓ Valid" if result['valid'] else "✗ Invalid"
        print(f"{status}: {name}")
        print(f"  SMARTS: {smarts}")
        if result['valid']:
            print(f"  Atoms: {result['num_atoms']}, Bonds: {result['num_bonds']}")
        if result['error']:
            print(f"  Error: {result['error']}")
        print()
    
    validator.print_summary()
    
    # Example 2: Validate from your document data
    print("\n2. Validating patterns from sample data:")
    print("-" * 70)
    
    sample_data = """(*).(*) disconnected_structures
[[2H],[3H],[13C],[13c],[14C],[14c],[15N],[15n]] isotopes
[A,a;![C,c,N,n,O,o,F,Si,P,S,s,Cl,Br,I,#1]] AtmCrp
[+1;#6,F,Cl,Br,I,#8,#16;!$([S+][O-])] poschatom
[[N;R0]=,#N] azodiaz
[N=[C;R0]=N] carbdiim
[N(C(=N)N)C(=N)N] biguan
OO peroxides
[O=C1OCC1] 4_membered_lactones
c1ccccc1 benzene"""
    
    validator_doc = SMARTSValidator()
    df = validator_doc.validate_from_text(sample_data)
    print(df.to_string(index=False))
    validator_doc.print_summary()
    
    # Show invalid patterns in detail
    invalid_patterns = [r for r in validator_doc.results if not r['valid']]
    if invalid_patterns:
        print("INVALID PATTERNS - DETAILED ERRORS:")
        print("=" * 70)
        for i, r in enumerate(invalid_patterns, 1):
            print(f"\n{i}. {r['name']}")
            print(f"   SMARTS: {r['smarts']}")
            print(f"   Error: {r['error']}")
    
    # Example 3: Test pattern matching
    print("\n\n3. Testing SMARTS pattern matching:")
    print("-" * 70)
    
    smarts = "c1ccccc1"  # Benzene pattern
    test_molecules = [
        ("c1ccccc1", "Benzene"),
        ("c1ccccc1O", "Phenol"),
        ("CCO", "Ethanol"),
        ("c1ccc(cc1)c2ccccc2", "Biphenyl"),
    ]
    
    for smiles, mol_name in test_molecules:
        result = validator.test_smarts_on_molecule(smarts, smiles)
        if result['matches']:
            print(f"✓ {mol_name}: Matches ({result['num_matches']} occurrence(s))")
        else:
            print(f"✗ {mol_name}: No match")
    
    # Example 4: Save results to file
    print("\n\n4. Saving results:")
    print("-" * 70)
    validator_doc.save_results("smarts_validation_results.csv")
    
    print("\n" + "=" * 70)
    print("Usage Examples:")
    print("=" * 70)
    print("""
# Validate single pattern
v = SMARTSValidator()
result = v.validate_smarts("[OH]", "hydroxyl")
print(result)

# Validate from file
df = v.validate_from_file("patterns.txt")
print(df)
v.print_summary()

# Test pattern against molecule
result = v.test_smarts_on_molecule("c1ccccc1", "c1ccccc1O")
print(result)

# Save results
v.save_results("output.csv")
""")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nValidation interrupted by user")
        sys.exit(0)
    except Exception as e:
        print(f"\n\nError: {e}")
        sys.exit(1)

import sys
import os
import json
import pytest
from flask import Flask

# Add the project root to the path so we can import the app
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from app import create_app

@pytest.fixture
def client():
    app = create_app()
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client

def test_pains_structure(client):
    """Get the baseline structure from PAINS endpoint"""
    smiles = ["c1ccccc1"] # Benzene (should pass)
    response = client.post('/api/v1/smarts_filter/get_filterpains', 
                           json={"SMILES": smiles, "Smile_Names": ["benzene"]})
    assert response.status_code == 200
    data = response.get_json()
    
    # Check top-level keys
    assert "results" in data
    assert "invalid" in data
    assert "all_pains_filters" in data # PAINS specific
    
    # Check result item structure
    result_item = data["results"][0]
    assert "name" in result_item
    assert "smiles" in result_item
    assert "failed" in result_item
    assert "reasons" in result_item
    assert "highlight_atoms" in result_item
    
    return data

def test_new_filters_structure(client):
    """Verify new filters match the baseline structure (mostly)"""
    filters = [
        ("get_filterblake", "Blake"),
        ("get_filterglaxo", "Glaxo"),
        ("get_filteroprea", "Oprea"),
        ("get_filteralarm", "Alarm")
    ]
    
    smiles = ["C=CC(=O)C"] # Methyl vinyl ketone (likely to hit some filters)
    
    for endpoint, name in filters:
        print(f"Testing {name} endpoint: {endpoint}")
        response = client.post(f'/api/v1/smarts_filter/{endpoint}', 
                               json={"SMILES": smiles, "Smile_Names": ["test_mol"]})
        
        assert response.status_code == 200, f"{name} endpoint failed with {response.status_code}"
        data = response.get_json()
        
        # Check top-level keys
        assert "results" in data, f"{name} missing 'results' key"
        assert "invalid" in data, f"{name} missing 'invalid' key"
        assert "all_smarts_filter" in data, f"{name} missing 'all_smarts_filter' key"
        # Verify it's a list
        assert isinstance(data["all_smarts_filter"], list), f"{name} 'all_smarts_filter' is not a list"
        
        # Check result item structure
        assert len(data["results"]) > 0
        result_item = data["results"][0]
        
        required_keys = ["name", "smiles", "failed", "reasons", "highlight_atoms"]
        for key in required_keys:
            assert key in result_item, f"{name} result item missing '{key}'"
            
        print(f"  -> {name} structure verified. Failed: {result_item['failed']}, Reasons: {result_item['reasons']}")

if __name__ == "__main__":
    # Manually run the tests if executed directly
    try:
        app = create_app()
        app.config['TESTING'] = True
        client = app.test_client()
        
        print("--- Verifying PAINS Structure ---")
        pains_data = test_pains_structure(client) # Pass the client object, not the fixture function
        print("PAINS structure OK.")
        
        print("\n--- Verifying New Filters ---")
        test_new_filters_structure(client)
        print("All new filters verified successfully!")
        
    except Exception as e:
        print(f"\nFAILED: {e}")
        import traceback
        traceback.print_exc()

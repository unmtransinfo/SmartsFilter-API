
import os
from rdkit.Chem import FilterCatalog
from utils.smarts_parser import SmartsFile

class CatalogLoader:
    _catalogs = {}

    @classmethod
    def get_catalog(cls, filter_name, file_paths):
        """
        Get or create a FilterCatalog for a specific filter type.
        
        Args:
            filter_name (str): Unique name for the filter (e.g., 'Blake', 'Glaxo')
            file_paths (list): List of absolute paths to .sma files
            
        Returns:
            FilterCatalog.FilterCatalog: The loaded catalog
        """
        if filter_name in cls._catalogs:
            return cls._catalogs[filter_name]

        print(f"Loading {filter_name} filters from {len(file_paths)} files...")
        
        # Create catalog from empty params, then add entries individually
        params = FilterCatalog.FilterCatalogParams()
        catalog = FilterCatalog.FilterCatalog(params)
        
        total_patterns = 0
        for path in file_paths:
            if not os.path.exists(path):
                print(f"Warning: SMARTS file not found: {path}")
                continue
            
            print(f"Parsing {path}...")
            sf = SmartsFile()
            try:
                # Parse using our ported parser
                sf.parse_file(path, strict=False)
                
                # Add valid patterns to catalog
                for smarts_obj in sf.smartses:
                    if smarts_obj.search is not None:
                        matcher = FilterCatalog.SmartsMatcher(
                            smarts_obj.name, 
                            smarts_obj.smarts,  # SMARTS string, not Mol object
                            1, 1
                        )
                        entry = FilterCatalog.FilterCatalogEntry(
                            smarts_obj.name,
                            matcher
                        )
                        catalog.AddEntry(entry)
                        total_patterns += 1
                        
            except Exception as e:
                print(f"Error loading {path}: {e}")

        cls._catalogs[filter_name] = catalog
        print(f"Loaded {total_patterns} patterns for {filter_name}")
        
        return catalog

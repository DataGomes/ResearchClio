"""
Export CLIO results as structured JSON for downstream applications
"""

import json
import os
import sys
from src.output_formatter import OutputFormatter

def export_run_to_json(run_dir, output_file=None, include_abstracts=False):
    """Export a CLIO run to structured JSON format"""
    
    # Load necessary files
    with open(os.path.join(run_dir, 'abstracts.json'), 'r') as f:
        abstracts = json.load(f)
    abstracts_dict = {a['pmid']: a for a in abstracts}
    
    # Load the raw hierarchy saved by Hierarchizer (not the formatted one)
    import glob
    hierarchy_files = glob.glob(os.path.join(run_dir, '*hierarchy*.json'))
    
    # Try to find the raw hierarchy data
    hierarchy = None
    for hf in hierarchy_files:
        with open(hf, 'r') as f:
            data = json.load(f)
            if 'levels' in data and 'base_clusters' in data:
                hierarchy = data
                break
    
    # If not found, reconstruct from available data
    if hierarchy is None:
        # Load cluster names
        with open(os.path.join(run_dir, 'cluster_names.json'), 'r') as f:
            cluster_names = json.load(f)
        
        # Create a basic hierarchy structure
        hierarchy = {
            'base_clusters': cluster_names,
            'levels': []  # Empty levels means single top cluster
        }
    
    # Create formatter
    formatter = OutputFormatter(abstracts_dict, run_dir=run_dir)
    
    # Generate structured JSON
    structured_data = formatter.format_json_output(hierarchy, include_full_abstracts=include_abstracts)
    
    # Save to file
    if output_file is None:
        output_file = os.path.join(run_dir, 'structured_hierarchy.json')
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(structured_data, f, indent=2, ensure_ascii=False)
    
    print(f"Exported structured JSON to: {output_file}")
    
    # Print summary of structure
    print("\nJSON Structure Summary:")
    print(f"- Total abstracts: {structured_data['metadata']['total_abstracts']}")
    print(f"- Hierarchy levels: {structured_data['metadata']['hierarchy_levels']}")
    print(f"- Base clusters: {structured_data['metadata']['base_clusters']}")
    
    if 'domains' in structured_data['hierarchy']:
        print(f"\nDomains ({len(structured_data['hierarchy']['domains'])}):")
        for domain in structured_data['hierarchy']['domains']:
            print(f"  - {domain['name']}: {domain['cluster_count']} clusters, {domain['abstract_count']} abstracts")
    
    return structured_data

def show_json_schema():
    """Show the JSON schema for downstream applications"""
    schema = {
        "metadata": {
            "generated_at": "ISO timestamp",
            "total_abstracts": "integer",
            "hierarchy_levels": "integer",
            "base_clusters": "integer"
        },
        "statistics": {
            "cluster_sizes": {
                "min": "integer",
                "max": "integer", 
                "mean": "float"
            }
        },
        "hierarchy": {
            "top_level": {
                "id": "string",
                "name": "string",
                "description": "string",
                "total_abstracts": "integer",
                "children_count": "integer"
            },
            "domains": [
                {
                    "name": "string (e.g., 'Medical Imaging & Radiology')",
                    "cluster_count": "integer",
                    "abstract_count": "integer",
                    "clusters": [
                        {
                            "id": "integer",
                            "name": "string",
                            "description": "string",
                            "abstract_count": "integer",
                            "abstract_ids": ["list of PMIDs"],
                            "abstracts": ["optional: full abstract data"]
                        }
                    ]
                }
            ],
            "base_clusters": [
                {
                    "id": "integer",
                    "name": "string",
                    "description": "string",
                    "abstract_count": "integer"
                }
            ]
        }
    }
    
    print("JSON Schema for CLIO Output:")
    print(json.dumps(schema, indent=2))

if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "--schema":
            show_json_schema()
        else:
            run_dir = sys.argv[1]
            include_abstracts = "--include-abstracts" in sys.argv
            export_run_to_json(run_dir, include_abstracts=include_abstracts)
    else:
        # Default to most recent run
        run_dir = "/Users/andregomes/Library/CloudStorage/GoogleDrive-andmagal@gmail.com/My Drive/DataGomes/clio/output/run_20250524_211451"
        export_run_to_json(run_dir)
"""
Create complete hierarchy view combining all available data
"""

import json
import os
import sys
import pickle

def create_complete_hierarchy(run_dir):
    """Combine all hierarchy data into one comprehensive view"""
    
    # Load all available data
    with open(os.path.join(run_dir, 'abstracts.json'), 'r') as f:
        abstracts = json.load(f)
    abstracts_dict = {a['pmid']: a for a in abstracts}
    
    with open(os.path.join(run_dir, 'cluster_names.json'), 'r') as f:
        cluster_names = json.load(f)
    
    with open(os.path.join(run_dir, 'clustering.pkl'), 'rb') as f:
        clustering = pickle.load(f)
    
    # Count abstracts per cluster
    from collections import Counter
    cluster_counts = Counter(clustering['labels'])
    
    # Check for hierarchy.json (the original output)
    hierarchy_file = os.path.join(run_dir, 'hierarchy.json')
    if os.path.exists(hierarchy_file):
        with open(hierarchy_file, 'r') as f:
            original_hierarchy = json.load(f)
    else:
        original_hierarchy = None
    
    # Build complete output
    lines = []
    lines.append("COMPLETE CLIO HIERARCHICAL CLUSTERING RESULTS")
    lines.append("=" * 80)
    lines.append("")
    lines.append(f"Total abstracts processed: {len(abstracts)}")
    lines.append(f"Base-level clusters created: {len(cluster_names)}")
    lines.append("")
    
    # Show the top-level cluster from original hierarchy
    if original_hierarchy and 'hierarchy' in original_hierarchy:
        lines.append("TOP-LEVEL CLUSTER (from LLM clustering):")
        lines.append("=" * 80)
        for cluster in original_hierarchy['hierarchy']:
            lines.append(f"▪ {cluster['name']}")
            lines.append(f"  {cluster['description']}")
            lines.append(f"  Total abstracts: {original_hierarchy['total_abstracts']}")
            lines.append("")
    
    # Show cluster statistics
    lines.append("CLUSTER STATISTICS:")
    lines.append("-" * 40)
    lines.append(f"Smallest cluster: {min(cluster_counts.values())} abstracts")
    lines.append(f"Largest cluster: {max(cluster_counts.values())} abstracts")
    lines.append(f"Average size: {sum(cluster_counts.values()) / len(cluster_counts):.1f} abstracts")
    lines.append("")
    
    # Show all base clusters
    lines.append("ALL BASE CLUSTERS (sorted by size):")
    lines.append("=" * 80)
    lines.append("")
    
    # Sort clusters by size
    sorted_clusters = sorted(cluster_names.items(), 
                           key=lambda x: cluster_counts.get(int(x[0]), 0), 
                           reverse=True)
    
    for rank, (cid, info) in enumerate(sorted_clusters, 1):
        count = cluster_counts.get(int(cid), 0)
        lines.append(f"{rank}. {info['name']} ({count} abstracts)")
        lines.append(f"   Cluster ID: {cid}")
        lines.append(f"   {info['description']}")
        
        # Get sample abstracts for this cluster
        cluster_abstracts = []
        for i, label in enumerate(clustering['labels']):
            if label == int(cid):
                pmid = clustering['pmids'][i]
                if pmid in abstracts_dict:
                    cluster_abstracts.append(abstracts_dict[pmid])
        
        # Show first 3 abstracts
        if cluster_abstracts:
            lines.append("   Sample abstracts:")
            for j, abstract in enumerate(cluster_abstracts[:3], 1):
                lines.append(f"     {j}) [{abstract['pmid']}] {abstract['title'][:70]}...")
        
        lines.append("")
    
    # Create complete JSON structure
    complete_json = {
        "metadata": {
            "total_abstracts": len(abstracts),
            "base_clusters": len(cluster_names),
            "hierarchy_levels": 1  # Since everything converged to 1 top cluster
        },
        "top_level_cluster": {
            "name": "Artificial Intelligence in Healthcare and Life Sciences",
            "description": original_hierarchy['hierarchy'][0]['description'] if original_hierarchy else "All abstracts converged into a single cluster",
            "total_abstracts": len(abstracts)
        },
        "base_clusters": []
    }
    
    # Add base clusters to JSON
    for cid, info in sorted_clusters:
        cluster_data = {
            "id": int(cid),
            "name": info['name'],
            "description": info['description'],
            "abstract_count": cluster_counts.get(int(cid), 0),
            "pmids": []
        }
        
        # Get PMIDs for this cluster
        for i, label in enumerate(clustering['labels']):
            if label == int(cid):
                cluster_data['pmids'].append(clustering['pmids'][i])
        
        complete_json['base_clusters'].append(cluster_data)
    
    # Save outputs
    text_output = '\n'.join(lines)
    
    # Save text version
    text_file = os.path.join(run_dir, 'complete_hierarchy.txt')
    with open(text_file, 'w', encoding='utf-8') as f:
        f.write(text_output)
    print(f"Saved complete hierarchy text to: {text_file}")
    
    # Save JSON version
    json_file = os.path.join(run_dir, 'complete_hierarchy.json')
    with open(json_file, 'w', encoding='utf-8') as f:
        json.dump(complete_json, f, indent=2, ensure_ascii=False)
    print(f"Saved complete hierarchy JSON to: {json_file}")
    
    return text_output, complete_json

if __name__ == "__main__":
    run_dir = sys.argv[1] if len(sys.argv) > 1 else "/Users/andregomes/Library/CloudStorage/GoogleDrive-andmagal@gmail.com/My Drive/DataGomes/clio/output/run_20250524_211451"
    create_complete_hierarchy(run_dir)
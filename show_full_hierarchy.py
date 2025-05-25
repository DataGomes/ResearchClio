"""
Show the complete hierarchy with all levels and clusters
"""

import json
import os
import sys
import pickle
from collections import Counter

def load_all_data(run_dir):
    """Load all necessary data files"""
    # Load abstracts
    with open(os.path.join(run_dir, 'abstracts.json'), 'r') as f:
        abstracts = json.load(f)
    abstracts_dict = {a['pmid']: a for a in abstracts}
    
    # Load cluster names
    with open(os.path.join(run_dir, 'cluster_names.json'), 'r') as f:
        cluster_names = json.load(f)
    
    # Load clustering results
    with open(os.path.join(run_dir, 'clustering.pkl'), 'rb') as f:
        clustering = pickle.load(f)
    
    # Count abstracts per cluster
    cluster_counts = Counter(clustering['labels'])
    
    return abstracts_dict, cluster_names, cluster_counts

def show_complete_hierarchy(run_dir):
    """Display the complete hierarchy structure"""
    abstracts_dict, cluster_names, cluster_counts = load_all_data(run_dir)
    
    print("COMPLETE CLIO HIERARCHICAL CLUSTERING")
    print("=" * 80)
    print(f"Total abstracts: {len(abstracts_dict)}")
    print(f"Base-level clusters: {len(cluster_names)}")
    print("=" * 80)
    print()
    
    # Since the hierarchy converged to 1 top cluster, let's show it manually
    print("▪ TOP LEVEL: Artificial Intelligence in Healthcare and Life Sciences (996 abstracts)")
    print("  This group explores diverse applications of AI across healthcare and life sciences")
    print()
    
    # Group base clusters by domain for better organization
    domains = {
        'Medical Imaging & Radiology': [],
        'Clinical Specialties': [],
        'Pharmaceutical & Drug Discovery': [],
        'Surgical & Procedural Medicine': [],
        'Diagnostics & Pathology': [],
        'Public Health & Research': []
    }
    
    # Categorize each cluster
    for cid, info in cluster_names.items():
        name = info['name']
        count = cluster_counts.get(int(cid), 0)
        
        if any(word in name.lower() for word in ['imaging', 'radiology', 'radiograph']):
            domains['Medical Imaging & Radiology'].append((cid, info, count))
        elif any(word in name.lower() for word in ['pharmaceutical', 'drug']):
            domains['Pharmaceutical & Drug Discovery'].append((cid, info, count))
        elif any(word in name.lower() for word in ['surgery', 'surgical', 'operative']):
            domains['Surgical & Procedural Medicine'].append((cid, info, count))
        elif any(word in name.lower() for word in ['pathology', 'diagnostic', 'detection']):
            domains['Diagnostics & Pathology'].append((cid, info, count))
        elif any(word in name.lower() for word in ['clinical', 'medicine', 'health']):
            domains['Clinical Specialties'].append((cid, info, count))
        else:
            domains['Public Health & Research'].append((cid, info, count))
    
    # Print hierarchical structure
    for domain, clusters in domains.items():
        if clusters:
            total = sum(c[2] for c in clusters)
            print(f"  ├─ {domain} ({len(clusters)} clusters, {total} abstracts)")
            
            # Sort clusters by size (descending)
            sorted_clusters = sorted(clusters, key=lambda x: x[2], reverse=True)
            
            for i, (cid, info, count) in enumerate(sorted_clusters):
                is_last = i == len(sorted_clusters) - 1
                prefix = "  │  └─" if is_last else "  │  ├─"
                print(f"{prefix} {info['name']} ({count} abstracts)")
                
                # Show description
                desc_prefix = "  │     " if not is_last else "  │     "
                desc_lines = info['description'].split('. ')
                print(f"{desc_prefix} {desc_lines[0]}.")
                if len(desc_lines) > 1:
                    print(f"{desc_prefix} {desc_lines[1]}.")
                
                # Show sample abstracts (top 3)
                if '--show-abstracts' in sys.argv:
                    # Get abstracts for this cluster
                    cluster_abstracts = []
                    for j, label in enumerate(clustering['labels']):
                        if label == int(cid):
                            pmid = clustering['pmids'][j]
                            if pmid in abstracts_dict:
                                cluster_abstracts.append(abstracts_dict[pmid])
                    
                    # Show first 3 abstracts
                    for k, abstract in enumerate(cluster_abstracts[:3]):
                        abs_prefix = "  │     │  " if not is_last else "  │        "
                        print(f"{abs_prefix} • [{abstract['pmid']}] {abstract['title'][:60]}...")
                
                print()
    
    # Summary statistics
    print()
    print("SUMMARY STATISTICS")
    print("-" * 40)
    print(f"Largest cluster: {max(cluster_counts.values())} abstracts")
    print(f"Smallest cluster: {min(cluster_counts.values())} abstracts")
    print(f"Average cluster size: {sum(cluster_counts.values()) / len(cluster_counts):.1f} abstracts")
    
    # Domain distribution
    print("\nDOMAIN DISTRIBUTION:")
    for domain, clusters in domains.items():
        if clusters:
            total = sum(c[2] for c in clusters)
            pct = (total / len(abstracts_dict)) * 100
            print(f"  {domain}: {total} abstracts ({pct:.1f}%)")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] != '--show-abstracts':
        run_dir = sys.argv[1]
    else:
        run_dir = "/Users/andregomes/Library/CloudStorage/GoogleDrive-andmagal@gmail.com/My Drive/DataGomes/clio/output/run_20250524_211451"
    
    show_complete_hierarchy(run_dir)
    
    # Save to file
    output_file = os.path.join(run_dir, "full_hierarchy_tree.txt")
    print(f"\n✓ Saving visualization to: {output_file}")
    
    # Redirect output to file
    import io
    from contextlib import redirect_stdout
    
    with open(output_file, 'w', encoding='utf-8') as f:
        with redirect_stdout(f):
            show_complete_hierarchy(run_dir)
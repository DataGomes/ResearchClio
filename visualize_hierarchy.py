"""
Visualize the complete hierarchy from a CLIO run
"""

import json
import sys
import os

def load_run_data(run_dir):
    """Load all data from a run directory"""
    # Load cluster names
    with open(os.path.join(run_dir, 'cluster_names.json'), 'r') as f:
        cluster_names = json.load(f)
    
    # Load hierarchy
    hierarchy_path = os.path.join(run_dir, 'hierarchy.json')
    if os.path.exists(hierarchy_path):
        with open(hierarchy_path, 'r') as f:
            hierarchy_data = json.load(f)
    else:
        hierarchy_data = None
    
    # Load abstracts to get counts
    with open(os.path.join(run_dir, 'abstracts.json'), 'r') as f:
        abstracts = json.load(f)
    
    # Load clustering results
    import pickle
    with open(os.path.join(run_dir, 'clustering.pkl'), 'rb') as f:
        clustering = pickle.load(f)
    
    return cluster_names, hierarchy_data, abstracts, clustering

def print_hierarchy_recursive(node, level=0, prefix="", cluster_counts=None, abstracts_dict=None):
    """Recursively print hierarchy nodes"""
    indent = "  " * level
    
    # Print current node
    if 'name' in node:
        abstract_count = ""
        if 'children' in node and node['children']:
            # Count abstracts in children
            if isinstance(node['children'][0], dict):
                # Children are clusters
                child_count = len(node['children'])
                abstract_count = f" ({child_count} sub-clusters)"
            else:
                # Children are abstract IDs
                abstract_count = f" ({len(node['children'])} abstracts)"
        
        print(f"{indent}{prefix}{node['name']}{abstract_count}")
        if 'description' in node:
            print(f"{indent}  {node['description'][:120]}...")
    
    # Print children
    if 'children' in node and node['children']:
        for i, child in enumerate(node['children']):
            if isinstance(child, dict):
                # Child is another cluster
                child_prefix = "├─ " if i < len(node['children']) - 1 else "└─ "
                print_hierarchy_recursive(child, level + 1, child_prefix, cluster_counts, abstracts_dict)
            else:
                # Child is an abstract ID (leaf node)
                if abstracts_dict and str(child) in abstracts_dict:
                    abstract = abstracts_dict[str(child)]
                    child_prefix = "├─ " if i < len(node['children']) - 1 else "└─ "
                    print(f"{indent}  {child_prefix}[{abstract['pmid']}] {abstract['title'][:80]}...")

def print_full_hierarchy(run_dir):
    """Print the complete hierarchy from the actual hierarchy.json file"""
    # Load data
    cluster_names, hierarchy_data, abstracts, clustering = load_run_data(run_dir)
    
    # Convert abstracts list to dict for easy lookup
    abstracts_dict = {a['pmid']: a for a in abstracts}
    
    # Count abstracts per cluster
    from collections import Counter
    cluster_counts = Counter(clustering['labels'])
    
    print("COMPLETE HIERARCHICAL CLUSTERING RESULTS")
    print("=" * 80)
    print(f"Total abstracts: {len(abstracts)}")
    print(f"Base-level clusters: {len(cluster_names)}")
    
    # Load the raw hierarchy file to get the actual structure
    import pickle
    hierarchy_file = os.path.join(run_dir, 'hierarchy.json')
    with open(hierarchy_file, 'r') as f:
        raw_hierarchy = json.load(f)
    
    # Check if we have the internal structure
    if 'levels' in raw_hierarchy and 'base_clusters' in raw_hierarchy:
        print(f"Hierarchy levels: {len(raw_hierarchy['levels']) + 1}")
        print("=" * 80)
        print()
        
        # Build complete tree from raw hierarchy
        if raw_hierarchy['levels']:
            # Start from top level
            top_level = raw_hierarchy['levels'][-1]['clusters']
            
            for cluster_id, cluster in top_level.items():
                # Create full tree structure
                tree = build_tree_from_hierarchy(cluster, raw_hierarchy, cluster_counts, abstracts_dict)
                print_hierarchy_recursive(tree, prefix="▪ ", cluster_counts=cluster_counts, abstracts_dict=abstracts_dict)
                print()
        else:
            # No hierarchy levels, just show base clusters grouped by theme
            print("\nNOTE: Hierarchical building resulted in a single top-level cluster.")
            print("Showing thematic grouping of base clusters instead:\n")
            print_themed_clusters(cluster_names, cluster_counts)
    else:
        # Fallback to themed view
        print("=" * 80)
        print()
        print_themed_clusters(cluster_names, cluster_counts)

def build_tree_from_hierarchy(cluster, raw_hierarchy, cluster_counts, abstracts_dict):
    """Build a complete tree structure from the hierarchy data"""
    tree = {
        'name': cluster['name'],
        'description': cluster.get('description', ''),
        'children': []
    }
    
    # Get children IDs
    children_ids = cluster.get('children', [])
    
    if children_ids:
        # Find which level these children are in
        child_clusters = None
        
        # Check each level from bottom to top
        for level_info in reversed(raw_hierarchy['levels'][:-1]):
            if any(str(cid) in level_info['clusters'] for cid in children_ids):
                child_clusters = level_info['clusters']
                break
        
        # If not in levels, check base clusters
        if child_clusters is None:
            child_clusters = raw_hierarchy['base_clusters']
        
        # Build children
        for child_id in children_ids:
            child_id_str = str(child_id)
            if child_id_str in child_clusters:
                # Recursively build child cluster
                child_tree = build_tree_from_hierarchy(
                    child_clusters[child_id_str], 
                    raw_hierarchy, 
                    cluster_counts, 
                    abstracts_dict
                )
                tree['children'].append(child_tree)
            elif child_id_str in raw_hierarchy['base_clusters']:
                # This is a base cluster
                base_cluster = raw_hierarchy['base_clusters'][child_id_str]
                # Add actual abstract count
                base_tree = {
                    'name': base_cluster['name'],
                    'description': base_cluster.get('description', ''),
                    'abstract_count': cluster_counts.get(int(child_id), 0)
                }
                tree['children'].append(base_tree)
    
    return tree

def print_themed_clusters(cluster_names, cluster_counts):
    """Print clusters organized by theme"""
    # Group clusters by theme (simple grouping based on keywords)
    themes = {
        'Imaging & Radiology': [],
        'Clinical Medicine': [],
        'Surgery & Procedures': [],
        'Oncology': [],
        'Neurology': [],
        'Cardiology': [],
        'Pharmaceutical': [],
        'Other Specialties': []
    }
    
    # Categorize clusters
    for cid, info in cluster_names.items():
        name = info['name'].lower()
        categorized = False
        
        if any(word in name for word in ['imaging', 'radiology', 'radiograph', 'image']):
            themes['Imaging & Radiology'].append((cid, info))
            categorized = True
        elif any(word in name for word in ['oncology', 'cancer', 'tumor']):
            themes['Oncology'].append((cid, info))
            categorized = True
        elif any(word in name for word in ['surgery', 'surgical', 'operative']):
            themes['Surgery & Procedures'].append((cid, info))
            categorized = True
        elif any(word in name for word in ['neuro', 'brain', 'cerebr']):
            themes['Neurology'].append((cid, info))
            categorized = True
        elif any(word in name for word in ['cardio', 'heart', 'ecg', 'electrocardiogram']):
            themes['Cardiology'].append((cid, info))
            categorized = True
        elif any(word in name for word in ['pharmaceutical', 'drug']):
            themes['Pharmaceutical'].append((cid, info))
            categorized = True
        elif any(word in name for word in ['clinical', 'medicine', 'health']):
            themes['Clinical Medicine'].append((cid, info))
            categorized = True
        
        if not categorized:
            themes['Other Specialties'].append((cid, info))
    
    # Print organized hierarchy
    for theme, clusters in themes.items():
        if clusters:
            print(f"▪ {theme.upper()} ({len(clusters)} clusters)")
            print("-" * 70)
            
            for cid, info in sorted(clusters, key=lambda x: x[1]['name']):
                count = cluster_counts[int(cid)]
                print(f"  └─ {info['name']} ({count} abstracts)")
                print(f"     {info['description'][:100]}...")
                print()
            
            print()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        run_dir = sys.argv[1]
    else:
        # Default to most recent run
        run_dir = "/Users/andregomes/Library/CloudStorage/GoogleDrive-andmagal@gmail.com/My Drive/DataGomes/clio/output/run_20250524_211451"
    
    # Capture output to both console and file
    import io
    from contextlib import redirect_stdout
    
    # Create a string buffer to capture output
    output_buffer = io.StringIO()
    
    # Print to both console and buffer
    class Tee:
        def __init__(self, *files):
            self.files = files
        def write(self, obj):
            for f in self.files:
                f.write(obj)
                f.flush()
        def flush(self):
            for f in self.files:
                f.flush()
    
    # Redirect stdout to both console and buffer
    original_stdout = sys.stdout
    sys.stdout = Tee(sys.stdout, output_buffer)
    
    # Run the visualization
    print_full_hierarchy(run_dir)
    
    # Restore stdout
    sys.stdout = original_stdout
    
    # Save to file
    output_file = os.path.join(run_dir, "complete_hierarchy_visualization.txt")
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(output_buffer.getvalue())
    
    print(f"\n✓ Visualization saved to: {output_file}")
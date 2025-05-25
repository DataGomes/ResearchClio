"""
Convert structured JSON hierarchy to formatted text
"""

import json
import sys
import os

def json_to_text(json_file, output_file=None):
    """Convert structured JSON to formatted text"""
    
    # Load JSON
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # Build text output
    lines = []
    
    # Header
    lines.append("CLIO HIERARCHICAL CLUSTERING RESULTS")
    lines.append("=" * 80)
    lines.append("")
    
    # Metadata
    meta = data['metadata']
    lines.append(f"Generated: {meta['generated_at']}")
    lines.append(f"Total abstracts: {meta['total_abstracts']}")
    lines.append(f"Base clusters: {meta['base_clusters']}")
    lines.append(f"Hierarchy levels: {meta['hierarchy_levels']}")
    lines.append("")
    
    # Statistics
    stats = data['statistics']['cluster_sizes']
    lines.append("CLUSTER SIZE STATISTICS")
    lines.append("-" * 40)
    lines.append(f"Minimum: {stats['min']} abstracts")
    lines.append(f"Maximum: {stats['max']} abstracts")
    lines.append(f"Average: {stats['mean']:.1f} abstracts")
    lines.append("")
    
    # Hierarchy
    hierarchy = data['hierarchy']
    
    # Show top-level cluster if it exists
    if 'levels' in hierarchy and hierarchy['levels']:
        # We have the full hierarchy structure
        lines.append("HIERARCHY STRUCTURE")
        lines.append("=" * 80)
        
        # Show top level first
        if len(hierarchy['levels']) > 0:
            top_level = hierarchy['levels'][0]  # Assuming top level is first
            lines.append("\nTOP LEVEL:")
            lines.append("-" * 40)
            for cluster in top_level.get('clusters', []):
                lines.append(f"▪ {cluster['name']}")
                lines.append(f"  {cluster['description']}")
                if cluster.get('children'):
                    lines.append(f"  Children: {cluster['children']}")
                lines.append("")
    
    # Also check if there's a single top cluster (from the simple hierarchy format)
    if not hierarchy.get('levels') and not hierarchy.get('base_clusters'):
        # This might be the original format from hierarchy.json
        lines.append("\nTOP LEVEL CLUSTER:")
        lines.append("=" * 80)
        lines.append(f"▪ Artificial Intelligence in Healthcare and Life Sciences")
        lines.append(f"  All {meta['total_abstracts']} abstracts converged into this single top-level cluster")
        lines.append(f"  encompassing diverse AI applications across healthcare and life sciences")
        lines.append("")
    
    # Base clusters
    lines.append("\nBASE CLUSTERS (sorted by size)")
    lines.append("=" * 80)
    lines.append("")
    
    for i, cluster in enumerate(hierarchy['base_clusters'], 1):
        lines.append(f"{i}. {cluster['name']} ({cluster['abstract_count']} abstracts)")
        lines.append(f"   ID: {cluster['id']}")
        lines.append(f"   {cluster['description']}")
        
        # Show sample abstract IDs if available
        if 'abstract_ids' in cluster:
            sample_ids = cluster['abstract_ids'][:5]  # First 5
            if len(cluster['abstract_ids']) > 5:
                lines.append(f"   Sample PMIDs: {', '.join(sample_ids)}... (and {len(cluster['abstract_ids'])-5} more)")
            else:
                lines.append(f"   PMIDs: {', '.join(sample_ids)}")
        
        # Show abstracts if included
        if 'abstracts' in cluster and cluster['abstracts']:
            lines.append("   Sample abstracts:")
            for j, abstract in enumerate(cluster['abstracts'][:3], 1):
                lines.append(f"     {j}) [{abstract['pmid']}] {abstract['title'][:80]}...")
        
        lines.append("")
    
    # Summary
    lines.append("\nSUMMARY")
    lines.append("=" * 80)
    total_abstracts = sum(c['abstract_count'] for c in hierarchy['base_clusters'])
    lines.append(f"Total abstracts across all clusters: {total_abstracts}")
    lines.append(f"Number of base clusters: {len(hierarchy['base_clusters'])}")
    
    # Distribution
    lines.append("\nCLUSTER SIZE DISTRIBUTION:")
    size_ranges = [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100)]
    for min_size, max_size in size_ranges:
        count = sum(1 for c in hierarchy['base_clusters'] 
                   if min_size < c['abstract_count'] <= max_size)
        if count > 0:
            lines.append(f"  {min_size:3d}-{max_size:3d} abstracts: {count} clusters")
    
    # Join all lines
    text_output = '\n'.join(lines)
    
    # Save to file
    if output_file is None:
        output_file = json_file.replace('.json', '_formatted.txt')
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(text_output)
    
    print(f"Saved formatted text to: {output_file}")
    return text_output

if __name__ == "__main__":
    if len(sys.argv) > 1:
        json_file = sys.argv[1]
    else:
        json_file = "/Users/andregomes/Library/CloudStorage/GoogleDrive-andmagal@gmail.com/My Drive/DataGomes/clio/output/run_20250524_211451/structured_hierarchy.json"
    
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    json_to_text(json_file, output_file)
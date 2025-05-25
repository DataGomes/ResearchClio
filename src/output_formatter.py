"""
Output Formatter Module
Formats hierarchical clustering results for display and export
"""

import os
import json
from typing import Dict, List, Optional
import pandas as pd
from datetime import datetime
import numpy as np

def convert_to_serializable(obj):
    """Convert numpy types to Python native types for JSON serialization"""
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(i) for i in obj]
    else:
        return obj

class OutputFormatter:
    def __init__(self, abstracts_data: Dict[str, Dict], run_dir: str = None):
        """Initialize with abstracts data for reference"""
        self.abstracts = abstracts_data
        self.run_dir = run_dir
    
    def format_hierarchy_text(self, hierarchy: Dict, 
                            indent: str = "  ",
                            show_abstracts: bool = False,
                            max_abstracts_per_cluster: int = 3) -> str:
        """Format hierarchy as indented text"""
        output = []
        output.append("HIERARCHICAL CLUSTERING OF AI ABSTRACTS")
        output.append("=" * 80)
        
        # Load clustering data to get abstract counts
        import pickle
        import os
        from collections import Counter
        
        run_dir = os.path.dirname(list(self.abstracts.values())[0].get('_source_file', ''))
        if run_dir and os.path.exists(os.path.join(run_dir, 'clustering.pkl')):
            with open(os.path.join(run_dir, 'clustering.pkl'), 'rb') as f:
                clustering = pickle.load(f)
            cluster_counts = Counter(clustering['labels'])
        else:
            cluster_counts = {}
        
        # Get statistics
        total_abstracts = len(self.abstracts)
        if 'base_clusters' in hierarchy:
            base_clusters = hierarchy['base_clusters']
            num_base_clusters = len(base_clusters)
        else:
            base_clusters = {}
            num_base_clusters = 0
        
        output.append(f"Total abstracts: {total_abstracts}")
        output.append(f"Base-level clusters: {num_base_clusters}")
        
        # Determine hierarchy depth
        if 'levels' in hierarchy:
            hierarchy_depth = len(hierarchy['levels']) + 1
        else:
            hierarchy_depth = 1
            
        output.append(f"Hierarchy levels: {hierarchy_depth}")
        output.append("=" * 80)
        output.append("")
        
        # If we have a proper hierarchy structure
        if 'levels' in hierarchy and hierarchy['levels']:
            # Get top level
            top_level = hierarchy['levels'][-1]['clusters']
            
            # For single top cluster, show complete hierarchy
            if len(top_level) == 1:
                cluster_id, cluster = list(top_level.items())[0]
                output.append(f"▪ TOP LEVEL: {cluster['name']} ({total_abstracts} abstracts)")
                output.append(f"  {cluster['description']}")
                output.append("")
                
                # Show all base clusters sorted by size
                sorted_base = sorted(base_clusters.items(), 
                                   key=lambda x: cluster_counts.get(int(x[0]), 0), 
                                   reverse=True)
                
                output.append(f"  BASE CLUSTERS ({len(sorted_base)} clusters):")
                output.append("")
                
                for cid, info in sorted_base:
                    count = cluster_counts.get(int(cid), 0)
                    output.append(f"  └─ {info['name']} ({count} abstracts)")
                    output.append(f"     {info['description']}")
                    output.append("")
            else:
                # Multiple top clusters - use original recursive format
                output.append(f"TOP LEVEL ({len(top_level)} clusters):")
                output.append("")
                
                for cluster_id, cluster in top_level.items():
                    self._format_cluster_recursive(
                        cluster, hierarchy, 0, output, indent,
                        show_abstracts, max_abstracts_per_cluster
                    )
                    output.append("")
        else:
            # No hierarchy levels, just show base clusters
            output.append("Note: Hierarchical clustering converged to a single group.")
            output.append("Showing all base clusters:")
            output.append("")
            
            # Sort by cluster size
            sorted_base = sorted(base_clusters.items(), 
                               key=lambda x: cluster_counts.get(int(x[0]), 0), 
                               reverse=True)
            
            for cid, info in sorted_base:
                count = cluster_counts.get(int(cid), 0)
                output.append(f"└─ {info['name']} ({count} abstracts)")
                output.append(f"   {info['description']}")
                output.append("")
        
        # Add summary statistics
        output.append("")
        output.append("SUMMARY STATISTICS")
        output.append("-" * 40)
        if cluster_counts:
            output.append(f"Largest cluster: {max(cluster_counts.values())} abstracts")
            output.append(f"Smallest cluster: {min(cluster_counts.values())} abstracts")
            output.append(f"Average cluster size: {sum(cluster_counts.values()) / len(cluster_counts):.1f} abstracts")
        
        return "\n".join(output)
    
    
    def _format_cluster_recursive(self, cluster: Dict, hierarchy: Dict,
                                level: int, output: List[str], indent: str,
                                show_abstracts: bool, max_abstracts: int):
        """Recursively format a cluster and its children"""
        # Format current cluster
        prefix = indent * level
        output.append(f"{prefix}▪ {cluster['name']}")
        output.append(f"{prefix}  {cluster['description']}")
        
        # Check if this cluster has children in lower levels
        children_ids = cluster.get('children', [])
        
        if children_ids:
            # Find which level these children are in
            child_clusters = None
            for level_info in reversed(hierarchy['levels'][:-1]):
                if any(cid in level_info['clusters'] for cid in map(str, children_ids)):
                    child_clusters = level_info['clusters']
                    break
            
            # If not in levels, check base clusters
            if child_clusters is None:
                child_clusters = hierarchy['base_clusters']
            
            # Format children
            for child_id in children_ids:
                child_id_str = str(child_id)
                if child_id_str in child_clusters:
                    child = child_clusters[child_id_str]
                    self._format_cluster_recursive(
                        child, hierarchy, level + 1, output, indent,
                        show_abstracts, max_abstracts
                    )
                elif show_abstracts and child_id_str in self.abstracts:
                    # This is a leaf node - show abstract
                    abstract = self.abstracts[child_id_str]
                    output.append(f"{indent * (level + 1)}• {abstract['title'][:80]}...")
                    output.append(f"{indent * (level + 1)}  PMID: {abstract['pmid']}")
    
    def format_json_output(self, hierarchy: Dict, 
                         include_full_abstracts: bool = False) -> Dict:
        """Format hierarchy as structured JSON with complete hierarchy information"""
        # Load clustering data if available
        import pickle
        import os
        from collections import Counter
        
        cluster_counts = {}
        cluster_to_abstracts = {}
        
        if self.run_dir and os.path.exists(os.path.join(self.run_dir, 'clustering.pkl')):
            with open(os.path.join(self.run_dir, 'clustering.pkl'), 'rb') as f:
                clustering = pickle.load(f)
            cluster_counts = Counter(clustering['labels'])
            
            # Map clusters to abstract IDs
            for pmid, label in zip(clustering['pmids'], clustering['labels']):
                if label not in cluster_to_abstracts:
                    cluster_to_abstracts[label] = []
                cluster_to_abstracts[label].append(pmid)
        
        # Build complete structure
        result = {
            "metadata": {
                "generated_at": datetime.now().isoformat(),
                "total_abstracts": len(self.abstracts),
                "hierarchy_levels": len(hierarchy.get('levels', [])) + 1 if 'levels' in hierarchy else 1,
                "base_clusters": len(hierarchy.get('base_clusters', {}))
            },
            "statistics": {
                "cluster_sizes": {
                    "min": min(cluster_counts.values()) if cluster_counts else 0,
                    "max": max(cluster_counts.values()) if cluster_counts else 0,
                    "mean": sum(cluster_counts.values()) / len(cluster_counts) if cluster_counts else 0
                }
            },
            "hierarchy": self._build_complete_json_hierarchy(
                hierarchy, cluster_counts, cluster_to_abstracts, include_full_abstracts
            )
        }
        
        return result
    
    def _build_complete_json_hierarchy(self, hierarchy: Dict, cluster_counts: Dict,
                                     cluster_to_abstracts: Dict,
                                     include_full_abstracts: bool) -> Dict:
        """Build complete JSON hierarchy with all levels and metadata"""
        json_hierarchy = {
            "levels": [],
            "base_clusters": []
        }
        
        # Get base clusters
        base_clusters = hierarchy.get('base_clusters', {})
        
        # Build the complete hierarchy from top to bottom
        if 'levels' in hierarchy and hierarchy['levels']:
            # We have multiple levels - build them all
            for level_idx, level_data in enumerate(reversed(hierarchy['levels'])):
                level_clusters = []
                for cluster_id, cluster in level_data['clusters'].items():
                    cluster_entry = {
                        "id": cluster_id,
                        "name": cluster['name'],
                        "description": cluster['description'],
                        "children": cluster.get('children', [])
                    }
                    level_clusters.append(cluster_entry)
                
                json_hierarchy["levels"].append({
                    "level": len(hierarchy['levels']) - level_idx,
                    "clusters": level_clusters
                })
        
        # Add base clusters with their data
        for cid, info in base_clusters.items():
            cluster_data = {
                "id": int(cid),
                "name": info['name'],
                "description": info['description'],
                "abstract_count": cluster_counts.get(int(cid), 0)
            }
            
            # Add abstract IDs or full abstracts
            if int(cid) in cluster_to_abstracts:
                abstract_ids = cluster_to_abstracts[int(cid)]
                if include_full_abstracts:
                    cluster_data["abstracts"] = [
                        {
                            "pmid": pmid,
                            "title": self.abstracts[pmid]['title'],
                            "abstract": self.abstracts[pmid]['abstract'],
                            "journal": self.abstracts[pmid].get('journal', ''),
                            "authors": self.abstracts[pmid].get('authors', [])
                        }
                        for pmid in abstract_ids if pmid in self.abstracts
                    ]
                else:
                    cluster_data["abstract_ids"] = abstract_ids
            
            json_hierarchy["base_clusters"].append(cluster_data)
        
        # Sort base clusters by size for easier navigation
        json_hierarchy["base_clusters"].sort(key=lambda x: x['abstract_count'], reverse=True)
        
        return json_hierarchy
        
        # Always include base clusters summary
        for cid, info in base_clusters.items():
            json_hierarchy["base_clusters"].append({
                "id": int(cid),
                "name": info['name'],
                "description": info['description'],
                "abstract_count": cluster_counts.get(int(cid), 0)
            })
        
        return json_hierarchy
    
    def _build_json_hierarchy(self, hierarchy: Dict, 
                            include_full_abstracts: bool) -> List[Dict]:
        """Build JSON representation of hierarchy"""
        if not hierarchy['levels']:
            return []
        
        # Start from top level
        top_level = hierarchy['levels'][-1]['clusters']
        result = []
        
        for cluster_id, cluster in top_level.items():
            cluster_json = self._build_cluster_json(
                cluster, hierarchy, include_full_abstracts
            )
            result.append(cluster_json)
        
        return result
    
    def _build_cluster_json(self, cluster: Dict, hierarchy: Dict,
                          include_full_abstracts: bool) -> Dict:
        """Build JSON for a single cluster and its children"""
        cluster_json = {
            "name": cluster['name'],
            "description": cluster['description'],
            "children": []
        }
        
        children_ids = cluster.get('children', [])
        
        if children_ids:
            # Find child clusters
            child_clusters = None
            for level_info in reversed(hierarchy['levels'][:-1]):
                if any(cid in level_info['clusters'] for cid in map(str, children_ids)):
                    child_clusters = level_info['clusters']
                    break
            
            if child_clusters is None:
                child_clusters = hierarchy['base_clusters']
            
            # Process children
            for child_id in children_ids:
                child_id_str = str(child_id)
                if child_id_str in child_clusters:
                    # Recursive call for sub-cluster
                    child_json = self._build_cluster_json(
                        child_clusters[child_id_str], hierarchy, include_full_abstracts
                    )
                    cluster_json['children'].append(child_json)
                elif child_id_str in self.abstracts:
                    # Leaf node - abstract
                    abstract = self.abstracts[child_id_str]
                    abstract_json = {
                        "pmid": abstract['pmid'],
                        "title": abstract['title']
                    }
                    if include_full_abstracts:
                        abstract_json.update({
                            "abstract": abstract['abstract'],
                            "journal": abstract['journal'],
                            "authors": abstract['authors']
                        })
                    cluster_json['children'].append(abstract_json)
        
        return cluster_json
    
    def export_to_csv(self, hierarchy: Dict, filepath: str = "output/clusters.csv"):
        """Export clusters to CSV format"""
        rows = []
        
        # Flatten hierarchy for CSV
        self._flatten_hierarchy_for_csv(hierarchy, rows)
        
        # Create DataFrame
        df = pd.DataFrame(rows)
        
        # Save to CSV
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        df.to_csv(filepath, index=False, encoding='utf-8')
        print(f"Exported clusters to {filepath}")
        
        return df
    
    def _flatten_hierarchy_for_csv(self, hierarchy: Dict, rows: List[Dict],
                                 parent_path: str = ""):
        """Flatten hierarchy for CSV export"""
        if hierarchy['levels']:
            # Process from top level
            top_level = hierarchy['levels'][-1]['clusters']
            
            for cluster_id, cluster in top_level.items():
                cluster_path = cluster['name']
                self._process_cluster_for_csv(
                    cluster, hierarchy, rows, cluster_path, 1
                )
    
    def _process_cluster_for_csv(self, cluster: Dict, hierarchy: Dict,
                               rows: List[Dict], path: str, level: int):
        """Process a cluster for CSV export"""
        children_ids = cluster.get('children', [])
        
        if children_ids:
            # Find child clusters
            child_clusters = None
            for level_info in reversed(hierarchy['levels'][:-1]):
                if any(cid in level_info['clusters'] for cid in map(str, children_ids)):
                    child_clusters = level_info['clusters']
                    break
            
            if child_clusters is None:
                child_clusters = hierarchy['base_clusters']
            
            # Process children
            for child_id in children_ids:
                child_id_str = str(child_id)
                if child_id_str in child_clusters:
                    child = child_clusters[child_id_str]
                    child_path = f"{path} > {child['name']}"
                    self._process_cluster_for_csv(
                        child, hierarchy, rows, child_path, level + 1
                    )
                elif child_id_str in self.abstracts:
                    # Leaf node - add row
                    abstract = self.abstracts[child_id_str]
                    rows.append({
                        'hierarchy_path': path,
                        'level': level,
                        'cluster_name': cluster['name'],
                        'cluster_description': cluster['description'],
                        'pmid': abstract['pmid'],
                        'title': abstract['title'],
                        'journal': abstract['journal']
                    })
    
    def save_outputs(self, hierarchy: Dict, output_dir: str = "output"):
        """Save all output formats"""
        os.makedirs(output_dir, exist_ok=True)
        
        # Save text format
        text_output = self.format_hierarchy_text(hierarchy)
        with open(os.path.join(output_dir, "hierarchy.txt"), 'w', encoding='utf-8') as f:
            f.write(text_output)
        
        # Save complete hierarchy JSON (new format with top cluster)
        complete_json = self._create_complete_hierarchy_json(hierarchy)
        complete_json = convert_to_serializable(complete_json)
        with open(os.path.join(output_dir, "complete_hierarchy.json"), 'w', encoding='utf-8') as f:
            json.dump(complete_json, f, indent=2, ensure_ascii=False)
        
        # Also save legacy format for compatibility
        json_output = self.format_json_output(hierarchy, include_full_abstracts=False)
        json_output = convert_to_serializable(json_output)
        with open(os.path.join(output_dir, "hierarchy.json"), 'w', encoding='utf-8') as f:
            json.dump(json_output, f, indent=2, ensure_ascii=False)
        
        # Save CSV format
        self.export_to_csv(hierarchy, os.path.join(output_dir, "clusters.csv"))
        
        print(f"\nAll outputs saved to {output_dir}/")
        print("- complete_hierarchy.json: Full hierarchy with top cluster")
        print("- hierarchy.txt: Human-readable text format")
        print("- hierarchy.json: Detailed structure with abstract IDs")
        print("- clusters.csv: Flattened CSV for analysis")
    
    def _create_complete_hierarchy_json(self, hierarchy: Dict) -> Dict:
        """Create complete hierarchy JSON with top cluster and all levels"""
        # Load clustering data for counts
        import pickle
        from collections import Counter
        
        cluster_counts = {}
        cluster_to_abstracts = {}
        
        if self.run_dir and os.path.exists(os.path.join(self.run_dir, 'clustering.pkl')):
            with open(os.path.join(self.run_dir, 'clustering.pkl'), 'rb') as f:
                clustering = pickle.load(f)
            cluster_counts = Counter(clustering['labels'])
            
            # Map clusters to abstracts
            for pmid, label in zip(clustering['pmids'], clustering['labels']):
                if label not in cluster_to_abstracts:
                    cluster_to_abstracts[label] = []
                cluster_to_abstracts[label].append(pmid)
        
        # Build complete structure
        result = {
            "metadata": {
                "total_abstracts": len(self.abstracts),
                "base_clusters": len(hierarchy.get('base_clusters', {})),
                "hierarchy_levels": len(hierarchy.get('levels', [])) + 1
            },
            "top_level_cluster": None,
            "intermediate_levels": [],
            "base_clusters": []
        }
        
        # Determine top cluster
        if 'levels' in hierarchy and hierarchy['levels']:
            # Get from highest level
            top_level = hierarchy['levels'][-1]['clusters']
            if len(top_level) == 1:
                cluster_id, cluster = list(top_level.items())[0]
                result["top_level_cluster"] = {
                    "name": cluster['name'],
                    "description": cluster['description'],
                    "total_abstracts": len(self.abstracts),
                    "direct_children": len(hierarchy.get('base_clusters', {}))
                }
            else:
                # Multiple top clusters
                result["top_level_cluster"] = {
                    "name": "Multiple Top-Level Clusters",
                    "description": f"The hierarchy has {len(top_level)} distinct top-level clusters",
                    "clusters": [
                        {
                            "id": cid,
                            "name": c['name'],
                            "description": c['description']
                        }
                        for cid, c in top_level.items()
                    ]
                }
        else:
            # No levels, single implicit top cluster
            result["top_level_cluster"] = {
                "name": "Artificial Intelligence in Healthcare and Life Sciences",
                "description": "All abstracts converged into a single top-level cluster encompassing diverse AI applications",
                "total_abstracts": len(self.abstracts),
                "direct_children": len(hierarchy.get('base_clusters', {}))
            }
        
        # Add intermediate levels if they exist
        if 'levels' in hierarchy:
            for i, level_data in enumerate(hierarchy['levels'][:-1]):  # Exclude top level
                level_info = {
                    "level": i + 1,
                    "clusters": []
                }
                for cid, cluster in level_data['clusters'].items():
                    level_info["clusters"].append({
                        "id": cid,
                        "name": cluster['name'],
                        "description": cluster['description'],
                        "children": cluster.get('children', [])
                    })
                result["intermediate_levels"].append(level_info)
        
        # Add base clusters with full info
        base_clusters = hierarchy.get('base_clusters', {})
        for cid, info in sorted(base_clusters.items(), 
                               key=lambda x: cluster_counts.get(int(x[0]), 0), 
                               reverse=True):
            cluster_data = {
                "id": int(cid),
                "name": info['name'],
                "description": info['description'],
                "abstract_count": cluster_counts.get(int(cid), 0),
                "pmids": cluster_to_abstracts.get(int(cid), [])
            }
            result["base_clusters"].append(cluster_data)
        
        return result


def test_output_formatter():
    """Test the output formatter"""
    print("Testing Output Formatter...")
    
    # Load test data
    with open('data/test_abstracts.json', 'r') as f:
        abstracts_list = json.load(f)
    
    # Convert to dict
    abstracts_dict = {a['pmid']: a for a in abstracts_list}
    
    # Load hierarchy
    with open('data/test_hierarchy.json', 'r') as f:
        hierarchy = json.load(f)
    
    # Initialize formatter
    formatter = OutputFormatter(abstracts_dict)
    
    # Test text output
    print("\n" + "="*50)
    print("TEXT OUTPUT:")
    print("="*50)
    text_output = formatter.format_hierarchy_text(hierarchy, max_abstracts_per_cluster=2)
    print(text_output)
    
    # Save all outputs
    formatter.save_outputs(hierarchy, "output/test")
    
    return True


if __name__ == "__main__":
    test_output_formatter()
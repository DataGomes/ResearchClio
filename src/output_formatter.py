"""
Output Formatter Module
Formats hierarchical clustering results for display and export
"""

import os
import json
from typing import Dict, List, Optional
import pandas as pd
from datetime import datetime

class OutputFormatter:
    def __init__(self, abstracts_data: Dict[str, Dict]):
        """Initialize with abstracts data for reference"""
        self.abstracts = abstracts_data
    
    def format_hierarchy_text(self, hierarchy: Dict, 
                            indent: str = "  ",
                            show_abstracts: bool = True,
                            max_abstracts_per_cluster: int = 3) -> str:
        """Format hierarchy as indented text"""
        output = []
        output.append("HIERARCHICAL CLUSTERING OF AI ABSTRACTS")
        output.append("=" * 50)
        output.append("")
        
        # Get top level clusters
        if hierarchy['levels']:
            top_level = hierarchy['levels'][-1]['clusters']
            output.append(f"TOP LEVEL ({len(top_level)} clusters):")
            output.append("")
            
            # Recursively format each top cluster
            for cluster_id, cluster in top_level.items():
                self._format_cluster_recursive(
                    cluster, hierarchy, 0, output, indent,
                    show_abstracts, max_abstracts_per_cluster
                )
                output.append("")  # Empty line between top clusters
        
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
        """Format hierarchy as structured JSON"""
        result = {
            "generated_at": datetime.now().isoformat(),
            "total_abstracts": len(self.abstracts),
            "hierarchy": self._build_json_hierarchy(hierarchy, include_full_abstracts)
        }
        return result
    
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
        
        # Save JSON format
        json_output = self.format_json_output(hierarchy, include_full_abstracts=True)
        with open(os.path.join(output_dir, "hierarchy.json"), 'w', encoding='utf-8') as f:
            json.dump(json_output, f, indent=2, ensure_ascii=False)
        
        # Save CSV format
        self.export_to_csv(hierarchy, os.path.join(output_dir, "clusters.csv"))
        
        print(f"\nAll outputs saved to {output_dir}/")
        print("- hierarchy.txt: Human-readable text format")
        print("- hierarchy.json: Structured JSON with full data")
        print("- clusters.csv: Flattened CSV for analysis")


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
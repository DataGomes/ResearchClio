"""
CLIO Pipeline for PubMed Abstracts
Main script to run the complete hierarchical clustering pipeline
"""

import os
import sys
import argparse
import json
import pickle
from datetime import datetime
from dotenv import load_dotenv

# Import all modules
from pubmed_collector import PubMedCollector
from embedder import AbstractEmbedder
from clusterer import AbstractClusterer
from cluster_namer import ClusterNamer
from hierarchizer import Hierarchizer
from output_formatter import OutputFormatter

load_dotenv()

class CLIOPipeline:
    def __init__(self, output_dir: str = "output"):
        """Initialize CLIO pipeline"""
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Initialize timestamp for this run
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.run_dir = os.path.join(output_dir, f"run_{self.timestamp}")
        os.makedirs(self.run_dir, exist_ok=True)
        
        print(f"CLIO Pipeline initialized. Output directory: {self.run_dir}")
    
    def run_pipeline(self, 
                    query: str = '("artificial intelligence" OR "machine learning" OR "deep learning") AND hasabstract',
                    max_abstracts: int = 100,
                    min_cluster_size: int = 5,
                    target_hierarchy_levels: int = 3,
                    top_k_clusters: int = 5,
                    use_cache: bool = True,
                    start_year: int = None,
                    language: str = "eng"):
        """Run the complete CLIO pipeline"""
        
        print("\n" + "="*60)
        print("CLIO PIPELINE FOR PUBMED ABSTRACTS")
        print("="*60)
        print(f"Query: {query}")
        print(f"Max abstracts: {max_abstracts}")
        print(f"Target hierarchy levels: {target_hierarchy_levels}")
        print(f"Top-level clusters: {top_k_clusters}")
        print("="*60 + "\n")
        
        # Step 1: Collect PubMed abstracts
        print("STEP 1: Collecting PubMed abstracts...")
        collector = PubMedCollector()
        
        # Build the full query with year and language filters
        full_query = query
        if start_year:
            full_query = f'{full_query} AND ("{start_year}"[Date - Publication] : "{start_year}"[Date - Publication])'
        if language:
            full_query = f'{full_query} AND {language}[Language]'
        
        abstracts_list = collector.collect_abstracts(
            query=full_query,
            max_results=max_abstracts, 
            use_cache=use_cache
        )
        
        if not abstracts_list:
            print("ERROR: No abstracts collected!")
            return
        
        # Save abstracts
        abstracts_file = os.path.join(self.run_dir, "abstracts.json")
        collector.save_abstracts(abstracts_list, abstracts_file)
        
        # Convert to dict for easy lookup
        abstracts_dict = {a['pmid']: a for a in abstracts_list}
        
        # Step 2: Generate embeddings
        print(f"\nSTEP 2: Generating embeddings for {len(abstracts_list)} abstracts...")
        embedder = AbstractEmbedder()
        embeddings, pmids = embedder.embed_abstracts(abstracts_list)
        
        # Save embeddings
        embeddings_file = os.path.join(self.run_dir, "embeddings.pkl")
        embedder.save_embeddings(embeddings, pmids, embeddings_file)
        
        # Step 3: Perform clustering
        print(f"\nSTEP 3: Clustering abstracts...")
        clusterer = AbstractClusterer(min_cluster_size=min_cluster_size)
        labels = clusterer.cluster_embeddings(embeddings)
        
        # Get cluster assignments
        cluster_abstracts = clusterer.get_cluster_abstracts(labels, pmids)
        
        # Save clustering results
        clustering_file = os.path.join(self.run_dir, "clustering.pkl")
        clusterer.save_clustering_results(labels, pmids, clustering_file)
        
        # Step 4: Generate cluster names
        print(f"\nSTEP 4: Generating names for {len(cluster_abstracts)} clusters...")
        namer = ClusterNamer()
        cluster_names = namer.name_all_clusters(
            cluster_abstracts, abstracts_dict, labels, embeddings, pmids
        )
        
        # Save cluster names
        names_file = os.path.join(self.run_dir, "cluster_names.json")
        namer.save_cluster_names(cluster_names, names_file)
        
        # Step 5: Build hierarchy
        print(f"\nSTEP 5: Building {target_hierarchy_levels}-level hierarchy...")
        hierarchizer = Hierarchizer()
        hierarchy = hierarchizer.build_hierarchy(
            base_clusters=cluster_names,
            target_levels=target_hierarchy_levels,
            top_k=top_k_clusters
        )
        
        # Save hierarchy
        hierarchy_file = os.path.join(self.run_dir, "hierarchy.json")
        hierarchizer.save_hierarchy(hierarchy, hierarchy_file)
        
        # Step 6: Format and export outputs
        print(f"\nSTEP 6: Formatting outputs...")
        formatter = OutputFormatter(abstracts_dict, run_dir=self.run_dir)
        formatter.save_outputs(hierarchy, self.run_dir)
        
        # Print summary
        self.print_summary(abstracts_list, cluster_abstracts, hierarchy)
        
        print(f"\nPipeline complete! All outputs saved to: {self.run_dir}")
        
        return hierarchy
    
    def print_summary(self, abstracts: list, clusters: dict, hierarchy: dict):
        """Print pipeline summary"""
        print("\n" + "="*60)
        print("PIPELINE SUMMARY")
        print("="*60)
        print(f"Total abstracts processed: {len(abstracts)}")
        print(f"Base-level clusters created: {len(clusters)}")
        
        # Cluster size distribution
        cluster_sizes = [len(pmids) for pmids in clusters.values()]
        print(f"Cluster sizes: min={min(cluster_sizes)}, max={max(cluster_sizes)}, "
              f"avg={sum(cluster_sizes)/len(cluster_sizes):.1f}")
        
        # Hierarchy summary
        print(f"\nHierarchy levels: {len(hierarchy['levels']) + 1}")
        for i, level in enumerate(hierarchy['levels']):
            print(f"  Level {i+1}: {len(level['clusters'])} clusters")
        
        # Sample output
        print("\nSample top-level clusters:")
        if hierarchy['levels']:
            top_level = hierarchy['levels'][-1]['clusters']
            for i, (cid, cluster) in enumerate(list(top_level.items())[:3]):
                print(f"  • {cluster['name']}")
                print(f"    {cluster['description'][:80]}...")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Run CLIO hierarchical clustering on PubMed abstracts"
    )
    parser.add_argument(
        "--max-abstracts", 
        type=int, 
        default=100,
        help="Maximum number of abstracts to collect (default: 100)"
    )
    parser.add_argument(
        "--min-cluster-size",
        type=int,
        default=5,
        help="Minimum abstracts per cluster (default: 5)"
    )
    parser.add_argument(
        "--hierarchy-levels",
        type=int,
        default=10,
        help="Maximum hierarchy levels to create (default: 10, algorithm will use fewer if appropriate)"
    )
    parser.add_argument(
        "--top-clusters",
        type=int,
        default=5,
        help="Number of top-level clusters (default: 5)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output",
        help="Output directory (default: output)"
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable caching of PubMed abstracts"
    )
    parser.add_argument(
        "--start-year",
        type=int,
        default=None,
        help="Filter abstracts from this year onwards (e.g., 2020)"
    )
    parser.add_argument(
        "--language",
        type=str,
        default="eng",
        help="Language filter (default: eng for English)"
    )
    
    args = parser.parse_args()
    
    # Run pipeline
    pipeline = CLIOPipeline(output_dir=args.output_dir)
    pipeline.run_pipeline(
        max_abstracts=args.max_abstracts,
        min_cluster_size=args.min_cluster_size,
        target_hierarchy_levels=args.hierarchy_levels,
        top_k_clusters=args.top_clusters,
        use_cache=not args.no_cache,
        start_year=args.start_year,
        language=args.language
    )


if __name__ == "__main__":
    main()
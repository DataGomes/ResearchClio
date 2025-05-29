"""
Clustering Module
Implements k-means clustering for abstract embeddings
"""

import os
import json
import numpy as np
import pickle
from typing import List, Dict, Tuple, Optional
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from collections import Counter
from tqdm import tqdm
import warnings

class AbstractClusterer:
    def __init__(self, min_cluster_size: int = 5):
        """Initialize clusterer with minimum cluster size"""
        self.min_cluster_size = min_cluster_size
        self.kmeans = None
        self.optimal_k = None
        
    def determine_optimal_k(self, embeddings: np.ndarray, 
                          min_k: int = 5, 
                          max_k: int = None,
                          method: str = 'hierarchical') -> int:
        """Determine optimal number of clusters"""
        n_samples = len(embeddings)
        
        # For hierarchical clustering, we need many more base clusters
        if method == 'hierarchical':
            # Dynamic scaling based on dataset size
            # For small datasets (<10k): ~50-100 abstracts per cluster
            # For medium datasets (10k-100k): ~100-200 abstracts per cluster  
            # For large datasets (>100k): ~200-500 abstracts per cluster
            
            if n_samples < 10000:
                # Small datasets: still relatively granular
                target_size = max(50, n_samples // 50)  # Cap at ~200 clusters for 10k docs
            elif n_samples < 50000:
                # Medium datasets (including 40k): target 150-300 clusters
                # For 40k documents: ~130-260 docs/cluster = ~150-300 clusters
                target_size = max(130, min(260, int(n_samples / 200)))
            else:
                # Large datasets: prevent too many clusters
                # For 100k documents: ~250 docs/cluster = ~400 clusters
                target_size = max(250, n_samples // 400)
            
            optimal_k = max(min_k, min(n_samples // target_size, n_samples // self.min_cluster_size))
            
            # Ensure we have enough clusters for meaningful hierarchy
            print(f"Using hierarchical method: {optimal_k} clusters for {n_samples} documents (target size: {target_size} docs/cluster)")
            self.optimal_k = optimal_k
            return optimal_k
        
        # Set max_k based on dataset size
        if max_k is None:
            max_k = min(int(np.sqrt(n_samples)), n_samples // self.min_cluster_size)
            max_k = max(max_k, min_k + 1)
        
        print(f"Determining optimal k between {min_k} and {max_k}...")
        
        if method == 'silhouette':
            scores = []
            k_values = list(range(min_k, max_k + 1))
            
            for k in tqdm(k_values):
                if k >= n_samples:
                    continue
                    
                kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
                labels = kmeans.fit_predict(embeddings)
                
                # Calculate silhouette score
                if len(np.unique(labels)) > 1:
                    score = silhouette_score(embeddings, labels, sample_size=min(1000, n_samples))
                    scores.append(score)
                else:
                    scores.append(-1)
            
            # Find k with highest silhouette score
            best_idx = np.argmax(scores)
            optimal_k = k_values[best_idx]
            
            print(f"Silhouette scores: {scores}")
            print(f"Optimal k: {optimal_k} (silhouette score: {scores[best_idx]:.4f})")
            
        else:  # Use sqrt heuristic
            optimal_k = int(np.sqrt(n_samples))
            optimal_k = max(min_k, min(optimal_k, max_k))
            print(f"Using sqrt heuristic: k = {optimal_k}")
        
        self.optimal_k = optimal_k
        return optimal_k
    
    def cluster_embeddings(self, embeddings: np.ndarray, 
                          k: Optional[int] = None) -> np.ndarray:
        """Perform k-means clustering on embeddings"""
        if k is None:
            k = self.determine_optimal_k(embeddings)
        
        print(f"Clustering {len(embeddings)} embeddings into {k} clusters...")
        
        self.kmeans = KMeans(
            n_clusters=k,
            random_state=42,
            n_init=10,
            max_iter=300
        )
        
        labels = self.kmeans.fit_predict(embeddings)
        
        # Print cluster statistics
        cluster_counts = Counter(labels)
        print(f"\nCluster distribution:")
        for cluster_id, count in sorted(cluster_counts.items()):
            print(f"  Cluster {cluster_id}: {count} abstracts")
        
        # Check for small clusters
        small_clusters = [cid for cid, count in cluster_counts.items() 
                         if count < self.min_cluster_size]
        if small_clusters:
            print(f"\nWarning: {len(small_clusters)} clusters have less than "
                  f"{self.min_cluster_size} abstracts")
        
        return labels
    
    def get_cluster_centers(self) -> np.ndarray:
        """Get cluster centers"""
        if self.kmeans is None:
            raise ValueError("Must run clustering first")
        return self.kmeans.cluster_centers_
    
    def get_cluster_abstracts(self, labels: np.ndarray, 
                            pmids: List[str]) -> Dict[int, List[str]]:
        """Group PMIDs by cluster"""
        cluster_abstracts = {}
        
        for pmid, label in zip(pmids, labels):
            if label not in cluster_abstracts:
                cluster_abstracts[label] = []
            cluster_abstracts[label].append(pmid)
        
        return cluster_abstracts
    
    def find_representative_abstracts(self, embeddings: np.ndarray, 
                                    labels: np.ndarray, 
                                    pmids: List[str],
                                    n_representatives: int = 5) -> Dict[int, List[str]]:
        """Find most representative abstracts for each cluster (closest to center)"""
        if self.kmeans is None:
            raise ValueError("Must run clustering first")
        
        representatives = {}
        centers = self.get_cluster_centers()
        
        for cluster_id in range(len(centers)):
            # Get indices of abstracts in this cluster
            cluster_indices = np.where(labels == cluster_id)[0]
            
            if len(cluster_indices) == 0:
                representatives[cluster_id] = []
                continue
            
            # Get embeddings for this cluster
            cluster_embeddings = embeddings[cluster_indices]
            cluster_pmids = [pmids[i] for i in cluster_indices]
            
            # Calculate distances to cluster center
            center = centers[cluster_id]
            distances = np.linalg.norm(cluster_embeddings - center, axis=1)
            
            # Get indices of closest abstracts
            n_reps = min(n_representatives, len(cluster_indices))
            closest_indices = np.argsort(distances)[:n_reps]
            
            representatives[cluster_id] = [cluster_pmids[i] for i in closest_indices]
        
        return representatives
    
    def save_clustering_results(self, labels: np.ndarray, pmids: List[str],
                              filepath: str = "data/clustering_results.pkl"):
        """Save clustering results"""
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        results = {
            'labels': labels,
            'pmids': pmids,
            'n_clusters': self.optimal_k,
            'cluster_centers': self.get_cluster_centers() if self.kmeans else None
        }
        
        with open(filepath, 'wb') as f:
            pickle.dump(results, f)
        
        print(f"Saved clustering results to {filepath}")
    
    def load_clustering_results(self, filepath: str) -> Tuple[np.ndarray, List[str]]:
        """Load clustering results"""
        with open(filepath, 'rb') as f:
            results = pickle.load(f)
        
        return results['labels'], results['pmids']


def test_clusterer():
    """Test the clustering module"""
    print("Testing Abstract Clusterer...")
    
    # Load embeddings from previous step
    with open('data/test_embeddings.pkl', 'rb') as f:
        data = pickle.load(f)
    
    embeddings = data['embeddings']
    pmids = data['pmids']
    
    # Initialize clusterer
    clusterer = AbstractClusterer(min_cluster_size=2)  # Lower threshold for test data
    
    # Perform clustering with automatic k selection
    labels = clusterer.cluster_embeddings(embeddings)
    
    # Get cluster assignments
    cluster_abstracts = clusterer.get_cluster_abstracts(labels, pmids)
    print(f"\nCluster assignments:")
    for cluster_id, cluster_pmids in cluster_abstracts.items():
        print(f"  Cluster {cluster_id}: {cluster_pmids}")
    
    # Find representative abstracts
    representatives = clusterer.find_representative_abstracts(embeddings, labels, pmids, n_representatives=2)
    print(f"\nRepresentative abstracts per cluster:")
    for cluster_id, rep_pmids in representatives.items():
        print(f"  Cluster {cluster_id}: {rep_pmids}")
    
    # Save results
    clusterer.save_clustering_results(labels, pmids, "data/test_clustering.pkl")
    
    # Test loading
    loaded_labels, loaded_pmids = clusterer.load_clustering_results("data/test_clustering.pkl")
    assert np.array_equal(labels, loaded_labels), "Loaded labels don't match!"
    assert pmids == loaded_pmids, "Loaded PMIDs don't match!"
    print("\nClustering save/load test passed!")
    
    return True


if __name__ == "__main__":
    test_clusterer()
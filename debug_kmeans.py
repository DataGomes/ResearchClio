#!/usr/bin/env python3
"""Debug K-means implementation issues"""

import numpy as np
import pickle
import os
from sklearn.cluster import KMeans
import warnings

# Show warnings but don't convert to errors
warnings.filterwarnings('default')

def test_kmeans_algorithms(embeddings, k=5):
    """Test different K-means algorithms"""
    print(f"\nTesting different K-means algorithms with k={k}...")
    
    algorithms = ['lloyd', 'elkan', 'auto']
    
    for algo in algorithms:
        print(f"\nAlgorithm: {algo}")
        try:
            kmeans = KMeans(
                n_clusters=k, 
                random_state=42, 
                n_init=10,
                algorithm=algo,
                max_iter=300
            )
            labels = kmeans.fit_predict(embeddings)
            print(f"  Success! Inertia: {kmeans.inertia_:.4f}")
            print(f"  Iterations: {kmeans.n_iter_}")
            print(f"  Cluster sizes: {np.bincount(labels)}")
            return labels
        except Exception as e:
            print(f"  Failed: {e}")
    
    return None

def test_minibatch_kmeans(embeddings, k=5):
    """Test MiniBatchKMeans as alternative"""
    from sklearn.cluster import MiniBatchKMeans
    
    print(f"\nTesting MiniBatchKMeans with k={k}...")
    try:
        kmeans = MiniBatchKMeans(
            n_clusters=k,
            random_state=42,
            batch_size=100,
            n_init=10
        )
        labels = kmeans.fit_predict(embeddings)
        print(f"  Success! Inertia: {kmeans.inertia_:.4f}")
        print(f"  Cluster sizes: {np.bincount(labels)}")
        return labels
    except Exception as e:
        print(f"  Failed: {e}")
        return None

def test_custom_distance(embeddings, k=5):
    """Test with precomputed distance matrix"""
    from sklearn.metrics.pairwise import cosine_distances
    from sklearn.cluster import AgglomerativeClustering
    
    print(f"\nTesting AgglomerativeClustering with cosine distance, k={k}...")
    try:
        # Compute cosine distance matrix
        distances = cosine_distances(embeddings)
        
        # Use hierarchical clustering
        clustering = AgglomerativeClustering(
            n_clusters=k,
            metric='precomputed',
            linkage='average'
        )
        labels = clustering.fit_predict(distances)
        print(f"  Success!")
        print(f"  Cluster sizes: {np.bincount(labels)}")
        return labels
    except Exception as e:
        print(f"  Failed: {e}")
        return None

def main():
    # Find the most recent embeddings file
    output_dirs = [d for d in os.listdir('output') if d.startswith('run_')]
    if not output_dirs:
        print("No output directories found!")
        return
    
    latest_dir = sorted(output_dirs)[-1]
    embeddings_path = os.path.join('output', latest_dir, 'embeddings.pkl')
    
    print(f"Loading embeddings from: {embeddings_path}")
    
    # Load embeddings
    with open(embeddings_path, 'rb') as f:
        data = pickle.load(f)
    
    embeddings = data['embeddings']
    
    print(f"Loaded {len(embeddings)} embeddings")
    print(f"Shape: {embeddings.shape}")
    print(f"Dtype: {embeddings.dtype}")
    
    # Test different approaches
    print("\n" + "="*50)
    print("TESTING DIFFERENT CLUSTERING APPROACHES")
    print("="*50)
    
    # Test K-means with different algorithms
    labels = test_kmeans_algorithms(embeddings, k=10)
    
    # Test MiniBatchKMeans
    if labels is None:
        labels = test_minibatch_kmeans(embeddings, k=10)
    
    # Test hierarchical clustering
    if labels is None:
        labels = test_custom_distance(embeddings, k=10)
    
    # Check sklearn version
    print("\n" + "="*50)
    print("ENVIRONMENT INFO")
    print("="*50)
    import sklearn
    print(f"Scikit-learn version: {sklearn.__version__}")
    
    import sys
    print(f"Python version: {sys.version}")
    
    # Try with float64
    print("\n" + "="*50)
    print("TESTING WITH FLOAT64")
    print("="*50)
    embeddings_64 = embeddings.astype(np.float64)
    test_kmeans_algorithms(embeddings_64, k=10)

if __name__ == "__main__":
    main()
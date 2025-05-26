#!/usr/bin/env python3
"""Debug script to investigate embedding and clustering issues"""

import numpy as np
import pickle
import os
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import warnings

# Enable all warnings
warnings.filterwarnings('error')

def check_embeddings(embeddings):
    """Check embeddings for numerical issues"""
    print(f"Embeddings shape: {embeddings.shape}")
    print(f"Embeddings dtype: {embeddings.dtype}")
    
    # Check for NaN values
    nan_count = np.sum(np.isnan(embeddings))
    print(f"NaN values: {nan_count}")
    
    # Check for infinite values
    inf_count = np.sum(np.isinf(embeddings))
    print(f"Inf values: {inf_count}")
    
    # Check for zero vectors
    zero_vectors = np.sum(np.all(embeddings == 0, axis=1))
    print(f"Zero vectors: {zero_vectors}")
    
    # Check value ranges
    print(f"Min value: {np.min(embeddings)}")
    print(f"Max value: {np.max(embeddings)}")
    print(f"Mean: {np.mean(embeddings)}")
    print(f"Std: {np.std(embeddings)}")
    
    # Check norms
    norms = np.linalg.norm(embeddings, axis=1)
    print(f"\nNorms - Min: {np.min(norms)}, Max: {np.max(norms)}, Mean: {np.mean(norms)}")
    
    # Check for duplicate embeddings
    unique_embeddings = np.unique(embeddings, axis=0)
    print(f"\nDuplicate embeddings: {len(embeddings) - len(unique_embeddings)}")
    
    return nan_count == 0 and inf_count == 0 and zero_vectors == 0

def test_kmeans_on_embeddings(embeddings, k=5):
    """Test K-means clustering on embeddings"""
    print(f"\nTesting K-means with k={k}...")
    
    try:
        # First try without any preprocessing
        print("1. Testing raw embeddings...")
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=1)
        labels = kmeans.fit_predict(embeddings)
        print(f"   Success! Unique labels: {np.unique(labels)}")
        
        # Try computing silhouette score
        if len(np.unique(labels)) > 1:
            score = silhouette_score(embeddings, labels)
            print(f"   Silhouette score: {score:.4f}")
    except Exception as e:
        print(f"   Failed: {e}")
        
        # Try with normalized embeddings
        print("\n2. Testing with L2 normalization...")
        from sklearn.preprocessing import normalize
        normalized = normalize(embeddings, norm='l2')
        
        # Check if normalization worked
        norms = np.linalg.norm(normalized, axis=1)
        print(f"   Norms after normalization - Min: {np.min(norms):.6f}, Max: {np.max(norms):.6f}")
        
        try:
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=1)
            labels = kmeans.fit_predict(normalized)
            print(f"   Success! Unique labels: {np.unique(labels)}")
            
            if len(np.unique(labels)) > 1:
                score = silhouette_score(normalized, labels)
                print(f"   Silhouette score: {score:.4f}")
        except Exception as e2:
            print(f"   Still failed: {e2}")
            
            # Try removing problematic vectors
            print("\n3. Testing after removing problematic vectors...")
            # Remove zero vectors
            non_zero_mask = ~np.all(embeddings == 0, axis=1)
            clean_embeddings = embeddings[non_zero_mask]
            print(f"   Removed {np.sum(~non_zero_mask)} zero vectors")
            
            # Remove vectors with very small norms
            norms = np.linalg.norm(clean_embeddings, axis=1)
            valid_mask = norms > 1e-6
            clean_embeddings = clean_embeddings[valid_mask]
            print(f"   Removed {np.sum(~valid_mask)} vectors with tiny norms")
            
            # Normalize
            clean_normalized = normalize(clean_embeddings, norm='l2')
            
            try:
                kmeans = KMeans(n_clusters=min(k, len(clean_normalized)), random_state=42, n_init=1)
                labels = kmeans.fit_predict(clean_normalized)
                print(f"   Success! Unique labels: {np.unique(labels)}")
            except Exception as e3:
                print(f"   Still failed: {e3}")

def main():
    # Find the most recent embeddings file
    output_dirs = [d for d in os.listdir('output') if d.startswith('run_')]
    if not output_dirs:
        print("No output directories found!")
        return
    
    latest_dir = sorted(output_dirs)[-1]
    embeddings_path = os.path.join('output', latest_dir, 'embeddings.pkl')
    
    print(f"Loading embeddings from: {embeddings_path}")
    
    if not os.path.exists(embeddings_path):
        print("Embeddings file not found!")
        return
    
    # Load embeddings
    with open(embeddings_path, 'rb') as f:
        data = pickle.load(f)
    
    embeddings = data['embeddings']
    pmids = data['pmids']
    
    print(f"\nLoaded {len(embeddings)} embeddings for {len(pmids)} PMIDs")
    
    # Check embeddings
    print("\n" + "="*50)
    print("CHECKING EMBEDDINGS")
    print("="*50)
    is_valid = check_embeddings(embeddings)
    
    # Test clustering
    print("\n" + "="*50)
    print("TESTING CLUSTERING")
    print("="*50)
    test_kmeans_on_embeddings(embeddings, k=5)
    
    # Test with different k values
    print("\n" + "="*50)
    print("TESTING WITH DIFFERENT K VALUES")
    print("="*50)
    for k in [2, 10, 20]:
        print(f"\n--- k={k} ---")
        test_kmeans_on_embeddings(embeddings, k=k)

if __name__ == "__main__":
    main()
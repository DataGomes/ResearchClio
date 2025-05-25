"""
Embedding Module
Generates embeddings for abstracts using Sentence-BERT (all-mpnet-base-v2)
"""

import os
import json
import numpy as np
from typing import List, Dict, Tuple
from sentence_transformers import SentenceTransformer
from tqdm import tqdm
import pickle

class AbstractEmbedder:
    def __init__(self, model_name: str = 'all-mpnet-base-v2'):
        """Initialize embedder with specified model"""
        print(f"Loading embedding model: {model_name}")
        
        # Set offline mode to use cached model
        import os
        os.environ['TRANSFORMERS_OFFLINE'] = '1'
        os.environ['HF_HUB_OFFLINE'] = '1'
        
        try:
            self.model = SentenceTransformer(model_name)
        except Exception as e:
            print(f"Warning: {e}")
            print("Attempting to load model without checking for updates...")
            # Force local loading
            self.model = SentenceTransformer(model_name, cache_folder=os.path.expanduser("~/.cache/huggingface/hub"))
            
        self.embedding_dim = self.model.get_sentence_embedding_dimension()
        print(f"Model loaded. Embedding dimension: {self.embedding_dim}")
        
    def create_text_for_embedding(self, abstract_data: Dict) -> str:
        """Combine title and abstract for embedding"""
        title = abstract_data.get('title', '')
        abstract = abstract_data.get('abstract', '')
        # Combine title and abstract with a separator
        combined_text = f"{title} [SEP] {abstract}"
        return combined_text
    
    def embed_abstracts(self, abstracts: List[Dict], batch_size: int = 32) -> Tuple[np.ndarray, List[str]]:
        """
        Generate embeddings for a list of abstracts
        Returns: (embeddings array, list of PMIDs)
        """
        print(f"Generating embeddings for {len(abstracts)} abstracts...")
        
        # Prepare texts and PMIDs
        texts = []
        pmids = []
        
        for abstract in abstracts:
            if abstract.get('abstract'):  # Only process if abstract exists
                texts.append(self.create_text_for_embedding(abstract))
                pmids.append(abstract['pmid'])
        
        # Generate embeddings
        embeddings = self.model.encode(
            texts,
            batch_size=batch_size,
            show_progress_bar=True,
            convert_to_numpy=True
        )
        
        print(f"Generated {len(embeddings)} embeddings")
        return embeddings, pmids
    
    def save_embeddings(self, embeddings: np.ndarray, pmids: List[str], 
                       filepath: str = "data/abstract_embeddings.pkl"):
        """Save embeddings and PMIDs to file"""
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        data = {
            'embeddings': embeddings,
            'pmids': pmids,
            'model_name': self.model.get_sentence_embedding_dimension(),
            'embedding_dim': self.embedding_dim
        }
        
        with open(filepath, 'wb') as f:
            pickle.dump(data, f)
            
        print(f"Saved embeddings to {filepath}")
        
    def load_embeddings(self, filepath: str = "data/abstract_embeddings.pkl") -> Tuple[np.ndarray, List[str]]:
        """Load embeddings from file"""
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
            
        return data['embeddings'], data['pmids']
    
    def compute_similarity(self, embedding1: np.ndarray, embedding2: np.ndarray) -> float:
        """Compute cosine similarity between two embeddings"""
        dot_product = np.dot(embedding1, embedding2)
        norm1 = np.linalg.norm(embedding1)
        norm2 = np.linalg.norm(embedding2)
        
        return dot_product / (norm1 * norm2)
    
    def find_similar_abstracts(self, query_pmid: str, embeddings: np.ndarray, 
                             pmids: List[str], top_k: int = 5) -> List[Tuple[str, float]]:
        """Find most similar abstracts to a query abstract"""
        # Find query index
        try:
            query_idx = pmids.index(query_pmid)
        except ValueError:
            print(f"PMID {query_pmid} not found")
            return []
        
        query_embedding = embeddings[query_idx]
        
        # Compute similarities
        similarities = []
        for i, (emb, pmid) in enumerate(zip(embeddings, pmids)):
            if i != query_idx:  # Skip self
                sim = self.compute_similarity(query_embedding, emb)
                similarities.append((pmid, sim))
        
        # Sort by similarity
        similarities.sort(key=lambda x: x[1], reverse=True)
        
        return similarities[:top_k]


def test_embedder():
    """Test the embedding module"""
    print("Testing Abstract Embedder...")
    
    # Load test abstracts
    with open('data/test_abstracts.json', 'r') as f:
        abstracts = json.load(f)
    
    # Initialize embedder
    embedder = AbstractEmbedder()
    
    # Generate embeddings
    embeddings, pmids = embedder.embed_abstracts(abstracts)
    
    print(f"\nEmbedding statistics:")
    print(f"Shape: {embeddings.shape}")
    print(f"Mean: {np.mean(embeddings):.4f}")
    print(f"Std: {np.std(embeddings):.4f}")
    
    # Save embeddings
    embedder.save_embeddings(embeddings, pmids, "data/test_embeddings.pkl")
    
    # Test similarity search
    if len(pmids) > 1:
        print(f"\nFinding similar abstracts to PMID {pmids[0]}...")
        similar = embedder.find_similar_abstracts(pmids[0], embeddings, pmids, top_k=3)
        
        for pmid, sim in similar:
            print(f"  PMID {pmid}: similarity = {sim:.4f}")
    
    # Test loading
    loaded_embeddings, loaded_pmids = embedder.load_embeddings("data/test_embeddings.pkl")
    assert np.array_equal(embeddings, loaded_embeddings), "Loaded embeddings don't match!"
    assert pmids == loaded_pmids, "Loaded PMIDs don't match!"
    print("\nEmbedding save/load test passed!")
    
    return True


if __name__ == "__main__":
    test_embedder()
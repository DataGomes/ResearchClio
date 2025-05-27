"""
Cluster Naming Module
Uses Claude API to generate descriptive names for clusters
"""

import os
import json
import numpy as np
from typing import List, Dict, Tuple, Optional
from anthropic import Anthropic
from dotenv import load_dotenv
import random
import pickle

load_dotenv()

class ClusterNamer:
    def __init__(self, model: str = None):
        """Initialize with Claude API"""
        self.api_key = os.getenv('ANTHROPIC_API_KEY')
        if not self.api_key:
            raise ValueError("ANTHROPIC_API_KEY not found in environment variables")
        
        self.client = Anthropic(api_key=self.api_key)
        self.model = model or os.getenv('CLAUDE_MODEL', 'claude-sonnet-4-20250514')
        print(f"Using Claude model: {self.model}")
        
    def sample_abstracts(self, cluster_pmids: List[str], 
                        all_abstracts: Dict[str, Dict],
                        n_samples: int = 10) -> List[Dict]:
        """Sample abstracts from cluster"""
        # Sample up to n_samples from cluster
        sampled_pmids = random.sample(cluster_pmids, min(n_samples, len(cluster_pmids)))
        sampled_abstracts = [all_abstracts[pmid] for pmid in sampled_pmids if pmid in all_abstracts]
        return sampled_abstracts
    
    def find_nearest_non_cluster_abstracts(self, cluster_id: int,
                                         labels: np.ndarray,
                                         embeddings: np.ndarray,
                                         pmids: List[str],
                                         all_abstracts: Dict[str, Dict],
                                         n_samples: int = 10) -> List[Dict]:
        """Find nearest abstracts not in the cluster"""
        # Get cluster and non-cluster indices
        cluster_indices = np.where(labels == cluster_id)[0]
        non_cluster_indices = np.where(labels != cluster_id)[0]
        
        if len(non_cluster_indices) == 0:
            return []
        
        # Calculate mean embedding for cluster
        cluster_embeddings = embeddings[cluster_indices]
        cluster_center = np.mean(cluster_embeddings, axis=0)
        
        # Calculate distances from non-cluster abstracts to cluster center
        non_cluster_embeddings = embeddings[non_cluster_indices]
        distances = np.linalg.norm(non_cluster_embeddings - cluster_center, axis=1)
        
        # Get nearest non-cluster abstracts
        n_nearest = min(n_samples, len(non_cluster_indices))
        nearest_indices = np.argsort(distances)[:n_nearest]
        nearest_pmids = [pmids[non_cluster_indices[i]] for i in nearest_indices]
        
        return [all_abstracts[pmid] for pmid in nearest_pmids if pmid in all_abstracts]
    
    def format_abstract_for_prompt(self, abstract: Dict) -> str:
        """Format abstract for Claude prompt"""
        title = abstract.get('title', 'No title')
        text = abstract.get('abstract', 'No abstract')[:500]  # Limit length
        return f"Title: {title}\nAbstract: {text}...\n"
    
    def generate_cluster_name(self, cluster_abstracts: List[Dict],
                            non_cluster_abstracts: List[Dict]) -> Dict:
        """Generate name and description for a cluster using Claude"""
        # Format abstracts for prompt
        cluster_texts = [self.format_abstract_for_prompt(a) for a in cluster_abstracts[:10]]
        non_cluster_texts = [self.format_abstract_for_prompt(a) for a in non_cluster_abstracts[:10]]
        
        prompt = f"""I have a cluster of scientific abstracts about artificial intelligence. Your task is to identify what makes this cluster unique and generate a name and description.

CLUSTER ABSTRACTS (these are all in the same cluster):
{'-' * 50}
{''.join(cluster_texts)}

NEARBY NON-CLUSTER ABSTRACTS (these are similar but NOT in the cluster):
{'-' * 50}
{''.join(non_cluster_texts)}

Based on what distinguishes the cluster abstracts from the nearby non-cluster abstracts, generate:

1. A descriptive name (maximum 10 words) that captures the unique theme or focus of this cluster
2. A 2-sentence description explaining what unifies the abstracts in this cluster

Respond with ONLY a JSON object (no other text, no explanations) in this exact format:
{{"name": "Your cluster name here", "description": "Your 2-sentence description here."}}

Example of a good response:
{{"name": "Deep Learning for Medical Image Analysis", "description": "This cluster focuses on using convolutional neural networks to analyze medical images. The applications range from disease detection to treatment planning across various imaging modalities."}}"""

        try:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=200,
                temperature=1.0,
                messages=[
                    {"role": "user", "content": prompt}
                ]
            )
            
            # Parse JSON response
            response_text = response.content[0].text.strip()
            
            # Clean up response if needed
            if response_text.startswith('```json'):
                response_text = response_text[7:]
            if response_text.endswith('```'):
                response_text = response_text[:-3]
            
            # Try to extract JSON from the response using regex
            # This handles cases where Claude adds extra text or returns partial JSON
            import re
            
            # First check if we have a partial response (missing opening brace)
            if response_text and not response_text.strip().startswith('{'):
                # Check if it looks like a JSON continuation
                if '"name"' in response_text or response_text.strip().startswith('"'):
                    response_text = '{' + response_text
            
            # Look for complete JSON object
            json_pattern = r'\{[^{}]*"name"[^{}]*:[^{}]*"[^"]*"[^{}]*,\s*"description"[^{}]*:[^{}]*"[^"]*"[^{}]*\}'
            json_match = re.search(json_pattern, response_text, re.DOTALL)
            if json_match:
                response_text = json_match.group()
            
            result = json.loads(response_text.strip())
            
            # Validate response
            if 'name' not in result or 'description' not in result:
                raise ValueError("Invalid response format")
            
            # Ensure name is not too long
            name_words = result['name'].split()
            if len(name_words) > 10:
                result['name'] = ' '.join(name_words[:10])
            
            return result
            
        except json.JSONDecodeError as e:
            print(f"Error parsing JSON response: {e}")
            print(f"Raw response: {response_text[:500] if 'response_text' in locals() else 'No response'}")
            
            # Try one more time with a stricter prompt
            try:
                retry_prompt = prompt + "\n\nIMPORTANT: Return ONLY the JSON object, nothing else."
                retry_response = self.client.messages.create(
                    model=self.model,
                    max_tokens=200,
                    temperature=0.5,  # Lower temperature for more consistent output
                    messages=[
                        {"role": "user", "content": retry_prompt},
                        {"role": "assistant", "content": "{\"name\": \""}
                    ]
                )
                
                retry_text = retry_response.content[0].text.strip()
                # Complete the JSON that we started
                retry_text = "{\"name\": \"" + retry_text
                
                result = json.loads(retry_text)
                return result
                
            except Exception as retry_e:
                print(f"Retry also failed: {retry_e}")
                raise RuntimeError(f"Failed to generate cluster name after retry: {e}")
                
        except Exception as e:
            print(f"Error generating cluster name: {e}")
            raise RuntimeError(f"Failed to generate cluster name using Claude API: {e}")
    
    def name_all_clusters(self, cluster_abstracts: Dict[int, List[str]],
                         all_abstracts: Dict[str, Dict],
                         labels: np.ndarray,
                         embeddings: np.ndarray,
                         pmids: List[str]) -> Dict[int, Dict]:
        """Generate names for all clusters"""
        cluster_names = {}
        
        print(f"Generating names for {len(cluster_abstracts)} clusters...")
        
        for cluster_id, cluster_pmids in cluster_abstracts.items():
            print(f"\nNaming cluster {cluster_id} ({len(cluster_pmids)} abstracts)...")
            
            # Sample cluster abstracts
            cluster_samples = self.sample_abstracts(cluster_pmids, all_abstracts)
            
            # Find nearest non-cluster abstracts
            non_cluster_samples = self.find_nearest_non_cluster_abstracts(
                cluster_id, labels, embeddings, pmids, all_abstracts
            )
            
            # Generate name and description
            cluster_info = self.generate_cluster_name(cluster_samples, non_cluster_samples)
            cluster_names[cluster_id] = cluster_info
            
            print(f"  Name: {cluster_info['name']}")
            print(f"  Description: {cluster_info['description'][:100]}...")
        
        return cluster_names
    
    def save_cluster_names(self, cluster_names: Dict[int, Dict],
                          filepath: str = "data/cluster_names.json"):
        """Save cluster names to file"""
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        # Convert numpy int keys to regular int for JSON serialization
        serializable_names = {int(k): v for k, v in cluster_names.items()}
        
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(serializable_names, f, indent=2, ensure_ascii=False)
        
        print(f"Saved cluster names to {filepath}")


def test_cluster_namer():
    """Test the cluster naming module"""
    print("Testing Cluster Namer...")
    
    # Load test data
    with open('data/test_abstracts.json', 'r') as f:
        abstracts_list = json.load(f)
    
    # Convert to dict for easy lookup
    all_abstracts = {a['pmid']: a for a in abstracts_list}
    
    # Load clustering results
    with open('data/test_clustering.pkl', 'rb') as f:
        clustering_data = pickle.load(f)
    
    labels = clustering_data['labels']
    pmids = clustering_data['pmids']
    
    # Load embeddings
    with open('data/test_embeddings.pkl', 'rb') as f:
        embedding_data = pickle.load(f)
    embeddings = embedding_data['embeddings']
    
    # Group abstracts by cluster
    cluster_abstracts = {}
    for pmid, label in zip(pmids, labels):
        if label not in cluster_abstracts:
            cluster_abstracts[label] = []
        cluster_abstracts[label].append(pmid)
    
    # Initialize namer
    namer = ClusterNamer()
    
    # Test naming a single cluster
    test_cluster_id = list(cluster_abstracts.keys())[0]
    test_cluster_pmids = cluster_abstracts[test_cluster_id]
    
    print(f"\nTesting single cluster naming (Cluster {test_cluster_id})...")
    cluster_samples = namer.sample_abstracts(test_cluster_pmids, all_abstracts, n_samples=3)
    non_cluster_samples = namer.find_nearest_non_cluster_abstracts(
        test_cluster_id, labels, embeddings, pmids, all_abstracts, n_samples=3
    )
    
    cluster_info = namer.generate_cluster_name(cluster_samples, non_cluster_samples)
    print(f"Generated name: {cluster_info['name']}")
    print(f"Description: {cluster_info['description']}")
    
    # Name all clusters
    print("\nNaming all clusters...")
    cluster_names = namer.name_all_clusters(
        cluster_abstracts, all_abstracts, labels, embeddings, pmids
    )
    
    # Save results
    namer.save_cluster_names(cluster_names, "data/test_cluster_names.json")
    
    return True


if __name__ == "__main__":
    test_cluster_namer()
"""
Hierarchizer Module
Builds hierarchical clustering structure following CLIO approach
"""

import os
import json
import numpy as np
from typing import List, Dict, Tuple, Optional
from sentence_transformers import SentenceTransformer
from sklearn.cluster import KMeans
from anthropic import Anthropic
from dotenv import load_dotenv
import pickle

load_dotenv()

class Hierarchizer:
    def __init__(self, model_name: str = 'all-mpnet-base-v2'):
        """Initialize Hierarchizer with embedding model and Claude API"""
        # Load embedding model
        self.embedding_model = SentenceTransformer(model_name)
        
        # Initialize Claude API
        self.api_key = os.getenv('ANTHROPIC_API_KEY')
        if not self.api_key:
            raise ValueError("ANTHROPIC_API_KEY not found")
        
        self.client = Anthropic(api_key=self.api_key)
        self.claude_model = os.getenv('CLAUDE_MODEL', 'claude-sonnet-4-20250514')
        
    def embed_cluster_descriptions(self, cluster_info: Dict[int, Dict]) -> np.ndarray:
        """Embed cluster names and descriptions"""
        texts = []
        cluster_ids = sorted(cluster_info.keys())
        
        for cid in cluster_ids:
            name = cluster_info[cid]['name']
            desc = cluster_info[cid]['description']
            # Combine name and description for embedding
            texts.append(f"{name}. {desc}")
        
        embeddings = self.embedding_model.encode(texts, show_progress_bar=True)
        return embeddings, cluster_ids
    
    def group_into_neighborhoods(self, embeddings: np.ndarray, 
                               target_per_neighborhood: int = None) -> np.ndarray:
        """Group clusters into neighborhoods for hierarchical processing"""
        n_clusters = len(embeddings)
        n_neighborhoods = max(1, n_clusters // target_per_neighborhood)
        
        if n_neighborhoods == 1:
            # All clusters in one neighborhood
            return np.zeros(n_clusters, dtype=int)
        
        # Use k-means to create neighborhoods
        kmeans = KMeans(n_clusters=n_neighborhoods, random_state=42)
        neighborhood_labels = kmeans.fit_predict(embeddings)
        
        return neighborhood_labels
    
    def propose_parent_clusters(self, child_clusters: List[Dict],
                              target_n_parents: int) -> List[Dict]:
        """Use Claude to propose parent clusters for a neighborhood"""
        # Format child clusters for prompt
        children_text = ""
        for i, child in enumerate(child_clusters):
            children_text += f"{i+1}. {child['name']}: {child['description']}\n"
        
        prompt = f"""I have {len(child_clusters)} clusters of AI research abstracts. I need to group them into approximately {target_n_parents} higher-level parent clusters.

Current clusters:
{children_text}

Please propose {target_n_parents} parent clusters that would logically group these clusters. Each parent should:
1. Have a descriptive name (max 10 words)
2. Have a 2-sentence description
3. Be thematically coherent

Respond with a JSON array of parent clusters in this format:
[
  {{"name": "Parent cluster name", "description": "Two sentence description.", "suggested_children": [1, 2, 3]}}
]

The suggested_children field should contain the numbers of child clusters that belong to this parent."""

        try:
            response = self.client.messages.create(
                model=self.claude_model,
                max_tokens=1000,
                temperature=1.0,
                messages=[{"role": "user", "content": prompt}]
            )
            
            # Debug the response
            if not response.content:
                raise RuntimeError("Claude API returned empty content")
            
            response_text = response.content[0].text.strip()
            
            if not response_text:
                raise RuntimeError("Claude API returned empty text")
            
            if response_text.startswith('```json'):
                response_text = response_text[7:]
            if response_text.endswith('```'):
                response_text = response_text[:-3]
            
            parent_clusters = json.loads(response_text.strip())
            return parent_clusters
            
        except json.JSONDecodeError as e:
            print(f"Error parsing Claude response as JSON: {e}")
            print(f"Raw response: {response_text[:500]}...")  # Show first 500 chars
            raise RuntimeError(f"Failed to parse Claude response as JSON: {e}")
        except Exception as e:
            print(f"Error proposing parent clusters: {e}")
            raise RuntimeError(f"Failed to generate parent clusters using Claude API: {e}")
    
    def deduplicate_clusters(self, all_proposed: List[Dict]) -> List[Dict]:
        """Deduplicate similar parent clusters across neighborhoods"""
        if len(all_proposed) <= 1:
            return all_proposed
        
        # Embed all proposed parents
        texts = [f"{p['name']}. {p['description']}" for p in all_proposed]
        embeddings = self.embedding_model.encode(texts)
        
        # Use cosine similarity to find duplicates
        unique_clusters = []
        used_indices = set()
        
        for i, cluster in enumerate(all_proposed):
            if i in used_indices:
                continue
                
            # Find similar clusters
            similarities = []
            for j in range(i+1, len(all_proposed)):
                if j not in used_indices:
                    sim = np.dot(embeddings[i], embeddings[j]) / (
                        np.linalg.norm(embeddings[i]) * np.linalg.norm(embeddings[j])
                    )
                    similarities.append((j, sim))
            
            # Merge very similar clusters (similarity > 0.9)
            merged_children = cluster['suggested_children'].copy()
            for j, sim in similarities:
                if sim > 0.9:  # High threshold - only merge very similar clusters
                    merged_children.extend(all_proposed[j]['suggested_children'])
                    used_indices.add(j)
            
            cluster['suggested_children'] = list(set(merged_children))
            unique_clusters.append(cluster)
            used_indices.add(i)
        
        return unique_clusters
    
    def assign_children_to_parents(self, child_clusters: Dict[int, Dict],
                                 parent_clusters: List[Dict],
                                 child_embeddings: np.ndarray,
                                 child_ids: List[int]) -> Dict[int, List[int]]:
        """Assign each child cluster to best parent using Claude"""
        parent_assignments = {i: [] for i in range(len(parent_clusters))}
        
        # For each child, find best parent
        for child_idx, child_id in enumerate(child_ids):
            child = child_clusters[child_id]
            
            # Create prompt for assignment
            parents_text = ""
            for i, parent in enumerate(parent_clusters):
                parents_text += f"{i}: {parent['name']} - {parent['description']}\n"
            
            prompt = f"""Given this child cluster:
Name: {child['name']}
Description: {child['description']}

Which parent cluster does it best belong to?

Parent clusters:
{parents_text}

Respond with just the number (0-{len(parent_clusters)-1}) of the best parent cluster."""

            try:
                response = self.client.messages.create(
                    model=self.claude_model,
                    max_tokens=10,
                    temperature=0.2,
                    messages=[{"role": "user", "content": prompt}]
                )
                
                parent_idx = int(response.content[0].text.strip())
                if 0 <= parent_idx < len(parent_clusters):
                    parent_assignments[parent_idx].append(child_id)
                else:
                    # Fallback to first parent
                    parent_assignments[0].append(child_id)
                    
            except:
                # Fallback: assign based on initial suggestions
                assigned = False
                for i, parent in enumerate(parent_clusters):
                    if child_idx + 1 in parent.get('suggested_children', []):
                        parent_assignments[i].append(child_id)
                        assigned = True
                        break
                if not assigned:
                    parent_assignments[0].append(child_id)
        
        return parent_assignments
    
    def regenerate_parent_names(self, parent_clusters: List[Dict],
                              parent_assignments: Dict[int, List[int]],
                              child_clusters: Dict[int, Dict]) -> List[Dict]:
        """Regenerate parent names based on actual assigned children"""
        updated_parents = []
        
        for parent_idx, parent in enumerate(parent_clusters):
            assigned_children = parent_assignments[parent_idx]
            
            if not assigned_children:
                continue
            
            # Get info about assigned children
            children_text = ""
            for child_id in assigned_children[:10]:  # Limit to 10 for prompt
                child = child_clusters[child_id]
                children_text += f"- {child['name']}: {child['description']}\n"
            
            prompt = f"""Based on these child clusters that have been grouped together:

{children_text}

Generate a parent cluster name and description that accurately represents this group.

Respond with JSON in this format:
{{"name": "Parent name (max 10 words)", "description": "Two sentence description."}}"""

            try:
                response = self.client.messages.create(
                    model=self.claude_model,
                    max_tokens=200,
                    temperature=1.0,
                    messages=[{"role": "user", "content": prompt}]
                )
                
                response_text = response.content[0].text.strip()
                if response_text.startswith('```json'):
                    response_text = response_text[7:]
                if response_text.endswith('```'):
                    response_text = response_text[:-3]
                
                parent_info = json.loads(response_text.strip())
                parent_info['children'] = assigned_children
                updated_parents.append(parent_info)
                
            except Exception as e:
                print(f"Error regenerating parent name: {e}")
                # Keep original
                parent['children'] = assigned_children
                updated_parents.append(parent)
        
        return updated_parents
    
    def build_hierarchy(self, base_clusters: Dict[int, Dict],
                       target_levels: int = 3,
                       top_k: int = 5) -> Dict:
        """Build complete hierarchy from base clusters to top level"""
        print(f"Building hierarchy (max {target_levels} levels)...")
        
        hierarchy = {
            'levels': [],
            'base_clusters': base_clusters
        }
        
        current_clusters = base_clusters
        current_level = 0
        
        # Continue building hierarchy until we reach target levels or can't reduce further
        while current_level < target_levels - 1 and len(current_clusters) > 1:
            print(f"\nLevel {current_level + 1}: {len(current_clusters)} clusters")
            
            # Embed current level
            embeddings, cluster_ids = self.embed_cluster_descriptions(current_clusters)
            
            # Calculate target number of parents
            n_current = len(current_clusters)
            levels_remaining = target_levels - current_level - 1
            
            # Each parent should have 5-10 children (sweet spot for exploration)
            min_children_per_parent = 5
            max_children_per_parent = 10
            
            # Calculate target based on desired children per parent
            min_parents = max(1, n_current // max_children_per_parent)
            max_parents = n_current // min_children_per_parent
            
            if n_current > 100:
                # For many clusters (e.g., 267), we need multiple levels
                # 267 -> 30-35 parents (8-9 children each)
                target_n_parents = max(min_parents, min(35, n_current // 8))
            elif n_current > 50:
                # For 50-100 clusters
                # E.g., 80 -> 10 parents (8 children each)
                target_n_parents = max(min_parents, min(10, n_current // 8))
            elif n_current > 30:
                # For 30-50 clusters
                # E.g., 40 -> 5 parents (8 children each)
                target_n_parents = max(min_parents, min(8, n_current // 7))
            elif n_current > 20:
                # For 20-30 clusters
                # E.g., 25 -> 4 parents (6-7 children each)
                target_n_parents = max(2, n_current // 7)
            elif n_current >= 10:
                # For 10-20 clusters
                # E.g., 15 -> 2 parents (7-8 children each)
                target_n_parents = max(2, n_current // 8)
            else:
                # For <10 clusters, create single top cluster
                target_n_parents = 1
            
            print(f"Target parent clusters: {target_n_parents}")
            
            # Group into neighborhoods
            # Create neighborhoods such that we don't exceed target parents
            # Each neighborhood will propose at least 1 parent
            target_per_neighborhood = max(1, n_current // target_n_parents)
            neighborhoods = self.group_into_neighborhoods(embeddings, target_per_neighborhood)
            n_neighborhoods = len(np.unique(neighborhoods))
            print(f"Created {n_neighborhoods} neighborhoods")
            
            # Propose parents for each neighborhood
            all_proposed = []
            for n in range(n_neighborhoods):
                mask = neighborhoods == n
                neighborhood_ids = [cluster_ids[i] for i in range(len(cluster_ids)) if mask[i]]
                neighborhood_clusters = [current_clusters[cid] for cid in neighborhood_ids]
                
                # Calculate proportional number of parents for this neighborhood
                n_in_neighborhood = len(neighborhood_clusters)
                n_parents_for_neighborhood = max(1, int(target_n_parents * n_in_neighborhood / n_current))
                
                proposed = self.propose_parent_clusters(neighborhood_clusters, n_parents_for_neighborhood)
                
                # Map suggested children back to actual cluster IDs
                for parent in proposed:
                    parent['suggested_children'] = [
                        neighborhood_ids[i-1] for i in parent.get('suggested_children', [])
                        if 0 < i <= len(neighborhood_ids)
                    ]
                
                all_proposed.extend(proposed)
            
            # Deduplicate
            parent_clusters = self.deduplicate_clusters(all_proposed)
            print(f"After deduplication: {len(parent_clusters)} parent clusters")
            
            # Assign children to parents
            parent_assignments = self.assign_children_to_parents(
                current_clusters, parent_clusters, embeddings, cluster_ids
            )
            
            # Regenerate parent names based on actual assignments
            parent_clusters = self.regenerate_parent_names(
                parent_clusters, parent_assignments, current_clusters
            )
            
            # Create parent cluster dict for next level
            next_level_clusters = {}
            for i, parent in enumerate(parent_clusters):
                next_level_clusters[i] = {
                    'name': parent['name'],
                    'description': parent['description'],
                    'children': parent['children']
                }
            
            # Store this level
            hierarchy['levels'].append({
                'level': current_level + 1,
                'clusters': next_level_clusters
            })
            
            # Move to next level
            current_clusters = next_level_clusters
            current_level += 1
        
        print(f"\nHierarchy complete with {len(hierarchy['levels'])} levels")
        return hierarchy
    
    def save_hierarchy(self, hierarchy: Dict, filepath: str = "data/hierarchy.json"):
        """Save hierarchy to file"""
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        # Convert numpy types to regular Python types for JSON serialization
        def convert_to_serializable(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {convert_to_serializable(k): convert_to_serializable(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_to_serializable(v) for v in obj]
            return obj
        
        serializable_hierarchy = convert_to_serializable(hierarchy)
        
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(serializable_hierarchy, f, indent=2, ensure_ascii=False)
        
        print(f"Saved hierarchy to {filepath}")


def test_hierarchizer():
    """Test the Hierarchizer module"""
    print("Testing Hierarchizer...")
    
    # Load cluster names from previous step
    with open('data/test_cluster_names.json', 'r') as f:
        cluster_names = json.load(f)
    
    # Convert string keys back to int
    cluster_names = {int(k): v for k, v in cluster_names.items()}
    
    # Initialize hierarchizer
    hierarchizer = Hierarchizer()
    
    # Build hierarchy (with just 2 levels for test data)
    hierarchy = hierarchizer.build_hierarchy(
        base_clusters=cluster_names,
        target_levels=2,
        top_k=2
    )
    
    # Save hierarchy
    hierarchizer.save_hierarchy(hierarchy, "data/test_hierarchy.json")
    
    # Print hierarchy structure
    print("\nHierarchy structure:")
    print(f"Base clusters: {len(hierarchy['base_clusters'])}")
    for level_info in hierarchy['levels']:
        level = level_info['level']
        clusters = level_info['clusters']
        print(f"Level {level}: {len(clusters)} clusters")
        for cid, cluster in clusters.items():
            print(f"  - {cluster['name']} ({len(cluster['children'])} children)")
    
    return True


if __name__ == "__main__":
    test_hierarchizer()
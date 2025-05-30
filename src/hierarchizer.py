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
    
    def _create_query_root_cluster(self, query: str, top_clusters: Dict[int, Dict]) -> Dict:
        """Create a root cluster based on the search query using Claude"""
        # Format top clusters for context
        clusters_text = ""
        for i, (cid, cluster) in enumerate(top_clusters.items()):
            clusters_text += f"{i+1}. {cluster['name']}: {cluster['description']}\n"
        
        prompt = f"""Given this PubMed search query used to collect biomedical research abstracts:
"{query}"

And these are the top-level clusters that were discovered from the PubMed abstracts:
{clusters_text}

Generate a name and description for the root cluster that encompasses all of these research areas.

Requirements:
1. The name should be descriptive and concise (max 10 words)
2. The description should be 2 sentences that capture the breadth of the collection
3. Both should reflect the original search intent while encompassing all discovered themes
4. Consider that this is a collection of biomedical/scientific literature from PubMed

Respond with JSON in this format:
{{"name": "Root cluster name", "description": "Two sentence description."}}"""

        try:
            response = self.client.messages.create(
                model=self.claude_model,
                max_tokens=200,
                temperature=0.7,
                messages=[{"role": "user", "content": prompt}]
            )
            
            response_text = response.content[0].text.strip()
            if response_text.startswith('```json'):
                response_text = response_text[7:]
            if response_text.endswith('```'):
                response_text = response_text[:-3]
            
            root_info = json.loads(response_text.strip())
            
            # Add children and root flag
            root_info['children'] = list(top_clusters.keys())
            root_info['is_root'] = True
            
            return root_info
            
        except Exception as e:
            print(f"Error generating root cluster name: {e}")
            # Re-raise the error - don't generate generic names
            raise RuntimeError(f"Failed to generate root cluster name: {e}")
    
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
    
    def balance_parent_children(self, parent_assignments: Dict[int, List[int]], 
                               parent_clusters: List[Dict],
                               child_clusters: Dict[int, Dict],
                               min_children: int = 2, 
                               max_children: int = 9,
                               ideal_children: int = 5) -> Tuple[Dict[int, List[int]], List[Dict]]:
        """Redistribute children to ensure balanced parent clusters using Claude"""
        print(f"\nBalancing parent-child relationships (target: {min_children}-{max_children} children per parent)")
        
        # Identify parents that need rebalancing
        overloaded = {p: children for p, children in parent_assignments.items() 
                      if len(children) > max_children}
        underloaded = {p: children for p, children in parent_assignments.items() 
                       if 0 < len(children) < min_children}
        empty_parents = [p for p, children in parent_assignments.items() if len(children) == 0]
        
        print(f"Found {len(overloaded)} overloaded parents, {len(underloaded)} underloaded parents")
        
        # Remove empty parents
        parent_clusters = [p for i, p in enumerate(parent_clusters) if i not in empty_parents]
        new_assignments = {i: parent_assignments[old_i] 
                          for i, old_i in enumerate(sorted(set(range(len(parent_assignments))) - set(empty_parents)))}
        
        # Split overloaded parents
        new_parent_id = len(parent_clusters)
        for parent_id, children in overloaded.items():
            if parent_id in empty_parents:
                continue
                
            # Map to new assignment index
            new_parent_idx = list(new_assignments.keys())[list(new_assignments.values()).index(children)]
            
            while len(new_assignments[new_parent_idx]) > max_children:
                # Take some children for a new parent
                split_size = min(ideal_children, len(new_assignments[new_parent_idx]) - ideal_children)
                split_children = new_assignments[new_parent_idx][-split_size:]
                new_assignments[new_parent_idx] = new_assignments[new_parent_idx][:-split_size]
                
                # Create new parent description using Claude
                children_text = ""
                for child_id in split_children[:10]:
                    child = child_clusters[child_id]
                    children_text += f"- {child['name']}: {child['description']}\n"
                
                prompt = f"""These child clusters need a parent cluster:

{children_text}

Generate a parent cluster name and description.

Respond with JSON in this format:
{{"name": "Parent name (max 10 words)", "description": "Two sentence description."}}"""

                try:
                    response = self.client.messages.create(
                        model=self.claude_model,
                        max_tokens=200,
                        temperature=0.7,
                        messages=[{"role": "user", "content": prompt}]
                    )
                    
                    response_text = response.content[0].text.strip()
                    if response_text.startswith('```json'):
                        response_text = response_text[7:]
                    if response_text.endswith('```'):
                        response_text = response_text[:-3]
                    
                    new_parent = json.loads(response_text.strip())
                    new_parent['children'] = split_children
                    parent_clusters.append(new_parent)
                    new_assignments[new_parent_id] = split_children
                    new_parent_id += 1
                    
                except Exception as e:
                    print(f"Error creating split parent: {e}")
                    # Fallback: just keep the split without creating new parent
                    break
        
        # Handle underloaded parents by merging or redistributing
        for parent_idx in list(underloaded.keys()):
            if parent_idx in empty_parents or parent_idx not in new_assignments:
                continue
                
            children = new_assignments.get(parent_idx, [])
            if len(children) >= min_children:  # Already fixed
                continue
                
            # Find a suitable parent to merge with (prefer smaller parents)
            candidates = [(p, len(cs)) for p, cs in new_assignments.items() 
                         if p != parent_idx and len(cs) + len(children) <= max_children]
            
            if candidates:
                # Merge with smallest suitable parent
                merge_target = min(candidates, key=lambda x: x[1])[0]
                new_assignments[merge_target].extend(children)
                new_assignments[parent_idx] = []
        
        # Remove empty assignments and renumber
        final_assignments = {}
        final_parents = []
        new_idx = 0
        
        for old_idx, (parent_idx, children) in enumerate(new_assignments.items()):
            if children:  # Only keep non-empty parents
                final_assignments[new_idx] = children
                if parent_idx < len(parent_clusters):
                    final_parents.append(parent_clusters[parent_idx])
                new_idx += 1
        
        print(f"After balancing: {len(final_parents)} parents")
        children_counts = [len(children) for children in final_assignments.values()]
        print(f"Children per parent: min={min(children_counts)}, max={max(children_counts)}, avg={sum(children_counts)/len(children_counts):.1f}")
        
        return final_assignments, final_parents
    
    def collapse_single_child_clusters(self, hierarchy: Dict) -> None:
        """Collapse clusters that have only one child by promoting grandchildren to direct children"""
        # Process each level from top to bottom
        for level_idx in range(len(hierarchy['levels']) - 1, -1, -1):
            level = hierarchy['levels'][level_idx]
            clusters = level['clusters']
            clusters_to_remove = []
            
            for cluster_id, cluster in list(clusters.items()):
                if len(cluster.get('children', [])) == 1:
                    single_child_id = cluster['children'][0]
                    
                    # Find the child cluster (could be in next level or base clusters)
                    child_cluster = None
                    child_level_idx = None
                    
                    # Check next level first
                    if level_idx > 0:
                        next_level = hierarchy['levels'][level_idx - 1]
                        if str(single_child_id) in next_level['clusters']:
                            child_cluster = next_level['clusters'][str(single_child_id)]
                            child_level_idx = level_idx - 1
                    
                    # If not found and we're at level 1, check base clusters
                    if child_cluster is None and level_idx == 0:
                        if single_child_id in hierarchy['base_clusters']:
                            child_cluster = hierarchy['base_clusters'][single_child_id]
                    
                    if child_cluster and child_level_idx is not None:
                        # Get grandchildren
                        grandchildren = child_cluster.get('children', [])
                        
                        if grandchildren:
                            print(f"Collapsing single-child cluster: '{cluster['name']}' -> '{child_cluster['name']}' -> {len(grandchildren)} grandchildren")
                            
                            # Promote grandchildren to be direct children of parent
                            cluster['children'] = grandchildren
                            
                            # Remove the intermediate child cluster
                            del hierarchy['levels'][child_level_idx]['clusters'][str(single_child_id)]
                    elif child_cluster is None and level_idx == 0:
                        # This is a parent with a single base cluster child
                        # In this case, we should remove this parent entirely
                        clusters_to_remove.append(cluster_id)
                        print(f"Removing single-child parent cluster: '{cluster['name']}' with base cluster child {single_child_id}")
            
            # Remove clusters marked for removal
            for cluster_id in clusters_to_remove:
                del clusters[cluster_id]
            
            # Update parent references in previous level if needed
            if level_idx < len(hierarchy['levels']) - 1:
                parent_level = hierarchy['levels'][level_idx + 1]
                for parent_cluster in parent_level['clusters'].values():
                    # Remove references to removed clusters
                    parent_cluster['children'] = [
                        child for child in parent_cluster.get('children', [])
                        if child not in clusters_to_remove
                    ]
    
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
                       query: str = None,
                       enforce_children_constraints: bool = True,
                       min_children: int = 2,
                       max_children: int = 9,
                       ideal_children: int = 5) -> Dict:
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
            
            # Simple calculation: target ideal_children per parent
            target_n_parents = max(1, n_current // ideal_children)
            
            # Ensure we don't create too few or too many parents
            min_allowed_parents = max(1, n_current // max_children)  # Each parent has at most max_children
            max_allowed_parents = n_current // min_children  # Each parent has at least min_children
            
            # Apply constraints
            target_n_parents = max(min_allowed_parents, min(max_allowed_parents, target_n_parents))
            
            # Special case: if we would create just 1-2 parents and we're not at the last level
            if target_n_parents <= 2 and levels_remaining > 1 and n_current > 10:
                # Adjust to create more parents for better hierarchy
                target_n_parents = max(3, n_current // 7)
            
            print(f"Target parent clusters: {target_n_parents}")
            
            # Stop building hierarchy if we would only create 1 parent cluster
            if target_n_parents == 1:
                print(f"Stopping hierarchy building - would only create 1 parent cluster")
                break
            
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
            
            # Balance parent-child relationships if requested
            if enforce_children_constraints:
                parent_assignments, parent_clusters = self.balance_parent_children(
                    parent_assignments, parent_clusters, current_clusters,
                    min_children=min_children, max_children=max_children, ideal_children=ideal_children
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
        
        # Collapse any single-child clusters
        if hierarchy['levels']:
            print("\nChecking for single-child clusters to collapse...")
            self.collapse_single_child_clusters(hierarchy)
        
        # Add query-based root cluster if query is provided
        if query and hierarchy['levels']:
            # Get the top level clusters
            top_level = hierarchy['levels'][-1]
            top_clusters = top_level['clusters']
            
            # Create root cluster based on query
            root_cluster = self._create_query_root_cluster(query, top_clusters)
            
            # Add as the final level
            hierarchy['levels'].append({
                'level': len(hierarchy['levels']) + 1,
                'clusters': {0: root_cluster}
            })
            
            print(f"Added query-based root cluster: {root_cluster['name']}")
        
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
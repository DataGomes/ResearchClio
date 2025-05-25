# CLIO Implementation for PubMed Abstract Clustering

## Overview
Implementation of CLIO clustering system for PubMed abstracts about artificial intelligence, based on the technical details from the CLIO paper.

## Architecture Components

### 1. Data Collection Module
- **API**: PubMed E-utilities API
- **Query**: "artificial intelligence"
- **Fields to extract**:
  - PMID (unique identifier)
  - Title
  - Abstract text
  - Publication date
  - Authors
  - Journal

### 2. Embedding Module
- **Model**: all-mpnet-base-v2 (via sentence-transformers)
- **Output**: 768-dimensional embeddings
- **Input**: Combined title + abstract text

### 3. Clustering Module
- **Algorithm**: k-means (scikit-learn implementation)
- **k selection**: 
  - For initial clustering: sqrt(n_abstracts) as starting point
  - Adjust to ensure meaningful cluster sizes (min 5-10 abstracts per cluster)
  - No privacy constraints needed (public data)

### 4. Cluster Naming Module
- **Model**: Claude 3.5 Haiku (for cost efficiency)
- **Temperature**: 1.0 for creative naming
- **Sampling strategy**:
  - Sample up to 10 abstracts from cluster
  - Sample up to 10 nearest non-cluster abstracts
- **Output**: 
  - Cluster name (≤10 words)
  - 2-sentence description

### 5. Hierarchizer Module
Based on CLIO's hierarchical clustering approach:

**Algorithm**:
1. **Embed cluster descriptions** using all-mpnet-base-v2
2. **Group into neighborhoods** via k-means (target ~40 clusters per neighborhood)
3. **For each neighborhood**:
   - Use Claude to propose higher-level groupings
   - Target hierarchy ratios: n_l/n_(l-1) based on desired levels
4. **Deduplicate** proposed clusters across neighborhoods
5. **Assign children** to parents using Claude's judgment
6. **Regenerate** parent names based on actual children
7. **Repeat** until reaching desired top level (e.g., 5-10 top themes)

### 6. Output Module
Hierarchical JSON structure:
```json
{
  "top_level": [
    {
      "name": "Machine Learning in Healthcare",
      "description": "...",
      "children": [
        {
          "name": "Diagnostic AI Systems",
          "description": "...",
          "abstracts": [
            {"pmid": "...", "title": "...", "abstract": "..."}
          ]
        }
      ]
    }
  ]
}
```

## Implementation Steps

### Phase 1: Data Collection
```python
# PubMed API wrapper
def fetch_pubmed_abstracts(query, max_results=1000):
    # Use Entrez E-utilities
    # Return list of {pmid, title, abstract}
```

### Phase 2: Embedding Generation
```python
from sentence_transformers import SentenceTransformer

model = SentenceTransformer('all-mpnet-base-v2')
embeddings = model.encode([f"{title} {abstract}" for title, abstract in data])
```

### Phase 3: Initial Clustering
```python
from sklearn.cluster import KMeans

# Dynamic k selection
k = int(np.sqrt(len(abstracts)))
k = max(k, 20)  # Minimum 20 clusters
k = min(k, len(abstracts) // 10)  # Max clusters = n/10

kmeans = KMeans(n_clusters=k, random_state=42)
clusters = kmeans.fit_predict(embeddings)
```

### Phase 4: Cluster Naming with Claude
```python
def generate_cluster_name(cluster_abstracts, non_cluster_abstracts):
    prompt = f"""
    I have a cluster of scientific abstracts about AI. 
    
    Cluster abstracts (sample):
    {format_abstracts(cluster_abstracts[:10])}
    
    Nearby non-cluster abstracts:
    {format_abstracts(non_cluster_abstracts[:10])}
    
    Generate:
    1. A descriptive name (≤10 words) that captures what makes this cluster unique
    2. A 2-sentence description of the cluster's focus
    
    Format as JSON: {{"name": "...", "description": "..."}}
    """
    # Call Claude API
```

### Phase 5: Hierarchical Clustering
```python
def build_hierarchy(clusters, target_levels=3):
    current_level = clusters
    hierarchy = []
    
    for level in range(target_levels - 1):
        # Embed current level descriptions
        embeddings = model.encode([c['description'] for c in current_level])
        
        # Group into neighborhoods
        n_neighborhoods = max(1, len(current_level) // 40)
        neighborhoods = KMeans(n_clusters=n_neighborhoods).fit_predict(embeddings)
        
        # For each neighborhood, use Claude to propose parents
        next_level = []
        for n in range(n_neighborhoods):
            neighborhood_clusters = [c for i, c in enumerate(current_level) if neighborhoods[i] == n]
            parents = propose_parent_clusters(neighborhood_clusters)
            next_level.extend(parents)
        
        # Deduplicate and assign children
        next_level = deduplicate_clusters(next_level)
        assign_children_to_parents(current_level, next_level)
        
        current_level = next_level
    
    return build_tree(clusters, hierarchy)
```

## Key Differences from Original CLIO

1. **No Privacy Mechanisms**: Working with public abstracts
2. **No UMAP Visualization**: Focus on hierarchical structure only
3. **Simplified Facet Extraction**: Direct to clustering (no conversation analysis)
4. **Domain-Specific**: Optimized for scientific abstracts vs conversations

## Cost Estimation

For 1000 PubMed abstracts:
- Embeddings: Free (local model)
- Initial cluster naming: ~$0.50 (Haiku)
- Hierarchy generation: ~$0.10 (Haiku)
- Total: ~$0.60

## Dependencies

```python
sentence-transformers
scikit-learn
numpy
pandas
anthropic
biopython  # For PubMed API
```

## Output Example

```
Top Level: AI in Medicine
├── Diagnostic Applications
│   ├── Medical Imaging AI
│   │   ├── Abstract: "Deep learning for chest X-ray..." (PMID: 12345)
│   │   └── Abstract: "CNN-based mammography screening..." (PMID: 12346)
│   └── Clinical Decision Support
│       └── Abstract: "AI-assisted diagnosis system..." (PMID: 12347)
└── Treatment Planning
    └── Personalized Medicine AI
        └── Abstract: "Machine learning for drug selection..." (PMID: 12348)
```
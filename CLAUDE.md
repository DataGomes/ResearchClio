# CLIO Implementation for PubMed Abstract Clustering

## Project Overview
Complete implementation of the CLIO clustering system for PubMed abstracts about artificial intelligence, based on the CLIO paper methodology.

## Virtual Environment Setup
**IMPORTANT**: Always use the virtual environment for this project:
```bash
# Activate virtual environment
source venv/bin/activate

# Deactivate when done
deactivate
```

## Implementation Status ✅
All components have been successfully implemented and tested:

1. **PubMed Data Collection** (`src/pubmed_collector.py`)
   - Fetches abstracts using PubMed E-utilities API
   - Handles rate limiting and SSL issues
   - Extracts title, abstract, PMID, authors, journal
   - **Caching**: Automatically caches abstracts for 7 days to avoid redundant API calls

2. **Embeddings** (`src/embedder.py`)
   - Uses Sentence-BERT (all-mpnet-base-v2) 
   - Generates 768-dimensional embeddings
   - Includes similarity computation functions
   - **Offline Mode**: Works with cached model to avoid download timeouts

3. **Clustering** (`src/clusterer.py`)
   - K-means clustering with automatic k selection
   - Silhouette score optimization
   - Representative abstract identification

4. **Cluster Naming** (`src/cluster_namer.py`)
   - Claude API integration (Haiku model)
   - Samples in-cluster and out-cluster abstracts
   - Generates descriptive names and 2-sentence descriptions

5. **Hierarchizer** (`src/hierarchizer.py`)
   - Multi-level hierarchy building
   - Neighborhood-based processing
   - Parent cluster generation and deduplication

6. **Output Formatter** (`src/output_formatter.py`)
   - Text format (hierarchical tree view)
   - JSON format (structured data)
   - CSV format (flattened for analysis)

## Running the Pipeline

### Quick Test (20 abstracts)
```bash
source venv/bin/activate
export TOKENIZERS_PARALLELISM=false  # Suppress warnings
python src/run_clio_pipeline.py --max-abstracts 20 --min-cluster-size 3 --hierarchy-levels 2 --top-clusters 3
```

### Production Run (1000 abstracts)
```bash
source venv/bin/activate
export TOKENIZERS_PARALLELISM=false
python src/run_clio_pipeline.py --max-abstracts 1000 --min-cluster-size 10 --hierarchy-levels 3 --top-clusters 5
```

### Command Line Options
- `--max-abstracts`: Number of abstracts to fetch (default: 100)
- `--min-cluster-size`: Minimum abstracts per cluster (default: 5)
- `--hierarchy-levels`: Number of hierarchy levels (default: 3)
- `--top-clusters`: Number of top-level clusters (default: 5)
- `--output-dir`: Output directory (default: output)
- `--no-cache`: Disable caching of PubMed abstracts (cache enabled by default)

## Output Structure
Each run creates a timestamped directory with:
```
output/run_YYYYMMDD_HHMMSS/
├── abstracts.json              # Raw PubMed data
├── embeddings.pkl              # Generated embeddings
├── clustering.pkl              # Clustering results
├── cluster_names.json          # Claude-generated names
├── complete_hierarchy.json     # Full hierarchy with top cluster and all levels
├── hierarchy.json              # Detailed structure with abstract IDs
├── hierarchy.txt               # Human-readable tree view
└── clusters.csv               # Flattened CSV for analysis
```

### Key Output: complete_hierarchy.json
This is the main output for downstream applications, containing:
- Top-level cluster(s) with name and description
- Intermediate hierarchy levels (if any)
- All base clusters sorted by size
- Complete PMID lists for each cluster

## Cost Estimation
- 100 abstracts: ~$0.10 (using Haiku)
- 1000 abstracts: ~$1.00
- 10,000 abstracts: ~$10.00

## Environment Variables
Set in `.env` file:
- `ANTHROPIC_API_KEY`: Your Claude API key
- `PUBMED_EMAIL`: Your email for PubMed API
- `CLAUDE_MODEL`: Model choice (default: claude-3-haiku-20240307)

## Troubleshooting
1. **SSL Certificate Errors**: Already handled in code
2. **Rate Limiting**: Code includes delays between API calls
3. **Memory Issues**: Process abstracts in smaller batches
4. **JSON Serialization**: Numpy types are converted automatically

## Example Output
```
TOP LEVEL (3 clusters):

▪ Artificial Intelligence in Medical Imaging and Diagnostics
  This cluster focuses on AI applications in radiology and pathology...
  ▪ Deep Learning for Medical Images
    Applications of CNNs in X-ray and MRI analysis...
  ▪ AI in Cancer Detection
    Machine learning for early cancer diagnosis...

▪ AI Applications in Drug Discovery
  This cluster covers the use of AI in pharmaceutical research...
  ▪ Molecular Design with ML
    AI-driven drug molecule optimization...
  ▪ Clinical Trial Prediction
    Using AI to predict drug efficacy...
```

## Behavioral Guidelines
- If any part of the instructions is not clear, ask for clarifications
- It is acceptable and recommended to ask questions to ensure accurate task completion
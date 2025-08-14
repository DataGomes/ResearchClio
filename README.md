# CLIO for PubMed: Hierarchical Clustering of Scientific Abstracts

This project implements the CLIO hierarchical clustering system for PubMed abstracts, based on the approach described in the CLIO paper. It automatically organizes scientific literature about artificial intelligence into a meaningful hierarchical structure.

## Features

- **Automatic Abstract Collection**: Fetches AI-related abstracts from PubMed
- **Advanced Embeddings**: Uses Sentence-BERT (all-mpnet-base-v2) for semantic understanding
- **Smart Clustering**: K-means with automatic optimization
- **AI-Powered Naming**: Claude generates descriptive cluster names
- **Hierarchical Organization**: Multi-level clustering for better organization
- **Multiple Output Formats**: Text, JSON, and CSV outputs

## Quick Start

### 1. Setup Environment

```bash
# Clone the repository
cd clio

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### 2. Configure API Keys

Create a `.env` file:
```
ANTHROPIC_API_KEY=your_claude_api_key_here
PUBMED_EMAIL=your_email@example.com
CLAUDE_MODEL=claude-3-haiku-20240307
```

### 3. Run the Pipeline

```bash
# Test with 20 abstracts
python src/run_clio_pipeline.py --max-abstracts 20

# Production run with 1000 abstracts
python src/run_clio_pipeline.py --max-abstracts 1000 --hierarchy-levels 3
```

## How It Works

1. **Data Collection**: Queries PubMed for AI-related abstracts
2. **Embedding**: Converts abstracts to 768-dimensional vectors
3. **Clustering**: Groups similar abstracts using k-means
4. **Naming**: Claude analyzes clusters and generates names
5. **Hierarchy**: Builds multi-level structure
6. **Output**: Generates human-readable and machine-readable formats

## Example Output

```
TOP LEVEL (3 clusters):

▪ Artificial Intelligence in Medical Imaging and Diagnostics
  This cluster focuses on AI applications in radiology...
  ▪ Deep Learning for Medical Images
    CNNs for X-ray and MRI analysis...
  ▪ AI in Cancer Detection
    Early cancer diagnosis using ML...

▪ AI Applications in Drug Discovery
  AI in pharmaceutical research...
  ▪ Molecular Design with ML
    Drug molecule optimization...
```

## Command Line Options

- `--max-abstracts`: Number of abstracts to process (default: 100)
- `--min-cluster-size`: Minimum cluster size (default: 5)
- `--hierarchy-levels`: Hierarchy depth (default: 3)
- `--top-clusters`: Top-level clusters (default: 5)
- `--output-dir`: Output directory (default: output)

## Project Structure

```
clio/
├── src/
│   ├── pubmed_collector.py    # PubMed API integration
│   ├── embedder.py            # Sentence-BERT embeddings
│   ├── clusterer.py           # K-means clustering
│   ├── cluster_namer.py       # Claude API naming
│   ├── hierarchizer.py        # Hierarchy building
│   ├── output_formatter.py    # Output generation
│   └── run_clio_pipeline.py   # Main pipeline
├── output/                    # Generated outputs
├── requirements.txt           # Python dependencies
├── .env                       # API keys (create this)
└── README.md                  # This file
```

## Cost Estimation

Using Claude 3 Haiku:
- 100 abstracts: ~$0.10
- 1,000 abstracts: ~$1.00
- 10,000 abstracts: ~$10.00

## Citation

Based on the CLIO paper methodology. This implementation focuses on PubMed abstracts without privacy mechanisms (as abstracts are public data).

## License

MIT License

"""
PubMed Data Collection Module
Fetches abstracts about artificial intelligence from PubMed
"""

import os
import time
import json
import ssl
import certifi
from typing import List, Dict
from Bio import Entrez
from dotenv import load_dotenv
from tqdm import tqdm
import urllib.request
import hashlib
from datetime import datetime, timedelta

# Fix SSL certificate issues
ssl._create_default_https_context = ssl._create_unverified_context

load_dotenv()

class PubMedCollector:
    def __init__(self, email: str = None, cache_dir: str = "cache"):
        """Initialize PubMed collector with email for API access"""
        self.email = email or os.getenv('PUBMED_EMAIL')
        if self.email:
            Entrez.email = self.email
        else:
            print("Warning: No email provided. PubMed rate limits will be stricter.")
        
        # Setup cache directory
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)
    
    def _get_cache_key(self, query: str, max_results: int) -> str:
        """Generate a cache key based on query parameters"""
        cache_string = f"{query}:{max_results}"
        return hashlib.md5(cache_string.encode()).hexdigest()
    
    def _get_cache_path(self, cache_key: str) -> str:
        """Get the full path to the cache file"""
        return os.path.join(self.cache_dir, f"pubmed_{cache_key}.json")
    
    def _is_cache_valid(self, cache_path: str, max_age_days: int = 7) -> bool:
        """Check if cache file exists and is recent enough"""
        if not os.path.exists(cache_path):
            return False
        
        # Check age of cache file
        file_time = datetime.fromtimestamp(os.path.getmtime(cache_path))
        age = datetime.now() - file_time
        
        return age < timedelta(days=max_age_days)
    
    def _load_from_cache(self, cache_path: str) -> List[Dict]:
        """Load abstracts from cache file"""
        with open(cache_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        print(f"Loaded {len(data)} abstracts from cache")
        return data
    
    def _save_to_cache(self, abstracts: List[Dict], cache_path: str):
        """Save abstracts to cache file"""
        with open(cache_path, 'w', encoding='utf-8') as f:
            json.dump(abstracts, f, indent=2, ensure_ascii=False)
        print(f"Saved {len(abstracts)} abstracts to cache")
            
    def search_pubmed(self, query: str, max_results: int = 100) -> List[str]:
        """Search PubMed and return list of PMIDs"""
        print(f"Searching PubMed for: {query}")
        
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results,
            sort="relevance"
        )
        
        results = Entrez.read(handle)
        handle.close()
        
        pmids = results["IdList"]
        print(f"Found {len(pmids)} articles")
        return pmids
    
    def fetch_abstracts(self, pmids: List[str]) -> List[Dict]:
        """Fetch full abstract data for given PMIDs"""
        abstracts = []
        
        print(f"Fetching {len(pmids)} abstracts...")
        
        # Fetch in batches of 10 to respect rate limits
        batch_size = 10
        
        for i in tqdm(range(0, len(pmids), batch_size)):
            batch = pmids[i:i + batch_size]
            
            try:
                handle = Entrez.efetch(
                    db="pubmed",
                    id=",".join(batch),
                    rettype="xml",
                    retmode="xml"
                )
                
                records = Entrez.read(handle)
                handle.close()
                
                for article in records['PubmedArticle']:
                    abstract_data = self._extract_abstract_data(article)
                    if abstract_data:
                        abstracts.append(abstract_data)
                
                # Be nice to NCBI servers
                time.sleep(0.5)
                
            except Exception as e:
                print(f"Error fetching batch {i//batch_size + 1}: {e}")
                continue
                
        print(f"Successfully fetched {len(abstracts)} abstracts with text")
        return abstracts
    
    def _extract_abstract_data(self, article: Dict) -> Dict:
        """Extract relevant data from PubMed article"""
        try:
            medline = article['MedlineCitation']
            article_info = medline['Article']
            
            # Get PMID
            pmid = str(medline['PMID'])
            
            # Get title
            title = article_info.get('ArticleTitle', '')
            
            # Get abstract text
            abstract_text = ""
            if 'Abstract' in article_info:
                abstract_parts = article_info['Abstract'].get('AbstractText', [])
                if isinstance(abstract_parts, list):
                    abstract_text = ' '.join(str(part) for part in abstract_parts)
                else:
                    abstract_text = str(abstract_parts)
            
            # Skip if no abstract
            if not abstract_text:
                return None
            
            # Get publication date
            pub_date = "Unknown"
            if 'DateCompleted' in medline:
                date = medline['DateCompleted']
                pub_date = f"{date.get('Year', '')}-{date.get('Month', '')}-{date.get('Day', '')}"
            
            # Get journal
            journal = article_info.get('Journal', {}).get('Title', 'Unknown')
            
            # Get authors
            authors = []
            if 'AuthorList' in article_info:
                for author in article_info['AuthorList']:
                    if isinstance(author, dict):
                        last = author.get('LastName', '')
                        first = author.get('ForeName', '')
                        if last:
                            authors.append(f"{last}, {first}".strip(', '))
            
            return {
                'pmid': pmid,
                'title': title,
                'abstract': abstract_text,
                'journal': journal,
                'pub_date': pub_date,
                'authors': authors[:5]  # Limit to first 5 authors
            }
            
        except Exception as e:
            print(f"Error extracting data: {e}")
            return None
    
    def collect_ai_abstracts(self, max_results: int = 100, use_cache: bool = True) -> List[Dict]:
        """Main method to collect AI-related abstracts"""
        # Search query focused on artificial intelligence
        query = '("artificial intelligence" OR "machine learning" OR "deep learning") AND hasabstract'
        
        # Check cache first
        if use_cache:
            cache_key = self._get_cache_key(query, max_results)
            cache_path = self._get_cache_path(cache_key)
            
            if self._is_cache_valid(cache_path):
                print(f"Using cached PubMed data (cache key: {cache_key})")
                return self._load_from_cache(cache_path)
            else:
                print("Cache not found or expired, fetching from PubMed...")
        
        # Get PMIDs
        pmids = self.search_pubmed(query, max_results)
        
        # Fetch full abstracts
        abstracts = self.fetch_abstracts(pmids)
        
        # Save to cache if caching is enabled
        if use_cache and abstracts:
            cache_key = self._get_cache_key(query, max_results)
            cache_path = self._get_cache_path(cache_key)
            self._save_to_cache(abstracts, cache_path)
        
        return abstracts
    
    def save_abstracts(self, abstracts: List[Dict], filename: str = "pubmed_abstracts.json"):
        """Save abstracts to JSON file"""
        # If filename contains path, use it directly
        if os.path.dirname(filename):
            output_path = filename
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
        else:
            output_path = os.path.join("data", filename)
            os.makedirs("data", exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(abstracts, f, indent=2, ensure_ascii=False)
        
        print(f"Saved {len(abstracts)} abstracts to {output_path}")
        return output_path


def test_collector():
    """Test the PubMed collector with a small sample"""
    print("Testing PubMed Collector...")
    
    collector = PubMedCollector()
    
    # Test with just 10 abstracts
    abstracts = collector.collect_ai_abstracts(max_results=10)
    
    if abstracts:
        print(f"\nSample abstract:")
        print(f"Title: {abstracts[0]['title']}")
        print(f"PMID: {abstracts[0]['pmid']}")
        print(f"Abstract: {abstracts[0]['abstract'][:200]}...")
        print(f"Journal: {abstracts[0]['journal']}")
        print(f"Authors: {', '.join(abstracts[0]['authors'])}")
        
        # Save test data
        collector.save_abstracts(abstracts, "test_abstracts.json")
        return True
    else:
        print("No abstracts collected!")
        return False


if __name__ == "__main__":
    test_collector()
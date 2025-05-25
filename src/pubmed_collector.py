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
    
    def _get_cache_key(self, query: str, max_results: int = None, year_month: str = None) -> str:
        """Generate a cache key based on query parameters"""
        if year_month:
            # For month-specific caching
            cache_string = f"{query}:{year_month}"
        else:
            # For full query caching
            cache_string = f"{query}:{max_results}"
        return hashlib.md5(cache_string.encode()).hexdigest()
    
    def _get_segment_cache_key(self, query: str) -> str:
        """Generate a cache key for a specific query segment (month/day)"""
        return hashlib.md5(query.encode()).hexdigest()
    
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
    
    def _fetch_segment_with_cache(self, query: str, max_results: int = 9995, use_cache: bool = True) -> List[Dict]:
        """Fetch abstracts for a specific query segment with caching"""
        if use_cache:
            cache_key = self._get_segment_cache_key(query)
            cache_path = self._get_cache_path(cache_key)
            
            if self._is_cache_valid(cache_path):
                print(f"  Using cached data for this segment")
                return self._load_from_cache(cache_path)
        
        # Not in cache, fetch from PubMed
        pmids = self.search_pubmed(query, max_results=max_results)
        if not pmids:
            return []
            
        abstracts = self.fetch_abstracts(pmids)
        
        # Cache this segment
        if use_cache and abstracts:
            cache_key = self._get_segment_cache_key(query)
            cache_path = self._get_cache_path(cache_key)
            self._save_to_cache(abstracts, cache_path)
        
        return abstracts
            
    def search_pubmed(self, query: str, max_results: int = 100) -> List[str]:
        """Search PubMed and return list of PMIDs"""
        print(f"Searching PubMed for: {query}")
        
        # First, get the total count
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=0,  # Just get count
            sort="relevance"
        )
        results = Entrez.read(handle)
        handle.close()
        
        total_count = int(results["Count"])
        print(f"Total articles available: {total_count}")
        
        # If we want more than 9999, we need to handle it differently
        all_pmids = []
        batch_size = 9995  # Slightly less than PubMed's limit of 9999 to be safe
        
        # PubMed ESearch can only retrieve first 9999 records total
        if max_results > 9999 or total_count > 9999:
            print(f"Warning: PubMed ESearch limited to first 9999 results.")
            print(f"Requested: {max_results}, Available: {total_count}")
            actual_max = min(9995, max_results)
        else:
            actual_max = min(max_results, total_count)
        
        # Just get what we can in one request (max 9995)
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=actual_max,
            sort="relevance"
        )
        
        results = Entrez.read(handle)
        handle.close()
        
        pmids = results["IdList"]
        print(f"Retrieved {len(pmids)} articles")
        
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
    
    def collect_abstracts_by_month(self, base_query: str, year: int, language: str = "eng", use_cache: bool = True) -> List[Dict]:
        """Collect abstracts month by month to bypass API limits"""
        all_abstracts = []
        
        for month in range(1, 13):
            month_str = f"{month:02d}"
            month_query = f'{base_query} AND ("{year}/{month_str}"[Date - Publication])'
            
            if language:
                month_query = f'{month_query} AND {language}[Language]'
            
            print(f"\nSearching {year}/{month_str}...")
            month_abstracts = self._fetch_segment_with_cache(month_query, max_results=9995, use_cache=use_cache)
            if month_abstracts:
                all_abstracts.extend(month_abstracts)
                print(f"Total collected so far: {len(all_abstracts)}")
            
            # Be nice to NCBI
            time.sleep(1)
        
        return all_abstracts
    
    def collect_abstracts(self, query: str, max_results: int = 100, use_cache: bool = True) -> List[Dict]:
        """Generic method to collect abstracts with custom query"""
        # Check if we need month-by-month collection
        # First, check how many results the query would return
        if max_results > 9995:
            # Test the query to see total count
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=0,
                sort="relevance"
            )
            results = Entrez.read(handle)
            handle.close()
            total_count = int(results["Count"])
            
            if total_count > 9995:
                print(f"Query has {total_count} results. Will search month by month to get all.")
                # Extract year from query if present
                import re
                year_match = re.search(r'"(\d{4})"\[Date - Publication\]', query)
                if year_match:
                    year = int(year_match.group(1))
                    # Remove year from base query to add month-specific dates
                    base_query = re.sub(r'AND \("\d{4}"\[Date - Publication\] : "\d{4}"\[Date - Publication\]\)', '', query).strip()
                    
                    # Check cache for full dataset
                    if use_cache:
                        cache_key = self._get_cache_key(query, max_results)
                        cache_path = self._get_cache_path(cache_key)
                        
                        if self._is_cache_valid(cache_path):
                            print(f"Using cached PubMed data (cache key: {cache_key})")
                            return self._load_from_cache(cache_path)
                    
                    # Collect by month
                    all_abstracts = []
                    for month in range(1, 13):
                        month_str = f"{month:02d}"
                        month_query = f'{base_query} AND ("{year}/{month_str}"[Date - Publication])'
                        
                        print(f"\nSearching {year}/{month_str}...")
                        
                        # First check how many results this month has
                        handle = Entrez.esearch(
                            db="pubmed",
                            term=month_query,
                            retmax=0,
                            sort="relevance"
                        )
                        month_results = Entrez.read(handle)
                        handle.close()
                        month_count = int(month_results["Count"])
                        
                        if month_count > 9995:
                            # Need to split by day ranges
                            print(f"Month {year}/{month_str} has {month_count} articles. Splitting by day ranges...")
                            
                            # Try first half and second half of month
                            for day_range in [(1, 15), (16, 31)]:
                                start_day, end_day = day_range
                                day_query = f'{base_query} AND ("{year}/{month_str}/{start_day:02d}"[Date - Publication] : "{year}/{month_str}/{end_day:02d}"[Date - Publication])'
                                
                                print(f"  Searching {year}/{month_str}/{start_day:02d}-{end_day:02d}...")
                                
                                # Check count for this range
                                handle = Entrez.esearch(
                                    db="pubmed",
                                    term=day_query,
                                    retmax=0,
                                    sort="relevance"
                                )
                                range_results = Entrez.read(handle)
                                handle.close()
                                range_count = int(range_results["Count"])
                                
                                if range_count > 9995:
                                    # Even half month is too much, go day by day
                                    print(f"    Range has {range_count} articles. Searching day by day...")
                                    for day in range(start_day, min(end_day + 1, 32)):
                                        day_query = f'{base_query} AND ("{year}/{month_str}/{day:02d}"[Date - Publication])'
                                        print(f"    Searching {year}/{month_str}/{day:02d}...")
                                        
                                        day_abstracts = self._fetch_segment_with_cache(day_query, max_results=9995, use_cache=use_cache)
                                        if day_abstracts:
                                            all_abstracts.extend(day_abstracts)
                                            print(f"    Total collected so far: {len(all_abstracts)}")
                                            
                                            if len(all_abstracts) >= max_results:
                                                print(f"Reached requested limit of {max_results}")
                                                return all_abstracts[:max_results]
                                        
                                        time.sleep(0.5)
                                else:
                                    # This range is OK, fetch it
                                    range_abstracts = self._fetch_segment_with_cache(day_query, max_results=9995, use_cache=use_cache)
                                    if range_abstracts:
                                        all_abstracts.extend(range_abstracts)
                                        print(f"  Total collected so far: {len(all_abstracts)}")
                                        
                                        if len(all_abstracts) >= max_results:
                                            print(f"Reached requested limit of {max_results}")
                                            return all_abstracts[:max_results]
                                    
                                    time.sleep(1)
                        else:
                            # Month has less than 9995, fetch normally
                            month_abstracts = self._fetch_segment_with_cache(month_query, max_results=9995, use_cache=use_cache)
                            if month_abstracts:
                                all_abstracts.extend(month_abstracts)
                                print(f"Total collected so far: {len(all_abstracts)}")
                                
                                if len(all_abstracts) >= max_results:
                                    print(f"Reached requested limit of {max_results}")
                                    return all_abstracts[:max_results]
                            
                            time.sleep(1)
                    
                    # Save to cache
                    if use_cache and all_abstracts:
                        self._save_to_cache(all_abstracts[:max_results], cache_path)
                    
                    return all_abstracts[:max_results]
        
        # Regular approach for queries under 9995 results
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
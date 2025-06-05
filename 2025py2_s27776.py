#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever
Extended script to retrieve and process genetic sequence records from NCBI GenBank.
"""

from Bio import Entrez, SeqIO
import time
import os
import ssl
import pandas as pd
import matplotlib.pyplot as plt

ssl._create_default_https_context = ssl._create_unverified_context

class NCBIRetriever:
    def __init__(self, email, api_key):
        self.email = email
        self.api_key = api_key
        self.min_len = None
        self.max_len = None
        self.organism_name = None
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid):
            """Search for all records associated with a taxonomic ID."""
            print(f"Searching for records with taxID: {taxid}")
            try:
                # Get taxonomy information first
                handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
                records = Entrez.read(handle)
                organism_name = records[0]["ScientificName"]
                print(f"Organism: {organism_name} (TaxID: {taxid})")
                # Search for records
                search_term = f"txid{taxid}[Organism]"
                handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
                search_results = Entrez.read(handle)
                count = int(search_results["Count"])
                if count == 0:
                    print(f"No records found for {organism_name}")
                    return None
                print(f"Found {count} records")
                # Store search results for later use
                self.webenv = search_results["WebEnv"]
                self.query_key = search_results["QueryKey"]
                self.count = count
                return count
            except Exception as e:
                print(f"Error searching TaxID {taxid}: {e}")
                return None

    def fetch_records(self, max_records=100):
        #  """Fetch a batch of records using the stored search results."""
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("Run search_taxid() first.")
            return []
        try:
            batch_size = min(max_records, 500)
            handle = Entrez.efetch(
                db="nucleotide", rettype="gb", retmode="text",
                retstart=0, retmax=batch_size,
                webenv=self.webenv, query_key=self.query_key
            )
            records = []
            for rec in SeqIO.parse(handle, "gb"):
                seqlen = len(rec.seq)
                if (self.min_len is None or seqlen >= self.min_len) and \
                   (self.max_len is None or seqlen <= self.max_len):
                    records.append(rec)
            return records
        except Exception as e:
            print(f"Error fetching records: {e}")
            return []
        
    def sequence_length_filtering(self, min_len=None, max_len=None):
        self.min_len = min_len
        self.max_len = max_len

    def generate_csv(self, records, filename):
        if not records:
            print("no records are written")
            return None
        data = {
            'accession number': [r.id for r in records],
            'sequence length': [len(r.seq) for r in records],
            'sequence description': [r.description for r in records]
        }
        df = pd.DataFrame(data).sort_values("length", ascending=False)
        df.to_csv(filename, index=False)
        print(f"saved CSV to {filename}")
        return df

    def visualize_data(self, df, filename):
        if df is None or df.empty:
            print("no data to visualize")
            return
        plt.figure(figsize=(10, 6))
        plt.plot(df['accession'], df['length'], marker='o')
        plt.xticks(rotation=90, fontsize=8)
        plt.ylabel('sequence Length')
        plt.xlabel('GenBank accession number')
        plt.title(f'sequence lengths for {self.organism_name}')
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()
        print(f"chart saved to {filename}")

def main():
    # Get user credentials
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")
    retriever = NCBIRetriever(email, api_key)

    min_len = input("Minimum sequence length (press Enter for none): ")
    max_len = input("Maximum sequence length (press Enter for none): ")
    retriever.sequence_length_filtering(
        int(min_len) if min_len else None,
        int(max_len) if max_len else None
    )

    taxid = input("Enter taxonomic ID (taxid) of the organism: ")
    count = retriever.search_taxid(taxid)
    if not count:
        print("No records found. Exiting.")
        return

    print("\nFetching sample records...")
    records = retriever.fetch_records(start=0, max_records=5)
    if not records:
        print("No matching records found after filtering.")
        return

    output_file = f"taxid_{taxid}_sample.gb"
    df = retriever.generate_csv(records, f"{output_file}_report.csv")
    retriever.visualize_data(df, f"{output_file}_plot.png")
    with open(f"{output_file}_records.gb", "w") as f:
        SeqIO.write(records, f, "genbank")
    print(f"GenBank records saved to {output_file}_records.gb")

if __name__ == "__main__":
    main()

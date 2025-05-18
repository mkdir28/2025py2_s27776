#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever
Basic script to connect to NCBI and retrieve genetic sequence records for a given taxonomic ID.
"""

from Bio import Entrez, SeqIO
import time
import os
import pandas as pd
import matplotlib.pyplot as plt


class NCBIRetriever:
    def __init__(self, email, api_key):
        """Initialize with NCBI credentials."""
        self.email = email
        self.api_key = api_key
        self.min_len = None
        self.max_len = None
        self.organism_name = None

        # Set up Entrez
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

    def fetch_records(self, start=0, max_records=10):
        """Fetch a batch of records using the stored search results."""
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []

        try:
            # Limit to prevent overloading
            batch_size = min(max_records, 500)

            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )

            # The raw GenBank record text
            records_text = handle.read()

            return records_text

        except Exception as e:
            print(f"Error fetching records: {e}")
            return ""

    def sequence_length_filter(self, min_len=None, max_len=None):
        self.max_len = max_len
        self.min_len = min_len

    def generate_csv(self, records, filename):
        if not records:
            print("No records to generate CSV report")
            return None

        # Create DataFrame
        data = {
            'Accession': [rec.id for rec in records],
            'Length': [len(rec.seq) for rec in records],
            'Description': [rec.description for rec in records]
        }

        df = pd.DataFrame(data).sort_values('Length', ascending=False)
        df.to_csv(filename, index=False)
        print(f"CSV report saved to {filename}")
        return df

    def visualize_data(self, df, filename):
        if df is None or df.empty:
            print("No data to visualize")
            return

        plt.figure(figsize=(10, 6))
        plt.plot(df['Accession'], df['Length'], marker='o')
        plt.xticks(rotation=90, fontsize=8)
        plt.ylabel('Sequence Length')
        plt.xlabel('GenBank Accession Number')
        plt.title(f'Sequence Lengths for {self.organism_name}')
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()
        print(f"Chart saved to {filename}")

def main():
    # Get user credentials
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")

    # Create retriever object
    retriever = NCBIRetriever(email, api_key)

    # Get sequence length filters
    min_len_input = input("Enter minimum sequence length (press Enter for none): ")
    max_len_input = input("Enter maximum sequence length (press Enter for none): ")

    min_len = int(min_len_input) if min_len_input else None
    max_len = int(max_len_input) if max_len_input else None
    retriever.sequence_length_filter(min_len, max_len)

    # Get taxid from user
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")

    # Search for records
    count = retriever.search_taxid(taxid)

    if not count:
        print("No records found. Exiting.")
        return

    # Fetch records
    print("\nFetching records...")
    records = retriever.fetch_records()
    if not records:
        print("No records to process after filtering. Exiting.")
        return

    # Generate outputs
    base_filename = f"taxid_{taxid}"
    df = retriever.generate_csv(records, f"{base_filename}_report.csv")
    retriever.visualize_data(df, f"{base_filename}_plot.png")

    # Save GenBank records
    with open(f"{base_filename}_records.gb", "w") as f:
        SeqIO.write(records, f, "genbank")
    print(f"GenBank records saved to {base_filename}_records.gb")

if __name__ == "__main__":
    main()


#!/usr/bin/env python3
"""
NCBI GenBank Data Retriever
Basic script to connect to NCBI and retrieve genetic sequence records for a given taxonomic ID.
"""

from Bio import Entrez, SeqIO
import time
import os
import pandas as pd
import csv
import matplotlib.pyplot as plt


class NCBIRetriever:
    def __init__(self, email, api_key):
        """Initialize with NCBI credentials."""
        self.email = email
        self.api_key = api_key
        self.max_len = None
        self.min_len = None

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

    def sequence_length_filter(self, max_len, min_len):
        self.max_len = max_len
        self.min_len = min_len

    def csv_report(self, records, filename):
        """Generate CSV report from records."""
        if not records:
            print("No records to generate CSV")
            return None

        data = {
            'Accession': [rec.id for rec in records],
            'Length': [len(rec.seq) for rec in records],
            'Description': [rec.description for rec in records]
        }
        df = pd.DataFrame(data)
        df.to_csv(filename, index=False)
        print(f"CSV report saved to {filename}")
        return df

    def data_visualisation(self, df, filename):
        """Generate and save a sequence length plot."""
        if df is None or df.empty:
            print("No data to plot")
            return

        plt.figure(figsize=(12, 6))
        plt.plot(df['Accession'], df['Length'], 'o-')
        plt.xticks(rotation=90)
        plt.xlabel('GenBank Accession Number')
        plt.ylabel('Sequence Length')
        plt.title(f'Sequence Lengths for {self.organism_name}')
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()
        print(f"Plot saved to {filename}")

def main():
    # Get user credentials
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")

    # Create retriever object
    retriever = NCBIRetriever(email, api_key)

        # Get sequence length filters
    min_len = input("Enter minimum sequence length (leave empty for no minimum): ")
    max_len = input("Enter maximum sequence length (leave empty for no maximum): ")
    if min_len:
        retriever.sequence_length_filter(int(max_len), int(min_len) if max_len else 1000000)

        # Get taxid from user
    taxid = input("Enter taxonomic ID (taxid) of the organism: ")

        # Search for records
    count = retriever.search_taxid(taxid)

    if not count:
        print("No records found. Exiting.")
        return

        # Fetch records
    print("\nFetching records...")
    records = retriever.fetch_records(start=0, max_records=min(100, count))

        # Generate outputs
    base_filename = f"taxid_{taxid}"
    df = retriever.csv_report(records, f"{base_filename}_report.csv")
    retriever.data_visualisation(df, f"{base_filename}_plot.png")

        # Save sample records
    output_file = f"{base_filename}_sample.gb"
    with open(output_file, "w") as f:
        SeqIO.write(records, f, "gb")
    print(f"Saved sample records to {output_file}")

if __name__ == "__main__":
    main()

#
# def main():
#     # Get user credentials
#     email = input("Enter your email address for NCBI: ")
#     api_key = input("Enter your NCBI API key: ")
#
#     # Create retriever object
#     retriever = NCBIRetriever(email, api_key)
#
#     max_len = input("Enter maximum sequence length (leave empty for no maximum): ")
#     min_len = input("Enter minimum sequence length (leave empty for no minimum): ")
#     if min_len:
#         retriever.sequence_length_filter(int(max_len), int(min_len) if max_len else 1000000)
#
#     # Get taxid from user
#     taxid = input("Enter taxonomic ID (taxid) of the organism: ")
#
#     # Search for records
#     count = retriever.search_taxid(taxid)
#
#     if not count:
#         print("No records found. Exiting.")
#         return
#
#     # Fetch first few records as a sample
#     print("\nFetching sample records...")
#     sample_records = retriever.fetch_records(start=0, max_records=5)
#
#     csv_file = f"taxid_{taxid}_filtered.csv"
#     img_file = f"taxid_{taxid}_plot.png"
#     df = csv_report(sample_records, csv_file)
#     make_plot(df, img_file)
#
#     print(f"CSV saved to {csv_file}")
#     print(f"Plot saved to {img_file}")
#     # Save to file
#     output_file = f"taxid_{taxid}_sample.gb"
#     with open(output_file, "w") as f:
#         f.write(sample_records)
#
#     print(f"Saved sample records to {output_file}")
#     print("\nNote: This is just a basic retriever. You need to extend its functionality!")
#
#
# if __name__ == "__main__":
#     main()


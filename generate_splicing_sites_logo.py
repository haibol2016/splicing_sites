#!/usr/bin/env python

import argparse
import os
from logomaker import Logo
import pandas as pd
import re
import gzip
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt


filter_biotype=[]
def get_exon_number_of_each_transcript(gtf_file):
    """
    Get the number of exons for each transcript from a GTF file.
    arguments:
    gtf_file -- a file object for the GTF file
    returns:
    A dictionary with transcript_id as keys and the number of exons as values.
    """
    exon_counts = {}
    for line in gtf_file:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] == "exon":
            attributes = fields[8]
            # Extract transcript_id from attributes
            transcript_id = None
            for attr in attributes.split(";"):
                if "transcript_id" in attr:
                    transcript_id = attr.split('"')[1]
                    break
            if transcript_id:
                if transcript_id not in exon_counts:
                    exon_counts[transcript_id] = 0
                exon_counts[transcript_id] += 1
    return exon_counts

def process_gtf(gtf_file):
    """
    Process a GTF file to extract unique exon information.
    exclude single-exon transcripts and return unique exons.
    arguments:
    gtf_file -- a file object for the GTF file
    returns:
    A list of unique exons in the format (chrom, start, end, strand).
    """

    exon_counts = get_exon_number_of_each_transcript(gtf_file)
    gtf_file.seek(0)  # Reset file pointer to the beginning after counting exons
    
    attr_p = re.compile(r';\s*')
    exons =[]
    for line in gtf_file:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] == "exon":
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            if strand != "+" and strand != "-":
                continue
            # Extract transcript_id from attributes
            transcript_id = None
            exon_number = None
            transcript_biotype=None
            for attr in re.split(attr_p, attributes):
                if "transcript_id" in attr:
                    transcript_id = attr.split('"')[1]
                if "exon_number" in attr:
                    # Ensembl GTF use exon_number "3", but Gencode GTF use exon_number 3.
                    if '"' in attr:
                        exon_number = int(attr.split('"')[1])
                    else:
                        exon_number = int(attr.split(' ')[1])
                if "transcript_biotype" in attr:
                    transcript_biotype = attr.split('"')[1]
                    if transcript_biotype in filter_biotype:
                        continue
                if transcript_id and exon_number and transcript_biotype:
                    break
            if not transcript_id or exon_counts[transcript_id] == 1:
                continue
            else:
                if exon_number == 1:
                    if strand == "+":
                        exons.append((chrom, None, end, strand)) # don't need the start of the first exon
                    else:
                        exons.append((chrom, start, None, strand))
                elif exon_number == exon_counts[transcript_id]: #don't need the end of the last exon
                    if strand == "+":
                        exons.append((chrom, start, None, strand))
                    else:
                        exons.append((chrom, None, end, strand))
                else:
                    exons.append((chrom, start, end, strand))
    unique_exons = set(exons)
    return list(unique_exons)

def parse_gtf(gtf_file):
    """Parse a GTF file to extract unique exon information. Exclude single-exon transcripts.
    arguments:
    gtf_file -- a GTF file path or file object
    returns:
    A list of unique exons in the format (chrom, start, end, strand).
    """
    exons =[]
    if gtf_file.endswith(".gz"):
        with gzip.open(gtf_file, "rt") as gtf:
            exons = process_gtf(gtf)
    else:
        with open(gtf_file, "r") as gtf:
            exons = process_gtf(gtf)
    
    # print(exons[0:20])
    return exons


def extract_flanking_sequences(exons, genome_fasta):
    """Extract flanking sequences for splicing donor and acceptor sites.
    arguments:
    exons -- a list of unique exons in the format (chrom, start, end, strand)
    genome_fasta -- a genome FASTA file path or file object
    returns:
    A dictionary with keys "donor" and "acceptor", each containing a list of flanking sequences.
    """
    if genome_fasta.endswith(".gz"):
        genome = SeqIO.to_dict(SeqIO.parse(gzip.open(genome_fasta, "rt"), "fasta"))
    else:
        genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    flanking_sequences = {"donor": [], "acceptor": []}
    
    # gtf is 1-based, biopython seq is 0-based, start and end should -1
    for chrom, start, end, strand in exons:
        if chrom not in genome:
            continue
        seq = genome[chrom].seq
        donor_site, acceptor_site = None, None
        if strand == "+":
            if start:
                acceptor_site = seq[(start - 21):(start + 2)]
            if end:
                donor_site = seq[(end - 3):(end + 9)]
        elif strand == "-":
            if start:
                donor_site = seq[(start - 10):(start + 2)].reverse_complement()

            if end:
                acceptor_site = seq[(end - 3):(end + 20)].reverse_complement()
        else:
            continue

        if donor_site:
            flanking_sequences["donor"].append(str(donor_site))
        if acceptor_site:
            flanking_sequences["acceptor"].append(str(acceptor_site))
    return flanking_sequences


def generate_sequence_logo(sequences, title, ax=None):
    """Generate a sequence logo from a list of sequences.
    arguments:
    sequences -- a list of sequences (strings)
    title -- title for the sequence logo
    ax -- matplotlib axes object to plot on (optional)
    returns:
    A sequence logo plot.
    """ 
    # Count the frequency of each nucleotide at each position
    counts = [Counter(pos) for pos in zip(*sequences)]
    df = pd.DataFrame(counts).fillna(0)
    df = df.div(df.sum(axis=1), axis=0)  # Normalize to probabilities
    # Create the sequence logo
    logo = Logo(df, shade_below=.5, fade_below=.5, ax=ax)
    logo.style_xticks(anchor=0, spacing=1)
    logo.ax.set_title(title)
    logo.ax.set_ylabel("Proportion")
    logo.ax.set_xlabel("Position")
    return logo


def generate_sequence_logo_for_splicing_sites(gtf_file, genome_fasta, filename='logo.pdf'):
    """Generate sequence logos for splicing donor and acceptor sites from a GTF file and a genome FASTA file.
    arguments:
    gtf_file -- a GTF file path or file object
    genome_fasta -- a genome FASTA file path or file object
    filename -- an output filename for the sequence logo (default: 'logo.pdf')
    returns:
    A sequence logo plot saved to a file.
    """
    # Parse the GTF file to get unique exons
    exons = parse_gtf(gtf_file)
    flanking_sequences = extract_flanking_sequences(exons, genome_fasta)

    # Generate sequence logos for donor and acceptor sites
    filename, file_extension = os.path.splitext(filename)
     
    # Generate sequence logos and save them to a single file
    plt.clf()
    # First logo 1/2 narrower (0.5x), second logo 2x wider (2x) -> ratio 0.5:2 = 1:4
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 2), gridspec_kw={'width_ratios': [1, 2]})
    plt.subplots_adjust(wspace=0.04)  # Adjust space between subplots

    # Generate donor logo on the first subplot
    donor_logo = generate_sequence_logo(flanking_sequences["donor"], "Donor Site Consensus Motif", ax=ax1)
    donor_logo.ax.set_xticklabels('%+d'%x for x in list(range(-3, 0)) + list(range(1, 10)))
    
    # Generate acceptor logo on the second subplot
    acceptor_logo = generate_sequence_logo(flanking_sequences["acceptor"], "Acceptor Site Consensus Motif", ax=ax2)
    acceptor_logo.ax.set_xticklabels('%+d'%x for x in list(range(-20, 0)) + list(range(1, 4)))
    # Remove y-axis label and ticks from the second plot only
    acceptor_logo.ax.set_ylabel("")
    acceptor_logo.ax.set_yticks([])
    acceptor_logo.ax.tick_params(left=False, labelleft=False)
    
    # Save both logos to a single file
    output_format = file_extension[1:] if file_extension else 'pdf'
    plt.savefig(filename + file_extension, format=output_format, bbox_inches='tight')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='generate consensus splicing site logos from a GTF file and a genome fasta file',
                                     usage='%(prog)s [-h] [--gtf GTF] [--genome-fasta GENOME_FASTA]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--gtf", required=True, metavar='',
                        help="an uncompressed or gzipped gtf file")
    parser.add_argument("--genome_fasta", required=True, metavar='',
                        help="an uncompressed or gzipped genome fasta file")
    parser.add_argument("--save_filename", metavar='', required =False,
                        help="an output filename for the sequence logo (default: logo.pdf), "
                        "the format is determined by the extension",
                        default='logo.pdf')

    args = parser.parse_args()
    print('Input gtf                  : ' + args.gtf)
    print('Input genome fasta         : ' + args.genome_fasta)
    print('Output filename            : ' + args.save_filename)

    gtf_file = args.gtf  # Replace with the actual GTF file path
    genome_fasta = args.genome_fasta  # Replace with the actual genome FASTA file path
    filename = args.save_filename  # Replace with the desired output filename
    if filename:
        generate_sequence_logo_for_splicing_sites(gtf_file, genome_fasta, filename=filename)
    else:
        generate_sequence_logo_for_splicing_sites(gtf_file, genome_fasta)


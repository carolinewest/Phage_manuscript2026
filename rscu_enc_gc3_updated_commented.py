
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ============================================================
# Calculate Relative Synonymous Codon Usage (RSCU)
# ============================================================
#
# PURPOSE:
# This script calculates RSCU values for every gene in a
# collection of fungal, bacterial, and phage coding sequences.
#
# RSCU measures how frequently a particular codon is used
# relative to the expected frequency if all synonymous codons
# (codons encoding the same amino acid) were used equally.
#
# RSCU interpretation:
#   RSCU = 1.0  : Codon is used as expected (no bias)
#   RSCU > 1.0  : Codon is used more frequently than expected
#   RSCU < 1.0  : Codon is used less frequently than expected
#   RSCU = 0.0  : Codon was not observed in the gene
#
# OUTPUT:
# Generates one CSV file per input folder.
# Each row represents an individual gene, and each codon
# has its own column containing its calculated RSCU value.
#
# GENETIC CODE:
# Uses the standard genetic code (NCBI translation table 1).
# Stop codons are excluded from RSCU calculations.
#
# ============================================================

import os
from Bio import SeqIO


# ============================================================
# 1. INPUT AND OUTPUT DIRECTORIES
# ============================================================

# List of directories containing FASTA files.
# Each directory represents a different dataset:
#
# fungal_cds:
#   Fungal coding DNA sequences.
#
# bacterial_ffn_split_contigs:
#   Bacterial coding sequences, separated by contig.
#
# phage_cds_cleaned_prokka:
#   Phage coding sequences predicted using Prokka.

input_folders = [
    "fungal_cds",
    "bacterial_ffn_split_contigs",
    "phage_cds_cleaned_prokka"
]

# Directory where all RSCU output files will be saved.
# Creates the directory if it does not already exist.

output_folder = "RSCU_results"
os.makedirs(output_folder, exist_ok=True)


# ============================================================
# 2. STANDARD GENETIC CODE
# ============================================================

# Dictionary mapping each DNA codon to its corresponding
# amino acid using the standard genetic code.
#
# Codons are represented using DNA nucleotides:
# A = Adenine
# T = Thymine
# C = Cytosine
# G = Guanine
#
# Amino acids are represented using their one-letter codes.
#
# "*" represents a stop codon.
#
# This dictionary is used to:
#   1. Identify the amino acid encoded by each codon.
#   2. Group synonymous codons together.
#   3. Calculate the expected frequency of each codon
#      based on the total usage of its amino acid.

aa_dict = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

# Create a sorted list of all sense codons (codons that encode
# amino acids).
#
# Stop codons are excluded because they are not used when
# calculating synonymous codon usage.
#
# The resulting list contains 61 codons.
#
# Sorting ensures that codon columns appear in a consistent
# order across all output files.

codon_list = sorted([c for c in aa_dict if aa_dict[c] != "*"])


# ============================================================
# 3. CALCULATE RSCU FOR AN INDIVIDUAL GENE
# ============================================================

def calc_rscu(seq):

    # Initialize the count of every sense codon to zero.
    # This ensures that codons not present in a gene are
    # still represented in the output.

    codon_counts = {c: 0 for c in codon_list}

    # Dictionary to store the total number of occurrences
    # of each amino acid across the gene.
    #
    # For example:
    #   F = total number of phenylalanine codons
    #   L = total number of leucine codons
    #
    # These totals are used to calculate expected codon usage.

    aa_totals = {}

    # Read the sequence in consecutive groups of three
    # nucleotides (codons).
    #
    # The range stops before the final incomplete codon,
    # if the sequence length is not divisible by three.

    for i in range(0, len(seq)-2, 3):

        # Extract the current three-nucleotide codon.

        codon = seq[i:i+3]

        # Only count codons present in the standard genetic
        # code dictionary of sense codons.
        #
        # Ambiguous or invalid codons and stop codons
        # are excluded from the calculation.

        if codon in codon_counts:

            # Add one to the observed count of this codon.

            codon_counts[codon] += 1

            # Identify the amino acid encoded by this codon.

            aa = aa_dict[codon]

            # Add one to the total count of this amino acid.
            # This combines all synonymous codons encoding
            # the same amino acid.

            aa_totals[aa] = aa_totals.get(aa, 0) + 1


    # --------------------------------------------------------
    # Calculate RSCU for each codon
    # --------------------------------------------------------

    rscu_vals = []

    # Iterate through all 61 sense codons in the same order
    # used for the output CSV columns.

    for codon in codon_list:

        # Identify the amino acid encoded by the current codon.

        aa = aa_dict[codon]

        # Identify all codons that encode the same amino acid.
        #
        # These are the synonymous codons used to calculate
        # the expected frequency of the current codon.

        syn_codons = [c for c in codon_list if aa_dict[c] == aa]

        # Calculate the expected count for this codon.
        #
        # Expected count =
        # Total occurrences of the amino acid /
        # Number of synonymous codons for that amino acid
        #
        # Assumes equal usage of all synonymous codons.

        expected = aa_totals.get(aa, 0) / len(syn_codons) if aa in aa_totals else 0

        # Calculate RSCU:
        #
        # RSCU = Observed codon count / Expected codon count
        #
        # If the amino acid does not occur in the gene,
        # its expected count is zero and RSCU is assigned 0.

        rscu = codon_counts[codon] / expected if expected > 0 else 0

        # Round the RSCU value to four decimal places
        # and append it to the output list.

        rscu_vals.append(round(rscu, 4))

    # Return the 61 RSCU values in the same order as codon_list.

    return rscu_vals


# ============================================================
# 4. PROCESS EACH INPUT DATASET
# ============================================================

# Process each input directory independently.
# Each directory will generate its own gene-level RSCU CSV.

for folder in input_folders:

    # Construct the output filename using the input folder name.
    #
    # Example:
    #   fungal_cds_gene_RSCU.csv

    gene_rscu_file = os.path.join(
        output_folder,
        f"{folder}_gene_RSCU.csv"
    )

    # --------------------------------------------------------
    # Create the output CSV and write the column headers.
    # --------------------------------------------------------

    # The first three columns identify each gene:
    #   file_name  = Name of the original FASTA file
    #   genome_id  = Genome identifier (currently file name)
    #   gene_id    = FASTA record identifier
    #
    # The remaining 61 columns contain the RSCU values
    # for each sense codon.

    with open(gene_rscu_file, "w") as out:

        out.write(
            "file_name,genome_id,gene_id," +
            ",".join(codon_list) + "\n"
        )


    # --------------------------------------------------------
    # Identify and process FASTA files in the directory.
    # --------------------------------------------------------

    # Examine every file in the current input directory.

    for file in os.listdir(folder):

        # Only process files ending in .fa or .ffn.
        # Other file types are skipped.

        if not (file.endswith(".fa") or file.endswith(".ffn")):
            continue

        # Construct the complete path to the FASTA file.

        file_path = os.path.join(folder, file)

        # Remove the file extension to obtain the filename.
        #
        # Example:
        #   genome_001.ffn becomes genome_001

        file_name_only = os.path.splitext(file)[0]

        # Assign the filename as the genome identifier.
        #
        # All genes within this FASTA file will share
        # the same genome_id.

        genome_id = file_name_only


        # ----------------------------------------------------
        # Read FASTA records and calculate gene-level RSCU.
        # ----------------------------------------------------

        # Open the output CSV in append mode so that results
        # from each FASTA file are added to the same CSV.

        with open(gene_rscu_file, "a") as out:

            # SeqIO.parse reads each FASTA record individually.
            #
            # Each record represents one coding sequence (gene).

            for record in SeqIO.parse(file_path, "fasta"):

                # Extract the FASTA record identifier.
                #
                # This is the first identifier in the FASTA
                # header, excluding any description after
                # whitespace.

                gene_id = record.id.split()[0]

                # Convert the nucleotide sequence to uppercase
                # to ensure consistent codon matching.

                seq = str(record.seq).upper()

                # Calculate RSCU values for this individual gene.

                rscu_vals = calc_rscu(seq)

                # Write one row to the output CSV.
                #
                # Each row contains:
                #   1. Original FASTA filename
                #   2. Genome identifier
                #   3. Gene identifier
                #   4. RSCU values for all 61 sense codons

                out.write(
                    f"{file_name_only},{genome_id},{gene_id}," +
                    ",".join(map(str, rscu_vals)) + "\n"
                )


    # Print a message after the current dataset is processed.

    print(
        f"Finished folder {folder}, "
        f"CSV saved in {output_folder}"
    )


# ============================================================
# END OF SCRIPT
# ============================================================

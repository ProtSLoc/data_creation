import csv
import pandas as pd
from Bio.Align import substitution_matrices
from Bio import pairwise2, SeqIO
from tqdm import tqdm
import concurrent.futures
import cupy as cp

# Load meta data
from meta import CSV_DATA_FILE_DIR, CSV_HPA_TESTSET_FILE_DIR

# Load the BLOSUM62 substitution matrix
substitution_matrix = substitution_matrices.load('BLOSUM62')  # use matrix names to load

def calculate_sequence_identity(seq1, seq2):
    """
    Calculate the sequence identity between two protein sequences.

    Parameters:
    - seq1 (str): The first protein sequence.
    - seq2 (str): The second protein sequence.

    Returns:
    - identity_percent (float): The percentage of sequence identity between the two sequences.
    """
    alignments = pairwise2.align.globalds(seq1, seq2, substitution_matrix, -0.5, -0.3)
    best_alignment = alignments[0]  # Assuming you want the best alignment
    aligned_seq1, aligned_seq2, _, _, _ = best_alignment

    # Calculate sequence identity
    seq_len = len(aligned_seq1)
    identical_positions = cp.sum(aligned_seq1 == aligned_seq2)
    identity_percent = (100 * identical_positions / ((len(aligned_seq1) + len(aligned_seq2)) / 2))

    return identity_percent

def load_sequences_from_csv(file_path, sequence_column):
    """
    Load protein sequences from a CSV file.

    Parameters:
    - file_path (str): The path to the CSV file.
    - sequence_column (str): The name of the column containing the protein sequences.

    Returns:
    - sequences (list): A list of protein sequences.
    """
    df = pd.read_csv(file_path)
    sequences = df[sequence_column].tolist()
    return sequences

def fasta_to_csv(input_file, output_file):
    """
    Convert a FASTA file to a CSV file.

    Parameters:
    - input_file (str): The path to the input FASTA file.
    - output_file (str): The path to the output CSV file.
    """
    with open(input_file, "r") as f:
        # Parse the FASTA file
        records = SeqIO.parse(f, "fasta")

        # Write to CSV
        with open(output_file, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)

            # Write header
            writer.writerow(["ID", "Sequence"])

            # Write sequences
            for record in records:
                writer.writerow([record.id, str(record.seq)])

def check_sequence_validity(seq1, dataset2_sequences, threshold):
    """
    Check if a protein sequence is valid by comparing it with a set of reference sequences.

    Parameters:
    - seq1 (str): The protein sequence to check.
    - dataset2_sequences (list): A list of reference protein sequences.
    - threshold (float): The minimum sequence identity threshold.

    Returns:
    - valid (bool): True if the sequence is valid, False otherwise.
    """
    for seq2 in dataset2_sequences:
        score = calculate_sequence_identity(seq1, seq2)
        if score > threshold:
            return False
    return True

def process_sequence(index, row, dataset2_sequences, threshold, results):
    valid = check_sequence_validity(row, dataset2_sequences, threshold)
    results[index] = valid

def main():
    # Load protein sequences from CSV files
    dataset1_sequences = load_sequences_from_csv(CSV_DATA_FILE_DIR, "Sequence")
    dataset2_sequences = load_sequences_from_csv(CSV_HPA_TESTSET_FILE_DIR, "fasta")

    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Submit tasks to the executor
        futures = {executor.submit(process_sequence, index, row, dataset2_sequences, 30, {}): index for index, row in enumerate(tqdm(dataset1_sequences))}
        # Wait for all tasks to complete
        concurrent.futures.wait(futures)

    # Retrieve results in order of completion
    results = {index: future.result() for index, future in futures.items()}
    print(results)

if __name__ == "__main__":
    main()

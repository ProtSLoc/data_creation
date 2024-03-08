import re
from Bio import  SeqIO
import pandas as pd
from meta import FASTA_FILE_DIR, CSV_DATA_FILE_DIR

# Define the keywords to search for
keywords = ["Membrane", "Cytoplasm", "Nucleus", "Extracellular", "Cell membrane", "Mitochondrion",
            "Plastid", "Endoplasmic reticulum", "Lysosome/Vacuole", "Golgi apparatus", "Peroxisome"]

def find_keywords(string, keywords):
    """
    Find the presence of keywords in a given string.

    Args:
        string (str): The input string to search for keywords.
        keywords (list): A list of keywords to search for.

    Returns:
        list: A list of 1s and 0s indicating the presence or absence of each keyword in the string.
    """
    present_keywords = []
    for keyword in keywords:
        if re.search(r'\b' + re.escape(keyword) + r'\b', string, re.IGNORECASE):
            present_keywords.append(1)
        else:
            present_keywords.append(0)
    return present_keywords

def fasta_to_dataframe(fasta_file):
    """
    Convert a FASTA file to a pandas DataFrame.

    Args:
        fasta_file (str): The path to the FASTA file.

    Returns:
        pandas.DataFrame: A DataFrame containing the ID, Sequence, and presence of keywords for each record in the FASTA file.
    """
    records = SeqIO.parse(fasta_file, "fasta")
    data = []
    for record in records:
        # Find which keywords are present in the sample string
        present_keywords = find_keywords(str(record.description), keywords)
        a=[record.id, str(record.seq)]
        a.extend(present_keywords)
        data.append(a)
    coloums=['ID', 'Sequence']
    coloums.extend(keywords)
    print(coloums)
    df = pd.DataFrame(data, columns=coloums)
    return df




df = fasta_to_dataframe(FASTA_FILE_DIR)
df.to_csv(CSV_DATA_FILE_DIR, index=False)
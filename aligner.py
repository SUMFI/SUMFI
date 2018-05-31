import os
from Bio import SeqIO


def needleman_wunsch(x, iterative_method, output_dir, proteins_dir):
    _, protein1, protein2 = x
    protein1_name, protein1_seq = read_sequence(protein1, proteins_dir)
    protein2_name, protein2_seq = read_sequence(protein2, proteins_dir)
    output = os.path.join(output_dir, '{}_{}_{}'.format(protein1_name, protein2_name, type(int)(iterative_method)))

    if iterative_method:
        _needleman_wunsch_iterative(protein1_seq, protein2_seq, output)
    else:
        _needleman_wunsch(protein1_seq, protein2_seq, output)


def smith_waterman(x, iterative_method, output_dir, proteins_dir):
    _, protein1, protein2 = x
    protein1_name, protein1_seq = read_sequence(protein1, proteins_dir)
    protein2_name, protein2_seq = read_sequence(protein2, proteins_dir)
    output = os.path.join(output_dir, '{}_{}_{}'.format(protein1_name, protein2_name, type(int)(iterative_method)))

    if iterative_method:
        _smith_waterman_iterative(protein1_seq, protein2_seq, output)
    else:
        _smith_waterman(protein1_seq, protein2_seq, output)


def read_sequence(protein_name, proteins_dir):
    filepath = os.path.join(proteins_dir, protein_name)
    with open(filepath, 'r') as f:
        seq = ''
        header = f.readline().strip()
        if not header.startswith('>'):
            raise ValueError('Not a FASTA file!')
        for line in f:
            if line[0] == '>':
                break
            seq += line.strip()

        return header[1:5], seq


def _needleman_wunsch(protein1_seq, protein2_seq, output):
    pass


def _smith_waterman(protein1_seq, protein2_seq, output):
    pass


def _needleman_wunsch_iterative(protein1_seq, protein2_seq, output):
    pass


def _smith_waterman_iterative(protein1_seq, protein2_seq, output):
    pass

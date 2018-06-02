import os

import config as cfg


def match_score(alpha, beta):
    if alpha == beta:
        return cfg.match_award
    elif alpha == '-' or beta == '-':
        return cfg.gap_penalty
    else:
        return cfg.mismatch_penalty


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

        return header, seq


def save_alignment(align1, align2, output):
    align1 = align1[::-1]  # reverse sequence 1
    align2 = align2[::-1]  # reverse sequence 2

    i, j = 0, 0

    # calcuate identity, score and aligned sequeces
    symbol = ''
    found = 0
    score = 0
    identity = 0
    for i in range(0, len(align1)):
        # if two AAs are the same, then output the letter
        if align1[i] == align2[i]:
            symbol = symbol + align1[i]
            identity = identity + 1
            score += match_score(align1[i], align2[i])

        # if they are not identical and none of them is gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-':
            score += match_score(align1[i], align2[i])
            symbol += ' '
            found = 0

        # if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':
            symbol += ' '
            score += cfg.gap_penalty

    identity = float(identity) / len(align1) * 100

    with open(output, 'w') as f:
        f.write('Identity = {}\n'.format(identity))
        f.write('Score = {}\n'.format(score))
        f.write('{}\n'.format(align1))
        f.write('{}\n'.format(symbol))
        f.write('{}\n'.format(align2))
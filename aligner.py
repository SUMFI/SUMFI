def needleman_wunsch(x, iterative_method, output_dir, proteins_dir):
    _, protein1, protein2 = x

    if iterative_method:
        _needleman_wunsch_iterative(protein1, protein2, iterative_method, output_dir, proteins_dir)
    else:
        _needleman_wunsch(protein1, protein2, iterative_method, output_dir, proteins_dir)


def smith_waterman(x, iterative_method, output_dir, proteins_dir):
    _, protein1, protein2 = x

    if iterative_method:
        _smith_waterman_iterative(protein1, protein2, iterative_method, output_dir, proteins_dir)
    else:
        _smith_waterman(protein1, protein2, iterative_method, output_dir, proteins_dir)


def _needleman_wunsch(protein1, protein2, iterative_method, output_dir, proteins_dir):
    pass


def _smith_waterman(protein1, protein2, iterative_method, output_dir, proteins_dir):
    pass


def _needleman_wunsch_iterative(protein1, protein2, iterative_method, output_dir, proteins_dir):
    pass


def _smith_waterman_iterative(protein1, protein2, iterative_method, output_dir, proteins_dir):
    pass

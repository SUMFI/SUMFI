import itertools
import os
import re

import click
import pandas as pd


@click.group()
def action():
    pass


@action.command()
@click.option('-p', '--proteins_dir', help='Path to folder with protein files', required=True)
@click.option('-o', '--output_dir', help='Path to output folder', required=True)
def prepare(proteins_dir, output_dir):
    _prepare(proteins_dir, output_dir)


def _prepare(proteins_dir, output_dir):
    list_files = os.listdir(path=proteins_dir)
    list_fasta = [file for file in list_files if re.match('.*\.fasta', file)]
    fasta_combinations = list(itertools.combinations(list_fasta, 2))
    protein_combinations_path = os.path.join(output_dir, 'protein_combinations.csv')
    pd.DataFrame(fasta_combinations).to_csv(protein_combinations_path)


@action.command()
def run():
    _run()


def _run():
    pass


@action.command()
@click.option('-p', '--proteins_dir', help='Path to folder with protein files', required=True)
@click.option('-o', '--output_dir', help='Path to output folder', required=True)
def prepare_run(proteins_dir_path):
    _prepare(proteins_dir_path)
    _run()


if __name__ == "__main__":
    action()

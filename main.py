import click


@click.group()
def action():
    pass


@action.command()
@click.option('-p', '--proteins_dir_path', help='Path to directory with fasta files', required=True)
def prepare(proteins_dir_path):
    _prepare(proteins_dir_path)


def _prepare(proteins_dir_path):
    pass


@action.command()
def run():
    _run()


def _run():
    pass


@action.command()
@click.option('-p', '--proteins_dir_path', help='Path to directory with fasta files', required=True)
def prepare_run(proteins_dir_path):
    _prepare(proteins_dir_path)
    _run()


if __name__ == "__main__":
    action()

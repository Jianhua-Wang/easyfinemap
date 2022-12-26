"""Console script for easy_finemap."""

from easyfinemap.ldref import LDRef

import typer

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
app = typer.Typer(context_settings=CONTEXT_SETTINGS, add_completion=False)

# TODO: 1. validate LD reference panel
# TODO: 2. validate GWAS summary statistics
# TODO: 3. identify lead SNPs
# TODO: 4. fine-mapping


@app.command()
def validate_ldref(
    ldref_path: str = typer.Argument(..., help="The path to the LD reference file."),
    outprefix: str = typer.Argument(..., help="The output prefix."),
    file_type: str = typer.Option("plink", "--file-type", "-f", help="The file type of the LD reference file."),
    mac: int = typer.Option(10, "--mac", "-m", help="The minor allele count threshold."),
    threads: int = typer.Option(1, "--threads", "-t", help="The number of threads."),
) -> None:
    """
    Validate the LD reference file.
    """
    ld = LDRef()
    ld.valid(ldref_path, outprefix, file_type, mac, threads)


@app.command()
def main():
    """Main entrypoint."""
    typer.echo("easy_finemap")
    typer.echo("=" * len("easy_finemap"))
    typer.echo("user-friendly pipeline for GWAS fine-mapping")


if __name__ == "__main__":
    app()

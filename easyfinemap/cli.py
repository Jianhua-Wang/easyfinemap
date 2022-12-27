"""Console script for easy_finemap."""

from pathlib import Path
import typer

import pandas as pd
from easyfinemap.ldref import LDRef
from easyfinemap.sumstat import SumStat

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
app = typer.Typer(context_settings=CONTEXT_SETTINGS, add_completion=False)

# 1. validate LD reference panel
# 2. validate GWAS summary statistics
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
    """Validate the LD reference file."""
    ld = LDRef()
    ld.valid(ldref_path, outprefix, file_type, mac, threads)


@app.command()
def validate_sumstats(
    sumstats_path: Path = typer.Argument(..., help="The path to the GWAS summary statistics file."),
    output: Path = typer.Argument(..., help="The output prefix."),
) -> None:
    """Validate the GWAS summary statistics file."""
    if sumstats_path.exists():
        sumstats = pd.read_csv(sumstats_path, sep="\t")
        typer.echo(f"Loaded {sumstats_path} successfully.")
        typer.echo(f"Number of SNPs: {sumstats.shape[0]}")
        typer.echo(f"Number of columns: {sumstats.shape[1]}")
        typer.echo(f"Columns: {list(sumstats.columns)}")
        valid_sumstats = SumStat(sumstats)
        valid_sumstats = valid_sumstats.standarize()
        if output.suffix == ".gz":
            valid_sumstats.to_csv(output, sep="\t", index=False, compression="gzip")
        else:
            valid_sumstats.to_csv(output, sep="\t", index=False)
        typer.echo(f"Saved the validated summary statistics to {output}.")
    else:
        typer.echo(f"Could not find {sumstats_path}.")
        raise FileNotFoundError(f"Could not find {sumstats_path}.")


if __name__ == "__main__":
    app()

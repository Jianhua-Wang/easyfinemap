"""Console script for easy_finemap."""

import logging
import sys
from enum import Enum
from pathlib import Path
from typing import List, Optional

import pandas as pd
import typer
from rich.console import Console

from easyfinemap import __version__
from easyfinemap.easyfinemap import EasyFinemap
from easyfinemap.ldref import LDRef
from easyfinemap.loci import Loci

# from easyfinemap.sumstat import SumStat

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
app = typer.Typer(context_settings=CONTEXT_SETTINGS, add_completion=False)

# 1. validate LD reference panel
# 2. validate GWAS summary statistics
# 3. identify lead SNPs
# 4. fine-mapping


class LociMethod(str, Enum):
    """The method to identify the lead SNPs."""

    distance = "distance"
    clumping = "clumping"
    conditional = "conditional"


class FinemapMethod(str, Enum):
    """The method to perform fine-mapping."""

    abf = "abf"
    finemap = "finemap"
    paintor = "paintor"
    caviarbf = "caviarbf"
    all = "all"


@app.callback(invoke_without_command=True, no_args_is_help=True)
def main(
    version: bool = typer.Option(False, '--version', '-V', help='Show version.'),
    verbose: bool = typer.Option(False, '--verbose', '-v', help='Show verbose info.'),
):
    """EasyFinemap: A user-friendly tool for fine-mapping."""
    console = Console()
    console.rule("[bold blue]EasyFinemap[/bold blue]")
    console.print(f"Version: {__version__}", justify="center")
    console.print("Author: Jianhua Wang", justify="center")
    console.print("Email: jianhua.mert@gmail.com", justify="center")
    if version:
        typer.echo(f'EasyFinemap version: {__version__}')
        raise typer.Exit()
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logging.info('Verbose mode is on.')
    else:
        logging.getLogger().setLevel(logging.INFO)


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


# @app.command()
# def validate_sumstats(
#     sumstats_path: Path = typer.Argument(..., help="The path to the GWAS summary statistics file."),
#     output: Path = typer.Argument(..., help="The output prefix."),
# ) -> None:
#     """Validate the GWAS summary statistics file."""
#     if sumstats_path.exists():
#         sumstats = pd.read_csv(sumstats_path, sep="\t")
#         typer.echo(f"Loaded {sumstats_path} successfully.")
#         typer.echo(f"Number of SNPs: {sumstats.shape[0]}")
#         typer.echo(f"Number of columns: {sumstats.shape[1]}")
#         typer.echo(f"Columns: {list(sumstats.columns)}")
#         valid_sumstats = SumStat(sumstats)
#         valid_sumstats = valid_sumstats.standarize()
#         if output.suffix == ".gz":
#             valid_sumstats.to_csv(output, sep="\t", index=False, compression="gzip")
#         else:
#             valid_sumstats.to_csv(output, sep="\t", index=False)
#         typer.echo(f"Saved the validated summary statistics to {output}.")
#     else:
#         logging.error(f"No such file of {sumstats_path}.")
#         sys.exit(1)


@app.command()
def get_loci(
    sumstats_path: Path = typer.Argument(..., help="The path to the GWAS summary statistics file."),
    output: str = typer.Argument(..., help="The output prefix."),
    sig_threshold: float = typer.Option(5e-8, "--sig-threshold", "-s", help="The significance threshold."),
    loci_extension: int = typer.Option(500, "--loci-extension", "-l", help="The extension from lead SNPs, in kb"),
    ldblock: str = typer.Option(None, "--ldblock", help="The path to the LD block file."),
    if_merge: bool = typer.Option(
        False, "--merge-loci", "-i", help="Whether to merge the loci, not recommanded for conditional mode."
    ),
    method: LociMethod = typer.Option(
        LociMethod.distance, "--method", "-m", help="The method to identify the lead SNPs."
    ),
    distance: int = typer.Option(50, "--distance", "-d", help="The distance threshold for distance method, in kb."),
    ldref: str = typer.Option(None, "--ldref", help="The path to the LD reference file."),
    clump_kb: int = typer.Option(500, "--clump-kb", "-c", help="The clumping window size, in kb."),
    clump_r2: float = typer.Option(0.1, "--clump-r2", "-r", help="The clumping r2 threshold."),
    sample_size: int = typer.Option(None, "--sample-size", "-n", help="The sample size for conditional method."),
    cojo_window_kb: int = typer.Option(10000, "--cojo-window-kb", help="The cojo window size, in kb."),
    cojo_collinear: float = typer.Option(0.9, "--cojo-collinear", help="The cojo collinear threshold."),
    diff_freq: float = typer.Option(0.2, "--diff-freq", help="The difference in allele frequency threshold."),
    use_ref_eaf: bool = typer.Option(False, "--use-ref-eaf", help="Whether to use the reference panel EAF."),
    only_use_sig_snps: bool = typer.Option(
        False, "--only-use-sig-snps", help="Whether to only use the significant SNPs."
    ),
    threads: int = typer.Option(1, "--threads", "-t", help="The number of threads."),
) -> None:
    """Get the loci from the GWAS summary statistics file."""
    if sumstats_path.exists():
        logging.info(f"Loading {sumstats_path}...")
        sumstats = pd.read_csv(sumstats_path, sep="\t")
        Loci().identify_indep_loci(
            sumstats=sumstats,
            sig_threshold=sig_threshold,
            loci_extend=loci_extension,
            ldblock=ldblock,
            if_merge=if_merge,
            outprefix=output,
            ldref=ldref,
            method=method,
            distance=distance,
            clump_kb=clump_kb,
            clump_r2=clump_r2,
            sample_size=sample_size,
            cojo_window_kb=cojo_window_kb,
            cojo_collinear=cojo_collinear,
            diff_freq=diff_freq,
            use_ref_EAF=use_ref_eaf,
            only_use_sig_snps=only_use_sig_snps,
            threads=threads,
        )
    else:
        logging.error(f"No such file of {sumstats_path}.")
        sys.exit(1)


@app.command()
def fine_mapping(
    sumstats_path: Path = typer.Argument(..., help="The path to the GWAS summary statistics file."),
    loci_path: Path = typer.Argument(..., help="The path to the loci file, generated by get-loci command."),
    lead_snps_path: Path = typer.Argument(..., help="The path to the lead SNPs file, generated by get-loci command."),
    methods: List[FinemapMethod] = typer.Option(..., "--methods", "-m", help="The methods to use."),
    outfile: str = typer.Argument(..., help="The output file."),
    var_prior: float = typer.Option(0.5, "--var-prior", help="The prior variance for the aBF method."),
    conditional: bool = typer.Option(False, "--conditional", "-c", help="Whether to use conditional mode."),
    sample_size: Optional[int] = typer.Option(
        None, "--sample-size", "-n", help="The sample size for conditional mode."
    ),
    ldref: str = typer.Option(None, "--ldref", help="The path to the LD reference file."),
    cond_snps_wind_kb: int = typer.Option(
        10000, "--cond-snps-wind-kb", help="The conditional SNPs window size, in kb."
    ),
    max_causal: int = typer.Option(1, "--max-causal", help="The maximum number of causal SNPs."),
    credible_threshold: Optional[float] = typer.Option(None, "--credible-threshold", help="The credible threshold."),
    credible_method: Optional[str] = typer.Option(
        None, "--credible-method", help="The fine-mapping method for credible set."
    ),
    use_ref_EAF: bool = typer.Option(False, "--use-ref-eaf", help="Whether to use the reference panel EAF."),
    threads: int = typer.Option(1, "--threads", "-t", help="The number of threads."),
) -> None:
    """Fine mapping."""
    if sumstats_path.exists() and loci_path.exists() and lead_snps_path.exists():
        sumstats = pd.read_csv(sumstats_path, sep="\t")
        loci = pd.read_csv(loci_path, sep="\t")
        lead_snps = pd.read_csv(lead_snps_path, sep="\t")
        EasyFinemap().finemap_all_loci(
            sumstats=sumstats,
            loci=loci,
            lead_snps=lead_snps,
            methods=methods,  # type: ignore
            outfile=outfile,
            var_prior=var_prior,
            conditional=conditional,
            sample_size=sample_size,
            ldref=ldref,
            cond_snps_wind_kb=cond_snps_wind_kb,
            max_causal=max_causal,
            credible_threshold=credible_threshold,
            credible_method=credible_method,
            use_ref_EAF=use_ref_EAF,
            threads=threads,
        )
    else:
        logging.error(f"No such file of {sumstats_path} or {loci_path} or {lead_snps_path}.")
        sys.exit(1)


if __name__ == "__main__":
    app()

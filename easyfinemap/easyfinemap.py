"""Perform Fine-mapping.

Perform fine-mapping for a locus using the following methods:
1. LD-free
    1.1. aBF
2. LD-based
    2.1. without annotation
        2.1.1. FINEMAP
        2.1.2. PAINTOR
        2.1.3. CAVIARBF
        TODO: 2.1.4. SuSiE
    2.2. with annotation
        TODO: 2.2.1. PolyFun+SuSiE
3. Support multiple causal variant
4. Support conditional mode
"""

import logging
import os
from pathlib import Path
from subprocess import PIPE, run
from typing import List, Optional

import numpy as np
import pandas as pd
from pathos.pools import _ProcessPool as Pool
from rich.progress import BarColumn, MofNCompleteColumn, Progress, TextColumn, TimeElapsedColumn

from easyfinemap.constant import ColName
from easyfinemap.ldref import LDRef

# from easyfinemap.sumstat import SumStat
from easyfinemap.tools import Tools
from easyfinemap.utils import io_in_tempdir, make_SNPID_unique


class EasyFinemap(object):
    """Main class."""

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger('EasyFinemap')
        tool = Tools()
        self.finemap = tool.finemap
        self.paintor = tool.paintor
        self.gcta = tool.gcta
        self.plink = tool.plink
        self.bcftools = tool.bcftools
        self.caviarbf = tool.caviarbf
        self.model_search = tool.model_search
        self.tmp_root = Path.cwd() / "tmp" / "easyfinemap"
        if not self.tmp_root.exists():
            self.tmp_root.mkdir(parents=True)

    def run_abf(self, sumstats: pd.DataFrame, var_prior: float = 0.2, max_causal: int = 1, **kwargs) -> pd.Series:
        """
        Run ABF.

        calculate the approximate Bayes factor (ABF) from BETA and SE, using the
        formula:
        SNP_BF = sqrt(SE/(SE + W^2))EXP(W^2/(SE + W^2)*(BETA^2/SE^2)/2)
        where W is variance prior, usually set to 0.15 for quantitative traits
        and 0.2 for binary traits.
        the posterior probability of each variant being causal is calculated
        using the formula:
        PP(causal) = SNP_BF / sum(all_SNP_BFs)

        Reference: Asimit, J. L. et al. Eur J Hum Genet (2016)

        Parameters
        ----------
        sumstats : pd.DataFrame
            Summary statistics.
        var_prior : float, optional
            Variance prior, by default 0.2, usually set to 0.15 for quantitative traits
            and 0.2 for binary traits.
        max_causal : int, optional
            Maximum number of causal variants, by default 1

        Returns
        -------
        pd.Series
            The result of ABF.
        """
        if max_causal > 1:
            raise NotImplementedError("ABF only support single causal variant.")
        df = sumstats.copy()
        df["W2"] = var_prior**2
        df["SNP_BF"] = np.sqrt((df[ColName.SE] ** 2 / (df[ColName.SE] ** 2 + df["W2"]))) * np.exp(
            df["W2"] / (df[ColName.BETA] ** 2 + df["W2"]) * (df[ColName.BETA] ** 2 / df[ColName.SE] ** 2) / 2
        )
        df[ColName.PP_ABF] = df["SNP_BF"] / df["SNP_BF"].sum()
        return pd.Series(data=df[ColName.PP_ABF].values, index=df[ColName.SNPID].tolist())

    @io_in_tempdir('./tmp/easyfinemap')
    def run_finemap(
        self,
        sumstats: pd.DataFrame,
        ld_matrix: str,
        sample_size: int,
        max_causal: int = 1,
        temp_dir: Optional[str] = None,
        **kwargs,
    ) -> pd.Series:
        """
        Run FINEMAP.

        Parameters
        ----------
        sumstats : pd.DataFrame
            Summary statistics.
        ld_matrix : str
            Path to LD matrix.
        sample_size : int
            Sample size.
        max_causal : int, optional
            Maximum number of causal variants, by default 1
        temp_dir : Optional[str], optional
            Path to tempdir, by default None

        Returns
        -------
        pd.Series
            The result of FINEMAP.
        """
        if ColName.MAF not in sumstats.columns:
            raise ValueError(f"{ColName.MAF} is required for FINEMAP.")
        finemap_input = sumstats.copy()
        finemap_input[ColName.MAF] = finemap_input[ColName.MAF].replace(0, 0.00001)
        finemap_input = finemap_input[
            [ColName.SNPID, ColName.CHR, ColName.BP, ColName.EA, ColName.NEA, ColName.MAF, ColName.BETA, ColName.SE]
        ]
        finemap_input.rename(
            columns={
                ColName.SNPID: "rsid",
                ColName.CHR: "chromosome",
                ColName.BP: "position",
                ColName.MAF: "maf",
                ColName.BETA: "beta",
                ColName.SE: "se",
                ColName.EA: "allele1",
                ColName.NEA: "allele2",
            },
            inplace=True,
        )
        finemap_input.to_csv(f"{temp_dir}/finemap.z", sep=" ", index=False)
        with open(f"{temp_dir}/finemap.master", "w") as f:
            master_content = [
                f"{temp_dir}/finemap.z",
                ld_matrix,
                f"{temp_dir}/finemap.snp",
                f"{temp_dir}/finemap.config",
                f"{temp_dir}/finemap.cred",
                f"{temp_dir}/finemap.log",
                str(sample_size),
            ]
            f.write("z;ld;snp;config;cred;log;n_samples\n")
            f.write(";".join(master_content))
        cmd = [
            self.finemap,
            "--sss",
            "--in-files",
            f"{temp_dir}/finemap.master",
            "--n-causal-snps",
            str(max_causal),
        ]
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        self.logger.debug(f"run FINEMAP: {' '.join(cmd)}")
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)
        else:
            # if max_causal == 1:
            finemap_res = pd.read_csv(f"{temp_dir}/finemap.snp", sep=" ", usecols=["rsid", "prob"])
            finemap_res = pd.Series(finemap_res["prob"].values, index=finemap_res["rsid"].values)  # type: ignore
            # else:
            #     raise NotImplementedError
            return finemap_res

    @io_in_tempdir('./tmp/easyfinemap')
    def run_paintor(
        self, sumstats: pd.DataFrame, ld_matrix: str, max_causal: int = 1, temp_dir: Optional[str] = None, **kwargs
    ):
        """
        Run PAINTOR.

        Parameters
        ----------
        sumstats : pd.DataFrame
            Summary statistics.
        ld_matrix : str
            Path to LD matrix.
        max_causal : int, optional
            Maximum number of causal variants, by default 1
        temp_dir : Optional[str], optional
            Path to tempdir, by default None

        Returns
        -------
        pd.Series
            The result of PAINTOR.
        """
        paintor_input = sumstats.copy()
        paintor_input["coding"] = 1  # TODO: support paintor annotation mode
        paintor_input["Zscore"] = paintor_input[ColName.BETA] / paintor_input[ColName.SE]
        input_prefix = "paintor.processed"
        paintor_input[[ColName.SNPID, ColName.CHR, ColName.BP, "Zscore"]].to_csv(
            f"{temp_dir}/{input_prefix}", sep=" ", index=False
        )
        paintor_input["coding"].to_csv(f"{temp_dir}/{input_prefix}.annotations", sep=" ", index=False, header=True)
        with open(f"{temp_dir}/{input_prefix}.input", "w") as f:
            f.write(input_prefix)
        ld_matrix_abs_path = os.path.abspath(ld_matrix)
        run(
            ['ln', "-s", ld_matrix_abs_path, f'{temp_dir}/{input_prefix}.ld'],
            stdout=PIPE,
            stderr=PIPE,
            universal_newlines=True,
        )
        cmd = [
            self.paintor,
            "-input",
            f"{temp_dir}/{input_prefix}.input",
            "-out",
            temp_dir,
            "-Zhead",
            "Zscore",
            "-LDname",
            "ld",
            "-enumerate",
            str(max_causal),
            "-in",
            temp_dir,
        ]
        self.logger.debug(f"run PAINTOR: {' '.join(cmd)}")
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)
        else:
            paintor_res = pd.read_csv(
                f"{temp_dir}/paintor.processed.results", sep=" ", usecols=["SNPID", "Posterior_Prob"]
            )
            paintor_res = pd.Series(paintor_res["Posterior_Prob"].values, index=paintor_res["SNPID"].tolist())
            return paintor_res

    @io_in_tempdir('./tmp/easyfinemap')
    def run_caviarbf(
        self, sumstats: pd.DataFrame, ld_matrix: str, max_causal: int = 1, temp_dir: Optional[str] = None, **kwargs
    ):
        """
        Run CAVIAR-BF.

        Parameters
        ----------
        sumstats : pd.DataFrame
            Summary statistics.
        ld_matrix : str
            Path to LD matrix.
        max_causal : int, optional
            Maximum number of causal variants, by default 1
        temp_dir : Optional[str], optional
            Path to tempdir, by default None

        Returns
        -------
        pd.Series
            The result of CAVIAR-BF.
        """
        caviar_input = sumstats.copy()
        caviar_input[ColName.Z] = caviar_input[ColName.BETA] / caviar_input[ColName.SE]
        caviar_input[[ColName.SNPID, ColName.Z]].to_csv(f"{temp_dir}/caviar.input", sep=" ", index=False, header=False)
        n_variants = caviar_input.shape[0]
        cmd = [
            self.caviarbf,
            "-z",
            f"{temp_dir}/caviar.input",
            "-r",
            ld_matrix,
            "-t",
            "0",
            "-a",
            "0.1281429",
            "-n",
            str(n_variants),
            "-c",
            str(max_causal),
            "-o",
            f"{temp_dir}/caviar.output",
        ]
        self.logger.debug(f"run CAVIAR-BF: {' '.join(cmd)}")
        run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        cmd = [
            self.model_search,
            "-i",
            f"{temp_dir}/caviar.output",
            "-m",
            str(n_variants),
            "-p",
            "0",
            "-o",
            f"{temp_dir}/caviar.prior0",
        ]
        self.logger.debug(f"run CAVIAR-BF: {' '.join(cmd)}")
        res = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        if res.returncode != 0:
            self.logger.error(res.stderr)
            raise RuntimeError(res.stderr)
        else:
            caviar_res = pd.read_csv(f"{temp_dir}/caviar.prior0.marginal", sep=" ", header=None)
            caviar_res.sort_values(by=0, inplace=True)  # type: ignore
            caviar_res = pd.Series(caviar_res[1].values, index=caviar_input[ColName.SNPID].tolist())
            return caviar_res

    def cond_sumstat(
        self,
        sumstats: pd.DataFrame,
        lead_snp: str,
        lead_snps: pd.DataFrame,
        ldref: str,
        sample_size: int,
        use_ref_EAF: bool = False,
        cond_snps_wind_kb: int = 1000,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Conditional sumstat.

        Parameters
        ----------
        sumstats : pd.DataFrame
            Summary statistics.
        lead_snp : str
            Lead SNP.
        ldref : str
            Path to LD reference.
        sample_size : int
            Sample size.
        use_ref_EAF : bool, optional
            Use reference EAF, by default False
        cond_snps_wind_kb : int, optional
            Conditional SNPs window in kb, by default 1000

        Returns
        -------
        pd.DataFrame
            Conditional sumstat.
        """
        if lead_snp is None:
            raise ValueError("Lead SNP is required for conditional finemapping")
        if lead_snps is None:
            raise ValueError("Lead SNPs are required for conditional finemapping")
        lead_snp_chr = lead_snps.loc[lead_snps[ColName.SNPID] == lead_snp, ColName.CHR].values[0]
        lead_snp_bp: int = lead_snps.loc[lead_snps[ColName.SNPID] == lead_snp, ColName.BP].values[0]  # type: ignore
        cond_snps = lead_snps[
            (lead_snps[ColName.CHR] == lead_snp_chr)
            & (lead_snps[ColName.BP] >= lead_snp_bp - cond_snps_wind_kb * 1000)
            & (lead_snps[ColName.BP] <= lead_snp_bp + cond_snps_wind_kb * 1000)
            & (lead_snps[ColName.SNPID] != lead_snp)
        ]
        if cond_snps.empty:
            self.logger.debug(f"No conditional SNPs found for {lead_snp}")
            cond_res = sumstats.copy()
            cond_res[ColName.COJO_BETA] = cond_res[ColName.BETA]
            cond_res[ColName.COJO_SE] = cond_res[ColName.SE]
            cond_res[ColName.COJO_P] = cond_res[ColName.P]
        else:
            ld = LDRef()
            cond_res = ld.cojo_cond(sumstats, cond_snps, ldref, sample_size, use_ref_EAF)  # type: ignore
        return cond_res

    def prepare_ld_matrix(
        self,
        sumstats: pd.DataFrame,
        ldref: str,
        outprefix: str,
        use_ref_EAF: bool = False,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Prepare LD matrix.

        Parameters
        ----------
        sumstats : pd.DataFrame
            Summary statistics.
        ldref : str
            Path to LD reference.
        outprefix : str
            Output prefix.
        use_ref_EAF : bool, optional
            Use reference EAF, by default False

        Returns
        -------
        pd.DataFrame
            LD matrix.
        """
        if ldref is None:
            raise ValueError("LD reference is required for LD-based finemapping")
        ld = LDRef()
        sumstats_ol = ld.intersect(sumstats, ldref, outprefix, use_ref_EAF)
        ld.make_ld(outprefix, outprefix)
        return sumstats_ol

    def get_credset(
        self,
        finemap_res: pd.DataFrame,
        max_causal: int,
        credible_threshold: Optional[float] = None,
        credible_method: Optional[str] = None,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Get credible set.

        Parameters
        ----------
        finemap_res : pd.DataFrame
            Finemap results.
        max_causal : int
            Maximum number of causal variants.
        credible_threshold : Optional[float], optional
            Credible threshold, by default None
        credible_method : Optional[str], optional
            Credible set method, by default None

        Returns
        -------
        pd.DataFrame
            Credible set.
        """
        if credible_threshold is None:
            return finemap_res
        else:
            credible_threshold = credible_threshold * max_causal
            if credible_method:
                pp_col = f"PP_{credible_method.upper()}"
                credible_set = finemap_res.sort_values(pp_col, ascending=False)
                credible_set = finemap_res.sort_values(by=pp_col, ascending=False)
                credible_set = credible_set[credible_set[pp_col].shift().fillna(0).cumsum() <= credible_threshold]
            else:
                raise ValueError("Must specify credible set method when credible threshold is specified")
        return credible_set.reset_index(drop=True)

    @io_in_tempdir('./tmp/easyfinemap')
    def finemap_locus(
        self,
        sumstats: pd.DataFrame,
        methods: List[str],
        lead_snp: str,
        conditional: bool = False,
        temp_dir: Optional[str] = None,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Finemap a locus.

        1. Check if LD is needed, abf, susie, polyfun+susie do not need LD.
        2. If LD is needed, intersect the locus with the LD reference and make the LD matrix.
        3. Run the finemapping method.
        4. Get the finemapping results.
        5. Merge the finemapping results with the input sumstats.
        6. Return the credible set or full summary statistics with posterior probabilities.

        Parameters
        ----------
        sumstats : pd.DataFrame
            Summary statistics.
        methods : List[str]
            Finemapping methods.
        lead_snp : str
            Lead SNP.
        conditional : bool, optional
            Conditional finemapping, by default False
        temp_dir : Optional[str], optional
            Temporary directory, by default None

        Returns
        -------
        pd.DataFrame
            Finemapping results.
        """
        if conditional:
            cond_res = self.cond_sumstat(sumstats=sumstats, lead_snp=lead_snp, **kwargs)
            fm_input = cond_res.copy()
            fm_input[ColName.BETA] = cond_res[ColName.COJO_BETA]
            fm_input[ColName.SE] = cond_res[ColName.COJO_SE]
            fm_input[ColName.P] = cond_res[ColName.COJO_P]
            out_sumstats = sumstats.merge(
                cond_res[[ColName.SNPID, ColName.COJO_BETA, ColName.COJO_SE, ColName.COJO_P]],
                on=ColName.SNPID,
                how="left",
            )
            max_causal = kwargs.get("max_causal", 1)
            if max_causal > 1:
                self.logger.warning("Conditional finemapping does not support multiple causal variants")
        else:
            fm_input = sumstats.copy()
            out_sumstats = sumstats.copy()

        allowed_methods = ["abf", "finemap", "paintor", "caviarbf"]
        if "all" in methods:
            methods = allowed_methods
        fm_input_ol = fm_input.copy()
        for method in methods:
            if method == "abf":
                abf_pp = self.run_abf(sumstats=fm_input, **kwargs)
                out_sumstats[ColName.PP_ABF] = out_sumstats[ColName.SNPID].map(abf_pp)
            elif method in ["finemap", "paintor", "caviarbf"]:
                ld_matrix = f"{temp_dir}/intersc.ld"
                if not os.path.exists(ld_matrix):
                    # TODO: reduce the number of SNPs when using paintor and caviarbf in multiple causal variant mode
                    fm_input_ol = self.prepare_ld_matrix(sumstats=fm_input, outprefix=f"{temp_dir}/intersc", **kwargs)
                if method == "finemap":
                    finemap_pp = self.run_finemap(sumstats=fm_input_ol, ld_matrix=ld_matrix, **kwargs)
                    out_sumstats[ColName.PP_FINEMAP] = out_sumstats[ColName.SNPID].map(finemap_pp)
                elif method == "paintor":
                    paintor_pp = self.run_paintor(sumstats=fm_input_ol, ld_matrix=ld_matrix, **kwargs)
                    out_sumstats[ColName.PP_PAINTOR] = out_sumstats[ColName.SNPID].map(paintor_pp)
                elif method == "caviarbf":
                    caviarbf_pp = self.run_caviarbf(sumstats=fm_input_ol, ld_matrix=ld_matrix, **kwargs)
                    out_sumstats[ColName.PP_CAVIARBF] = out_sumstats[ColName.SNPID].map(caviarbf_pp)
            else:
                raise ValueError(f"Method {method} is not supported")

        credible_set = self.get_credset(finemap_res=out_sumstats, **kwargs)
        credible_set[ColName.LEAD_SNP] = lead_snp
        return credible_set

    def finemap_all_loci(
        self,
        sumstats: pd.DataFrame,
        loci: pd.DataFrame,
        lead_snps: pd.DataFrame,
        methods: List[str],
        var_prior: float = 0.2,
        conditional: bool = False,
        sample_size: Optional[int] = None,
        ldref: Optional[str] = None,
        cond_snps_wind_kb: int = 1000,
        max_causal: int = 1,
        credible_threshold: Optional[float] = None,
        credible_method: Optional[str] = None,
        use_ref_EAF: bool = False,
        outfile: Optional[str] = None,
        threads: int = 1,
    ):
        """
        Perform finemapping for all loci.

        Parameters
        ----------
        sumstats : pd.DataFrame
            Summary statistics.
        loci : pd.DataFrame
            Loci.
        lead_snps : pd.DataFrame
            Lead SNPs.
        methods : List[str]
            Finemapping methods.
        var_prior : float, optional
            Variance prior, by default 0.2
        conditional : bool, optional
            Conditional finemapping, by default False
        sample_size : Optional[int], optional
            Sample size, by default None
        ldref : Optional[str], optional
            LD reference, by default None
        cond_snps_wind_kb : int, optional
            Conditional SNPs window, by default 1000
        max_causal : int, optional
            Maximum number of causal variants, by default 1
        credible_threshold : Optional[float], optional
            Credible threshold, by default None
        credible_method : Optional[str], optional
            Credible method, by default None
        use_ref_EAF : bool, optional
            Use reference EAF, by default False
        outfile : Optional[str], optional
            Output file, by default None
        threads : int, optional
            Number of threads, by default 1
        """
        sumstats = make_SNPID_unique(sumstats)
        if credible_threshold and credible_method is None and methods != ["all"] and len(methods) == 1:
            credible_method = methods[0]
        kwargs_list = []
        for chrom, start, end, lead_snp in loci[[ColName.CHR, ColName.START, ColName.END, ColName.LEAD_SNP]].values:
            locus_sumstats = sumstats.loc[
                (sumstats[ColName.CHR] == chrom) & (sumstats[ColName.BP] >= start) & (sumstats[ColName.BP] <= end)
            ]
            kwargs = {
                "sumstats": locus_sumstats,
                "lead_snp": lead_snp,
                "lead_snps": lead_snps,
                "methods": methods,
                "var_prior": var_prior,
                "conditional": conditional,
                "sample_size": sample_size,
                "ldref": ldref.format(chrom=chrom) if ldref else None,
                "cond_snps_wind_kb": cond_snps_wind_kb,
                "max_causal": max_causal,
                "credible_threshold": credible_threshold,
                "credible_method": credible_method,
                "use_ref_EAF": use_ref_EAF,
            }
            kwargs_list.append(kwargs)
        ef = EasyFinemap()
        output = []
        with Progress(
            TextColumn("{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TimeElapsedColumn(),
            auto_refresh=True,
        ) as progress:
            with Pool(threads) as p:
                task = progress.add_task("Perform Fine-mapping...", total=len(loci))
                results = [p.apply_async(ef.finemap_locus, kwds=kwargs) for kwargs in kwargs_list]
                for res in results:
                    progress.update(task, advance=1)
                    progress.refresh()
                    output.append(res.get())
        output_df = pd.concat(output, ignore_index=True)
        if outfile:
            output_df.to_csv(outfile, sep="\t", index=False)
        else:
            return output_df

"""Perform Fine-mapping.

Perform fine-mapping for a locus using the following methods:
1. LD-free
2. LD-based
    2.1. without annotation
        2.1.1. FINEMAP
        2.1.2. PAINTOR
        2.1.3. CAVIARBF
        2.1.4. SuSiE
    2.2. with annotation
        2.2.1. PolyFun+SuSiE
"""

import logging
from pathlib import Path
from subprocess import PIPE, run
from typing import List, Optional

import pandas as pd
import numpy as np

from easyfinemap.constant import ColName
from easyfinemap.ldref import LDRef
from easyfinemap.loci import Loci
from easyfinemap.sumstat import SumStat
from easyfinemap.tools import Tools
from easyfinemap.utils import get_significant_snps, io_in_tempdir


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
        self.tmp_root = Path.cwd() / "tmp" / "finemapping"
        if not self.tmp_root.exists():
            self.tmp_root.mkdir(parents=True)

    def run_abf(self, sumstats: pd.DataFrame, var_prior: float = 0.2) -> pd.Series:
        """
        Run ABF.

        calculate the approximate Bayes factor (ABF) from BETA and SE, using the
        formula:
        SNP_BF = sqrt(BETA/(BETA + W^2))EXP(W^2/(BETA + W^2)*(BETA^2/SE^2)/2)
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

        Returns
        -------
        pd.Series
            The result of ABF.
        """
        df = sumstats.copy()
        df["W2"] = var_prior**2
        df[ColName.BETA] = df[ColName.BETA].abs()
        df["SNP_BF"] = np.sqrt((df[ColName.BETA] / (df[ColName.BETA] + df["W2"]))) * np.exp(
            df["W2"] / (df[ColName.BETA] + df["W2"]) * (df[ColName.BETA] ** 2 / df[ColName.SE] ** 2) / 2
        )
        df[ColName.PP_ABF] = df["SNP_BF"] / df["SNP_BF"].sum()
        return pd.Series(data=df[ColName.PP_ABF].values, index=df[ColName.SNPID].tolist())

    @io_in_tempdir('./tmp/finemapping')
    def run_finemap(
        self,
        sumstats: pd.DataFrame,
        ld_matrix: str,
        sample_size: int,
        max_causal: int = 1,
        temp_dir: Optional[str] = None,
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
            finemap_res = pd.read_csv(f"{temp_dir}/finemap.config", sep=" ", usecols=["config", "prob"])
            finemap_res = pd.Series(finemap_res["prob"].values, index=finemap_res["config"].values)  # type: ignore
            return finemap_res

    def run_paintor(self, sumstat: pd.DataFrame, ldref: LDRef, out: str):
        """Run PAINTOR."""
        self.logger.info("Running PAINTOR")
        raise NotImplementedError

    def cond_sumstat(
        self,
        sumstats: pd.DataFrame,
        lead_snp: str,
        lead_snps: pd.DataFrame,
        ldref: str,
        sample_size: int,
        use_ref_EAF: bool = False,
        cond_snps_wind_kb: int = 1000,
    ) -> pd.DataFrame:
        """Conditional sumstat."""
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
    ) -> pd.DataFrame:
        """Prepare LD matrix."""
        if ldref is None:
            raise ValueError("LD reference is required for LD-based finemapping")
        ld = LDRef()
        sumstats_ol = ld.intersect(sumstats, ldref, outprefix, use_ref_EAF)
        ld.make_ld(outprefix, outprefix)
        return sumstats_ol

    def get_credset(
        self,
        finemap_res: pd.DataFrame,
        credible_threshold: Optional[float] = None,
        credible_method: Optional[str] = None,
    ) -> pd.DataFrame:
        """Get credible set."""
        if credible_threshold is None:
            return finemap_res
        else:
            if credible_method:
                credible_set = finemap_res.sort_values(credible_method, ascending=False)
                credible_set = credible_set[credible_set[credible_method].cumsum() >= credible_threshold]
            else:
                raise ValueError("Must specify credible set method when credible threshold is specified")
        return credible_set

    @io_in_tempdir('./tmp/finemapping')
    def finemap_locus(
        self,
        sumstats: pd.DataFrame,
        methods: List[str],
        var_prior: float = 0.2,
        conditonal: bool = False,
        sample_size: Optional[int] = None,
        max_causal: int = 1,
        credible_threshold: Optional[float] = None,
        credible_set_method: Optional[str] = None,
        use_ref_EAF: bool = False,
        out: Optional[str] = None,
        temp_dir: Optional[str] = None,
        **kwargs,
    ) -> Optional[pd.DataFrame]:
        """
        Finemap a locus.

        1. Check if LD is needed, abf, susie, polyfun+susie do not need LD.
        2. If LD is needed, intersect the locus with the LD reference and make the LD matrix.
        3. Run the finemapping method.
        4. Get the finemapping results.
        5. Merge the finemapping results with the input sumstats.
        6. Return the credible set or full summary statistics with posterior probabilities.
        """
        if conditonal:
            cond_res = self.cond_sumstat(**kwargs)
            fm_input = cond_res.copy()
            fm_input[ColName.BETA] = cond_res[ColName.COJO_BETA]
            fm_input[ColName.SE] = cond_res[ColName.COJO_SE]
            fm_input[ColName.P] = cond_res[ColName.COJO_P]
            out_sumstats = sumstats.merge(
                cond_res[[ColName.SNPID, ColName.COJO_BETA, ColName.COJO_SE, ColName.COJO_P]],
                on=ColName.SNPID,
                how="left",
            )
        else:
            fm_input = sumstats.copy()
            out_sumstats = sumstats.copy()

        allowed_methods = ["abf", "finemap"]
        if methods == ["all"]:
            methods = allowed_methods
        for method in methods:
            if method == "abf":
                if max_causal > 1:
                    raise ValueError("ABF only supports single causal variant")
                abf_pp = self.run_abf(fm_input, var_prior)
                out_sumstats[ColName.PP_ABF] = out_sumstats[ColName.SNPID].map(abf_pp)
            if method in ["finemap", "paintor"]:
                if sample_size is None:
                    raise ValueError("Sample size is required for conditional finemapping")
                if ColName.MAF not in fm_input.columns and not use_ref_EAF:
                    raise ValueError("MAF is required for conditional finemapping, please set use_ref_EAF to True")
                fm_input_ol = self.prepare_ld_matrix(sumstats=fm_input, outprefix=f"{temp_dir}/intersc", **kwargs)
                if method == "finemap":
                    if max_causal > 1:
                        raise NotImplementedError
                    finemap_pp = self.run_finemap(fm_input_ol, f"{temp_dir}/intersc.ld", sample_size, max_causal)
                    out_sumstats[ColName.PP_FINEMAP] = out_sumstats[ColName.SNPID].map(finemap_pp)
                elif method == "paintor":
                    raise NotImplementedError
            else:
                raise ValueError(f"Method {method} is not supported")

        if credible_threshold is not None:
            if credible_set_method is None:
                raise ValueError("Credible set method is required")
            if credible_set_method == "abf":
                credible_set = out_sumstats[out_sumstats[ColName.PP_ABF] >= credible_threshold]
            elif credible_set_method == "finemap":
                credible_set = out_sumstats[out_sumstats[ColName.PP_FINEMAP] >= credible_threshold]
            else:
                raise ValueError(f"Method {credible_set_method} is not supported")
            credible_set.to_csv(out, sep="\t", index=False)

    def finemap_all_loci(
        self,
        sumstats: pd.DataFrame,
        methods: List[str],
        var_prior: float = 0.2,
        lead_snp: Optional[str] = None,
        lead_snps: Optional[pd.DataFrame] = None,
        conditonal: bool = False,
        sample_size: Optional[int] = None,
        ldref: Optional[str] = None,
        cond_snps_wind_kb: int = 1000,
        max_causal: int = 1,
        credible_threshold: Optional[float] = None,
        credible_set_method: Optional[str] = None,
        use_ref_EAF: bool = False,
        out: Optional[str] = None,
        temp_dir: Optional[str] = None,
    ):
        """Perform finemapping for all loci."""
        raise NotImplementedError

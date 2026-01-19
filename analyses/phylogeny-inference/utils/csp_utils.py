import anndata as ad
import gzip
import numpy as np
import os
import pandas as pd
import pysam
import scipy as sp
from scipy import io
from scipy import sparse

def csp_load_data(data_dir, prefix="", is_genotype=False, is_gzip=True):
    """
    Loads cellSNP output with support for file prefixes and gzip.
    """
    suffix = ".gz" if is_gzip else ""

    # Construct filenames with prefix
    base_vcf_fn = f"{prefix}cellSNP.base.vcf{suffix}"
    cell_vcf_fn = f"{prefix}cellSNP.cells.vcf{suffix}"
    samples_fn = f"{prefix}cellSNP.samples.tsv{suffix}"
    ad_fn = f"{prefix}cellSNP.tag.AD.mtx{suffix}"
    dp_fn = f"{prefix}cellSNP.tag.DP.mtx{suffix}"
    oth_fn = f"{prefix}cellSNP.tag.OTH.mtx{suffix}"

    base_vcf_path = os.path.join(data_dir, base_vcf_fn)

    # Fallback checking
    if not os.path.exists(base_vcf_path) and is_gzip:
        # Try without gz if default fails
        if os.path.exists(os.path.join(data_dir, f"{prefix}cellSNP.base.vcf")):
            base_vcf_path = os.path.join(data_dir, f"{prefix}cellSNP.base.vcf")
            suffix = ""
            is_gzip = False

    base_vcf_comment, base_vcf = csp_load_vcf(base_vcf_path)

    cell_vcf, cell_vcf_comment = None, None
    if is_genotype:
        cell_vcf_comment, cell_vcf = csp_load_vcf(os.path.join(data_dir, cell_vcf_fn))

    samples = csp_load_samples(os.path.join(data_dir, samples_fn))
    AD_mtx = csp_load_matrix(os.path.join(data_dir, ad_fn))
    DP_mtx = csp_load_matrix(os.path.join(data_dir, dp_fn))
    OTH_mtx = csp_load_matrix(os.path.join(data_dir, oth_fn))

    adata = ad.AnnData(
        X = AD_mtx,
        obs = base_vcf,
        var = samples)

    adata.uns["is_gzip"] = is_gzip
    adata.uns["is_genotype"] = is_genotype
    adata.uns["base_vcf_comment"] = base_vcf_comment

    if cell_vcf is None:
        adata.uns["cell_vcf_comment"] = None
    else:
        adata.obsm["cell_vcf"] = cell_vcf
        adata.uns["cell_vcf_comment"] = cell_vcf_comment

    adata.layers["DP"] = DP_mtx
    adata.layers["OTH"] = OTH_mtx

    return adata


def csp_load_matrix(fn):
    mtx = None
    try:
        mtx = sp.io.mmread(fn)
    except:
        mtx = io.mmread(fn)
    mtx = mtx.toarray()    # convert from sparse matrix to ndarray to support slicing.
    return mtx


def csp_load_samples(fn):
    # pandas handles .gz automatically via 'infer' or checking extension
    df = pd.read_csv(fn, header=None, sep="\t")
    df.columns = ["cell"]
    return df


def csp_load_vcf(fn):
    # load comment
    fp = None
    if fn.endswith(".gz") or fn.endswith(".GZ"):
        fp = gzip.open(fn, "rt")
    else:
        fp = open(fn, "r")

    comment = ""
    pre_line = None
    for line in fp:
        if not line or line[0] != "#":
            break
        pre_line = line
        comment += line

    fp.close()

    if not pre_line:
        raise IOError(f"Header parsing failed for {fn}")
    assert len(pre_line) > 6
    assert pre_line[:6] == "#CHROM"

    # load content
    content = pd.read_csv(fn, sep="\t", header=None, comment="#")
    content.columns = pre_line.strip().split("\t")
    content.rename(columns={'#CHROM': 'CHROM'}, inplace=True)

    return (comment, content)


def csp_save_data(adata, out_dir, prefix=""):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    is_gzip = adata.uns["is_gzip"]
    is_genotype = adata.uns["is_genotype"]
    suffix = ".gz" if is_gzip else ""

    # Construct filenames
    base_vcf_fn = os.path.join(out_dir, f"{prefix}cellSNP.base.vcf{suffix}")
    cell_vcf_fn = os.path.join(out_dir, f"{prefix}cellSNP.cells.vcf{suffix}")
    samples_fn = os.path.join(out_dir, f"{prefix}cellSNP.samples.tsv{suffix}")
    ad_fn = os.path.join(out_dir, f"{prefix}cellSNP.tag.AD.mtx{suffix}")
    dp_fn = os.path.join(out_dir, f"{prefix}cellSNP.tag.DP.mtx{suffix}")
    oth_fn = os.path.join(out_dir, f"{prefix}cellSNP.tag.OTH.mtx{suffix}")

    # Prepare OBS for VCF saving
    vcf_df = adata.obs.copy()
    if 'CHROM' in vcf_df.columns:
        vcf_df.rename(columns={'CHROM': '#CHROM'}, inplace=True)

    csp_save_vcf(vcf_df,
        comment = adata.uns["base_vcf_comment"],
        fn = base_vcf_fn,
        is_gzip = is_gzip)

    if is_genotype:
        csp_save_vcf(adata.obsm["cell_vcf"],
            comment = adata.uns["cell_vcf_comment"],
            fn = cell_vcf_fn,
            is_gzip = is_gzip)

    csp_save_samples(adata.var, samples_fn)

    csp_save_matrix(adata.X, ad_fn)
    csp_save_matrix(adata.layers["DP"], dp_fn)
    csp_save_matrix(adata.layers["OTH"], oth_fn)


def csp_save_matrix(mtx, fn):
    mtx = sparse.csr_matrix(mtx)

    # Check if we need to compress
    if fn.endswith(".gz"):
        # Define a temporary raw filename (remove .gz)
        raw_fn = fn[:-3]

        # 1. Write uncompressed data to the raw filename
        io.mmwrite(raw_fn, mtx)

        # 2. Compress the raw file to the final destination (fn)
        with open(raw_fn, 'rb') as f_in:
            with gzip.open(fn, 'wb') as f_out:
                f_out.writelines(f_in)

        # 3. Cleanup the raw file
        if os.path.exists(raw_fn):
            os.remove(raw_fn)
    else:
        # If not gzipped, just write directly
        io.mmwrite(fn, mtx)


def csp_save_samples(df, fn):
    # Handle compression based on filename extension
    if fn.endswith(".gz"):
        df.to_csv(fn, sep="\t", header=False, index=False, compression='gzip')
    else:
        df.to_csv(fn, sep="\t", header=False, index=False)


def csp_save_vcf(df, comment, fn, is_gzip=True):
    # Using a temp file to write the DF content then merging with comment into the final file
    df_file = fn + ".tmp"
    df.to_csv(df_file, sep="\t", header=False, index=False)

    fp = pysam.BGZFile(fn, "w") if is_gzip else open(fn, "w")
    fp.write(comment.encode())

    with open(df_file, "r") as df_fp:
        for line in df_fp:
            fp.write(line.encode())
    fp.close()

    os.remove(df_file)

#!/usr/bin/env python
# coding: utf-8

import sys
import os
import shutil
import pandas as pd
import numpy as np
import fastcluster
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score
from pathlib import Path

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['figure.max_open_warning'] = 100
import matplotlib.pyplot as plt
import seaborn as sns

# Increase recursion depth carefully (clustermap may recurse)
sys.setrecursionlimit(3000)

def detect_matrix_type(df: pd.DataFrame):
    """Detect matrix type: binary (0/1) or ternary (0/1/3)."""
    vals = np.unique(df.values)
    try:
        vals_f = vals.astype(float)
    except Exception:
        raise ValueError("Matrix contains non-numeric values.")

    unique_set = set(np.round(vals_f, 6))
    if unique_set.issubset({0.0, 1.0}): return "binary"
    if unique_set.issubset({0.0, 1.0, 3.0}): return "ternary"
    if np.any((vals_f > 0.0) & (vals_f < 1.0)):
        raise ValueError("One-hot encoded/fractional matrices are NOT supported. Input must be binary (0/1) or ternary (0/1/3).")
    raise ValueError(f"Matrix contains unexpected values: {unique_set}. Expected binary or ternary.")

def weighted_nmf(V, W, rank, max_iter=500, tol=1e-4):
    np.random.seed(42)
    n, m = V.shape
    H = np.random.rand(rank, m)
    U = np.random.rand(n, rank)
    H = H / (np.linalg.norm(H, ord='fro') + 1e-10)
    U = U / (np.linalg.norm(U, ord='fro') + 1e-10)

    for i in range(max_iter):
        num_H = U.T @ (W * V)
        den_H = U.T @ (W * (U @ H)) + 1e-10
        H *= (num_H / den_H)
        num_U = (W * V) @ H.T
        den_U = ((W * (U @ H)) @ H.T) + 1e-10
        U *= (num_U / den_U)
        recon = U @ H
        err = np.linalg.norm(W * (V - recon), ord='fro')
        if err < tol: break
    return U, H

def optimal_k_selection(V, W, k_min=2, k_max=10):
    best_k = k_min
    best_score = -np.inf
    if k_max < k_min: k_max = k_min

    for k in range(k_min, k_max + 1):
        try:
            U, H = weighted_nmf(V, W, k)
            labels = np.argmax(U, axis=1)
            if len(np.unique(labels)) < 2: continue
            silhouette = silhouette_score(U, labels, metric='cosine')
            penalty = np.log(k)

            compactness = 0.0
            unique_labels = np.unique(labels)
            valid_clusters = [i for i in unique_labels if len(U[labels == i]) > 1]
            if valid_clusters:
                compactness = np.mean([
                    np.linalg.norm(U[labels == i] - U[labels == i].mean(axis=0))
                    for i in valid_clusters
                ])

            hybrid_score = silhouette * penalty + 0.1 * compactness
            if hybrid_score > best_score:
                best_score = hybrid_score
                best_k = k
        except Exception:
            continue
    return best_k

def save_heatmap(data, clone_ids, narrow_clone_ids, output_path, row_cluster):
    sorted_idx = np.lexsort((pd.factorize(narrow_clone_ids)[0], clone_ids))
    sorted_data = data.iloc[sorted_idx]
    sorted_clone_ids = np.array(clone_ids)[sorted_idx]
    sorted_narrow_clone_ids = np.array(narrow_clone_ids)[sorted_idx]

    sorted_annotations = pd.DataFrame({
        'Clone ID': sorted_clone_ids,
        'Narrow Clone ID': sorted_narrow_clone_ids
    }, index=sorted_data.index)
    sorted_annotations['Clone ID'] = sorted_annotations['Clone ID'].astype(str)
    sorted_annotations['Narrow Clone ID'] = sorted_annotations['Narrow Clone ID'].astype(str)

    unique_broad = sorted(sorted_annotations['Clone ID'].unique())
    unique_narrow = sorted(sorted_annotations['Narrow Clone ID'].unique())
    palette_broad = sns.color_palette("husl", len(unique_broad))
    palette_narrow = sns.color_palette("husl", len(unique_narrow))

    broad_to_color = dict(zip(unique_broad, palette_broad))
    narrow_to_color = dict(zip(unique_narrow, palette_narrow))

    row_colors_df = pd.DataFrame({
        'Clone ID': sorted_annotations['Clone ID'].map(broad_to_color).fillna('white'),
        'Narrow Clone ID': sorted_annotations['Narrow Clone ID'].map(narrow_to_color).fillna('white')
    }, index=sorted_data.index)

    linkage_matrix = linkage(sorted_data, method='ward') if row_cluster else None
    sns.clustermap(
        sorted_data, row_linkage=linkage_matrix, row_colors=row_colors_df,
        col_cluster=True, row_cluster=row_cluster, figsize=(10, 10),
        cmap="viridis", xticklabels=False, yticklabels=False
    )
    plt.savefig(output_path)
    plt.close()

def run_hierarchical(args, scratch_output_dir, SNV_CALLING_DIR, MODULE_RESULTS_DIR, MODULE_PLOTS_DIR):
    print(f"\n--- Starting Hierarchical SNP Clustering for {args.sample} ---")

    input_file = SNV_CALLING_DIR / args.project / f"{args.sample}.{args.cell_type}.cellSNP.snp.mtx.gz"

    if not input_file.exists():
        raise FileNotFoundError(f"Required input matrix not found: {input_file}")

    print(f"Reading matrix: {input_file}")

    # --- CORRECTION: REMOVED .transpose() ---
    # Input file format is: Rows = Cells, Columns = SNPs.
    # We must keep Cells as Rows (V = Cells x SNPs) for clustering to group Cells.
    df = pd.read_csv(input_file, sep='\t', index_col=0)

    df = df.apply(pd.to_numeric, errors='coerce').fillna(0.0)

    matrix_type = detect_matrix_type(df)
    print(f"Detected matrix type: {matrix_type}")

    num_cells = df.shape[0]
    if num_cells == 0:
        print("No cells available. Exiting.")
        return

    if matrix_type == "binary":
        V = df.astype(float).values
        W = np.ones_like(V)
    elif matrix_type == "ternary":
        V = df.replace(3, 0).astype(float).values
        W = (df != 3).astype(float).values

    # Determine K limits
    args_max_k = args.max_cluster
    max_k_snv = min(args_max_k, max(2, num_cells - 1)) if num_cells < args_max_k else args_max_k

    print("Performing broad hierarchical clustering...")
    broad_k = args.broad_k

    if num_cells <= broad_k:
         best_k = args.force_k if args.force_k else optimal_k_selection(V, W, args.min_k, max_k_snv)
         U, _ = weighted_nmf(V, W, best_k)
         labels = np.argmax(U, axis=1)
         clone_ids = np.ones(num_cells, dtype=int)
         narrow_clone_ids = [f"1_{x+1}" for x in labels]
    else:
        Z = fastcluster.linkage(V, method='ward', metric='euclidean')
        clone_ids = fcluster(Z, broad_k, criterion='maxclust')
        narrow_clone_ids = np.empty(num_cells, dtype=object)

        for broad_id in np.unique(clone_ids):
            clade_indices = np.where(clone_ids == broad_id)[0]
            if len(clade_indices) < args.min_k:
                narrow_clone_ids[clade_indices] = [f"{broad_id}_1"] * len(clade_indices)
                continue

            V_clade = V[clade_indices, :]
            W_clade = W[clade_indices, :]

            k_clade_max = min(len(clade_indices) - 1, max_k_snv)

            if args.force_k:
                k_clade = args.force_k
            else:
                k_clade = optimal_k_selection(V_clade, W_clade, args.min_k, k_clade_max)

            if k_clade < 1: k_clade = 1

            U_clade, _ = weighted_nmf(V_clade, W_clade, k_clade)
            clade_labels = np.argmax(U_clade, axis=1) + 1

            for i, narrow_id in enumerate(clade_labels):
                narrow_clone_ids[clade_indices[i]] = f"{broad_id}_{narrow_id}"

    out_tsv = scratch_output_dir / "clusters.tsv.gz"
    out_png = scratch_output_dir / "cluster.heatmap.png"

    out_df = pd.DataFrame({
        'cell_id': df.index.tolist(),
        'clone_id': clone_ids,
        'narrow_clone_id': narrow_clone_ids
    })
    out_df.to_csv(out_tsv, sep='\t', index=False, compression='gzip')
    print(f"Saved clusters to: {out_tsv}")

    try:
        print("Generating heatmap...")
        df_visual = pd.DataFrame(V, index=df.index, columns=df.columns)
        save_heatmap(df_visual, clone_ids, narrow_clone_ids, out_png, row_cluster=False)
        print(f"Saved heatmap to: {out_png}")
    except Exception as e:
        print(f"Heatmap generation failed: {e}")

    if args.final:
        _promote_files(args, out_tsv, out_png, MODULE_RESULTS_DIR, MODULE_PLOTS_DIR)

def _promote_files(args, out_tsv, out_png, MODULE_RESULTS_DIR, MODULE_PLOTS_DIR):
    print("\n--- Promoting Hierarchical Files to Module Directories ---")
    res_dest_dir = MODULE_RESULTS_DIR / args.project
    plot_dest_dir = MODULE_PLOTS_DIR / args.project
    if not res_dest_dir.exists(): os.makedirs(res_dest_dir)
    if not plot_dest_dir.exists(): os.makedirs(plot_dest_dir)

    dest_png = plot_dest_dir / f"{args.sample}.{args.cell_type}.clusters.hierarchical.png"
    if out_png.exists():
        shutil.copy(out_png, dest_png)
        print(f"Copied {out_png.name} -> {dest_png}")

    dest_tsv = res_dest_dir / f"{args.sample}.{args.cell_type}.clones.hierarchical.tsv.gz"
    if out_tsv.exists():
        shutil.copy(out_tsv, dest_tsv)
        print(f"Copied {out_tsv.name} -> {dest_tsv}")

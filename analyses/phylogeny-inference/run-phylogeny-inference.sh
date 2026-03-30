#!/bin/bash
# Eric Wafula
# 2025

# Set this so the whole loop stops if there is an error
set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

printf '\nStart TNBC sample data phylogeny inference...\n'

# Sample data - TNBC5
printf '\n-- Sample phylogeny inference...\n'
python3 01-phylogeny-inference.py \
  --project TNBC \
  --sample TNBC5 \
  --cell_type select-cell-types \
  --method iqtree \
  --threads 16
  
python3 01-phylogeny-inference.py \
  --project TNBC \
  --sample TNBC5 \
  --cell_type select-cell-types \
  --method fasttree \
  --threads 16  
  
  
printf '\n-- Filter VCF and matrices for manifold clustering.\n'
python3 02-filter-variants.py \
  --project TNBC \
  --sample TNBC5 \
  --cell_type select-cell-types  

printf '\n-- Sample manifold clonal clustering...\n'
# initial iteration to determine optimnal k
python3 03-snp-clustering.py \
  --project TNBC  \
  --sample TNBC5 \
  --cell_type select-cell-types \
  --method manifold \
  --max_cluster 10
  
# final interation  
python3 03-snp-clustering.py \
  --project TNBC  \
  --sample TNBC5 \
  --cell_type select-cell-types \
  --method manifold \
  --max_cluster 4 \
  --final

printf '\n-- Sample hierarchical clonal clustering...\n'
# initial iteration to determine optimnal k
python3 03-snp-clustering.py \
  --project TNBC \
  --sample TNBC5 \
  --cell_type select-cell-types \
  --method hierarchical \
  --max_cluste 10
  
# final interation  
python3 03-snp-clustering.py \
  --project TNBC \
  --sample TNBC5 \
  --cell_type select-cell-types \
  --method hierarchical \
  --broad_k 3 \
  --final  
  
printf '\n-- Sample annotations  and clonal clusters phylogeny maping...\n'

### Manifold clustering
# clonal clusters
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project TNBC \
  --tree_method fasttree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.5 \ # default
  --merge_clusters '0:1:2:3'

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project TNBC \
  --tree_method fasttree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell lineage annotations
  --cell_category TNBC-TNBC5-annotation_categories.tsv.gz # cell lineage annotations mappings 
  --concordance_clean \
  --min_clone_frac 0.5 \ # default
  --merge_clusters '0:1:2:3'
  
# clonal clusters
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project TNBC \
  --tree_method iqtree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.5 \ # default
  --merge_clusters '0:1:2:3'

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project TNBC \
  --tree_method iqtree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell lineage annotations
  --cell_category TNBC-TNBC5-annotation_categories.tsv.gz # cell lineage annotations mappings 
  --concordance_clean \
  --min_clone_frac 0.5 \ # default
  --merge_clusters '0:1:2:3'
  
### Hierarchical clustering
# clonal clusters
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project TNBC \
  --tree_method fasttree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.5 \ # default
  --merge_clusters '1:2:3'

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project TNBC \
  --tree_method fasttree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell lineage annotations
  --cell_category TNBC-TNBC5-annotation_categories.tsv.gz # cell lineage annotations mappings 
  --concordance_clean \
  --min_clone_frac 0.5 \ # default
  --merge_clusters '1:2:3'
  
# clonal clusters
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project TNBC \
  --tree_method iqtree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.5 \ # default
  --merge_clusters '1:2:3'

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project TNBC \
  --tree_method iqtree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell lineage annotations
  --cell_category TNBC-TNBC5-annotation_categories.tsv.gz # cell lineage annotations mappings 
  --concordance_clean \
  --min_clone_frac 0.5 \ # default
  --merge_clusters '1:2:3'
  
printf '\nStart SPECTRUM sample data phylogeny inference...\n'

# Sample data - MERGED 
printf '\n-- Sample phylogeny inference...\n'
python3 01-phylogeny-inference.py \
  --project SPECTRUM \
  --sample MERGED \
  --cell_type Ovarian-cancer \
  --method iqtree \
  --threads 16
  
python3 01-phylogeny-inference.py \
  --project SPECTRUM \
  --sample MERGED \
  --cell_type Ovarian-cancers \
  --method fasttree \
  --threads 16  
  
printf '\n-- Filter VCF and matrices for manifold clustering.\n'
python3 02-filter-variants.py \
  --project SPECTRUM \
  --sample MERGED \
  --cell_type Ovarian-cancer  

printf '\n-- Sample manifold clonal clustering...\n'
# initial iteration to determine optimnal k
python3 03-snp-clustering.py \
  --project SPECTRUM \
  --sample MERGED \
  --cell_type Ovarian-cancer \
  --method manifold \
  --max_cluster 10 
  
# final interation  
python3 03-snp-clustering.py \
  --project SPECTRUM \
  --sample MERGED \
  --cell_type Ovarian-cancer \
  --method manifold \
  --max_cluster 5 \
  --final
  
printf '\n-- Sample hierarchical clonal clustering...\n'
# initial iteration to determine optimnal k
python3 03-snp-clustering.py \
  --project SPECTRUM \
  --sample MERGED \
  --cell_type Ovarian-cancer \
  --method hierarchical \
  -max_cluster 10 
  
# final interation  
python3 03-snp-clustering.py \
  --project SPECTRUM \
  --sample MERGED \
  --cell_type Ovarian-cancer \
  --method hierarchical \
  --broad_k 5 \
  --final   
  
printf '\n-- Sample annotations  and clonal clusters phylogeny maping...\n'

### Manifold clustering
# clonal clusters
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project SPECTRUM \
  --tree_method fasttree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.5 \
  --merge_clusters '1:4' 

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project SPECTRUM \
  --tree_method fasttree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Sample \
  --concordance_clean \
  --min_clone_frac 0.5 \
  --merge_clusters '1:4'
  
# clonal clusters
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project SPECTRUM \
  --tree_method iqtree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.6 \
  --merge_clusters '1:4' 

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project SPECTRUM \
  --tree_method iqtree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Sample \
  --concordance_clean \
  --min_clone_frac 0.6 \
  --merge_clusters '1:4'
  
### Hierarchical clustering
# clonal clusters
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project SPECTRUM \
  --tree_method fasttree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.5 \
  --merge_clusters '1:4' 

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project SPECTRUM \
  --tree_method fasttree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Sample \
  --concordance_clean \
  --min_clone_frac 0.5 \
  --merge_clusters '1:4'
  
# clonal clusters
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project SPECTRUM \
  --tree_method iqtree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.5 \
  --merge_clusters '1:4' 

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project SPECTRUM \
  --tree_method iqtree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Sample \
  --concordance_clean \
  --min_clone_frac 0.5 \
  --merge_clusters '1:4'


printf '\nStart AML sample data phylogeny inference...\n'

# Sample data - LE1
printf '\n-- Sample phylogeny inference...\n'
python3 01-phylogeny-inference.py \
  --project AML \
  --sample LE1 \
  --cell_type all-cell-types \
  --method iqtree \
  --threads 16
  
python3 01-phylogeny-inference.py \
  --project AML \
  --sample LE1 \
  --cell_type all-cell-types \
  --method fasttree \
  --threads 16  
  
printf '\n-- Filter VCF and matrices for manifold clustering.\n'
python3 02-filter-variants.py \
  --project AML \
  --sample LE1 \
  --cell_type all-cell-types  

printf '\n-- Sample manifold clonal clustering...\n'
# initial iteration to determine optimnal k
python3 03-snp-clustering.py \
  --project AML \
  --sample LE1 \
  --cell_type all-cell-types \
  --method manifold \
  --max_cluster 10 
  
# final interation  
python3 03-snp-clustering.py \
  --project AML \
  --sample LE1 \
  --cell_type all-cell-types \
  --method manifold \
  --max_cluster 3 \
  --final
  
printf '\n-- Sample hierarchical clonal clustering...\n'
# initial iteration to determine optimnal k
python3 03-snp-clustering.py \
  --project AML \
  --sample LE1 \
  --cell_type all-cell-types \
  --method hierarchical \
  -max_cluster 10 
  
# final interation  
python3 03-snp-clustering.py \
  --project AML \
  --sample LE1 \
  --cell_type all-cell-types \
  --method hierarchical \
  --broad_k 3 \
  --final   
  
printf '\n-- Sample annotations  and clonal clusters phylogeny maping...\n'

### Manifold clustering
# clonal clusters
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method fasttree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.5 

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method fasttree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell lineage annotations
  --cell_category AML-LE1-annotation_categories.tsv.gz # cell lineage annotations mappings
  --concordance_clean \
  --min_clone_frac 0.5
  
# cell ploidy
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method fasttree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell ploidy inference 
  --cell_category AML-LE1-annotation_ploidy.tsv.gz # cell ploidy inference mappings
  --concordance_clean \
  --min_clone_frac 0.5
  
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method iqtree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.5 

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method iqtree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell lineage annotations
  --cell_category AML-LE1-annotation_categories.tsv.gz # cell lineage annotations mappings
  --concordance_clean \
  --min_clone_frac 0.5
  
# cell ploidy
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method iqtree \
  --cluster_method manifold \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell ploidy inference 
  --cell_category AML-LE1-annotation_ploidy.tsv.gz # cell ploidy inference mappings
  --concordance_clean \
  --min_clone_frac 0.5
  
### Hierarchical clustering
# clonal clusters
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method fasttree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.5 

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method fasttree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell lineage annotations
  --cell_category AML-LE1-annotation_categories.tsv.gz # cell lineage annotations mappings
  --concordance_clean \
  --min_clone_frac 0.5
  
# cell ploidy
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method fasttree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell ploidy inference 
  --cell_category AML-LE1-annotation_ploidy.tsv.gz # cell ploidy inference mappings
  --concordance_clean \
  --min_clone_frac 0.5
  
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method iqtree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Clone \
  --concordance_clean \
  --min_clone_frac 0.5 

# cell type lineages
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method iqtree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell lineage annotations
  --cell_category AML-LE1-annotation_categories.tsv.gz # cell lineage annotations mappings
  --concordance_clean \
  --min_clone_frac 0.5
  
# cell ploidy
Rscript --vanilla 04-infer-clonal-phylogeny.R \
  --project AML \
  --tree_method iqtree \
  --cluster_method hierarchical \
  --refseq \
  --rescale \
  --annot_type Lineage \ # cell ploidy inference 
  --cell_category AML-LE1-annotation_ploidy.tsv.gz # cell ploidy inference mappings
  --concordance_clean \
  --min_clone_frac 0.5

# remove Rplots files created ggtree
# need to figure out avoid this unintended pdfs being created by ggtree
rm -rf Rplots.pdf  

printf '\nPhylogeny inference Done...\n'

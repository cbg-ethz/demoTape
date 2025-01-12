#!/usr/bin/env Rscript

library(argparser)

p <- arg_parser('Run cloneAlign.')
p <- add_argument(p, '--exp', help='Input expression file (str)')
p <- add_argument(p, '--cnv', help='Input copy number file (str)')
p <- add_argument(p, '--mcpc', default=100,
    help = 'Minimal counts per cell, otherwise removed cell (int)')
p <- add_argument(p, '--out_file', help='output log file (str)')

argv <- parse_args(p)

library(tensorflow)
use_condaenv("r_env")
library(clonealign)

# library(tibble)
#
# Rewrite compute_correlations function
# compute_correlations <- function(Y, L, clones) {
#   unassigned <- clones == "unassigned"
#   Y <- Y[!unassigned,]
#  
#   ## Scale Y expression
#   Y <- scale(Y)
#  
#   clones <- clones[!unassigned]
#  
#   print(as_tibble(Y))
#   print(as_tibble(L))
#   print(as_tibble(clones))
#
#   sapply(seq_len(ncol(Y)), function(i) {
#     y <- Y[,i]
#     x <- L[i, clones]
#
#     suppressWarnings({
#       cor(x,y)
#     })
#   })
# }
# assignInNamespace('compute_correlations', compute_correlations, getNamespace('clonealign'))


exp_data <- t(as.matrix(read.csv(argv$exp, row.names=1, header=TRUE)))
cnv_data <- as.matrix(read.csv(argv$cnv, row.names=1, header=TRUE))

cnv_data[cnv_data == 0] <- 0.1

ca_data <- preprocess_for_clonealign(
    exp_data,
    cnv_data,
    min_counts_per_gene = 20,
    min_counts_per_cell=argv$mcpc, # 100
    remove_outlying_genes=TRUE,
    nmads=10,
    max_copy_number=6,
    remove_genes_same_copy_number=TRUE # TRUE
)

cal <- run_clonealign(
    ca_data$gene_expression_data,
    ca_data$copy_number_data,
    initial_shrinks = c(0.1, 1, 5, 10), # crashes for 0 for 1 samples for whatever reason
    n_repeats=3,
)

df <- as.data.frame(
    cal$clone,
    row.names=rownames(ca_data$gene_expression_data),
    nm="barcode"
)

write.table(
    df,
    file=argv$out_file,
    sep="\t",
    row.names=TRUE,
)
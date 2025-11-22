args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript fig6_simple_tree.R input_treefile output_pdf")
}

in_tree <- args[1]
out_pdf <- args[2]

suppressPackageStartupMessages(library(ape))

tree <- read.tree(in_tree)

# ladderize for a nice left-to-right shape
tree <- ladderize(tree)

pdf(out_pdf, width = 6, height = 7)
par(mar = c(1, 1, 1, 1))
plot(tree,
     cex = 0.5,        # tip label size
     no.margin = TRUE)
add.scale.bar()
dev.off()
cat("Saved simple tree figure to", out_pdf, "\n")

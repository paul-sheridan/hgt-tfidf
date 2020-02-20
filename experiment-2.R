########################################################################################################################
#
# Description: This R script compares tp, tpidf, and the hypergeometric test on NYSK documents in the document
#  summarization scenario.
#
# Variables of note:
#  term[i] = i'th term in vocabulary
#  doc[[j]] = j'th doc in collection (unique terms)
#  N_bar = number of documents in the collection
#  N = total number of terms in the collection (including multiplicities) 
#  n_bar[j] = number of unique terms in doc[[j]]
#  n[j] = number of terms in doc[[j]] (including multiplicities)
#  K_bar[i] = number of docs in which term[i] occurs at least once
#  K[i] = number of occurrences of term[i] in the collection
#  k[[j]][i] = number of occurrences of term[i] in doc[[j]]
#  T = number of unique terms in the collection
#  tp_score_lst[[j]][i] = tp score of term[i] in doc[[j]]
#  tpidf_score_lst[[j]][i] = tpidf score of term[i] in doc[[j]]
#  hgt_score_lst[[j]][i] = hypergeometric test score of term[i] in doc[[j]]
#
########################################################################################################################


## preliminaries
rm(list = ls()) # deletes any already initialized variables in R session
library(ggplot2)
library(magicaxis)
library(rjson)
working_dir <- getwd()
input_dir <- paste0(working_dir, "/stats/")
plots_output_dir <- paste0(working_dir, "/plots/")
setwd(working_dir)
set.seed(63429) # initialize random seed


## set plot colours 
gg_color_hue = function(n, alpha = 1) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1 : n]
}
cols <- c("grey25", gg_color_hue(n = 3))
tcols <- c("grey25", gg_color_hue(n = 3, alpha = 0.5))
gray25 <- cols[1]; red <- cols[2]; green <- cols[3]; blue <- cols[4]
tred <- cols[2]; tgreen <- tcols[3]; tblue <- tcols[4]


## load some NYSK data summary statistics
load(paste0(input_dir, "N_bar-uppercase.Rdata"))
load(paste0(input_dir, "K_bar-uppercase.Rdata"))
load(paste0(input_dir, "n_bar-lowercase.Rdata"))
load(paste0(input_dir, "doc.Rdata"))
load(paste0(input_dir, "tp_score_lst.Rdata"))
load(paste0(input_dir, "tpidf_score_lst.Rdata"))
load(paste0(input_dir, "hgt_score_lst.Rdata"))


## find over-represented terms in documents overlap
top <- 10
M <- max(K_bar)
p10_mat <- mat.or.vec(nr = N_bar, nc = 4)
colnames(p10_mat) <- c("HGT-RND", "HGT-TP", "HGT-TPIDF", "TP-TPIDF")

for (j in 1 : N_bar) {
  myterms <- doc[[j]]
  tp_ranked_terms <- myterms[sort(rank(-tp_score_lst[[j]], ties.method = "random"), index.return = TRUE)$ix]
  tpidf_ranked_terms <- myterms[sort(rank(-tpidf_score_lst[[j]], ties.method = "random"), index.return = TRUE)$ix]
  hgt_ranked_terms <- myterms[sort(rank(-hgt_score_lst[[j]], ties.method = "random"), index.return = TRUE)$ix]

  if(n_bar[j] <= top) {
    p10_mat[j, "HGT-RND"] <- top
  } else {
    p10_mat[j, "HGT-RND"] <- top * top / n_bar[j]
  }

  p10_mat[j, "HGT-TP"] <- length(intersect(hgt_ranked_terms[1 : top], tp_ranked_terms[1 : top]))
  p10_mat[j, "HGT-TPIDF"] <- length(intersect(hgt_ranked_terms[1 : top], tpidf_ranked_terms[1 : top]))
  p10_mat[j, "TP-TPIDF"] <- length(intersect(tp_ranked_terms[1 : top], tpidf_ranked_terms[1 : top]))
}

apply(p10_mat, 2, mean)
# > apply(p10_mat, 2, mean)
#   HGT-RND    HGT-TP HGT-TPIDF  TP-TPIDF 
# 0.5837463 4.1779100 8.4683812 3.1569907

apply(p10_mat, 2, sd)
# > apply(p10_mat, 2, sd)
#   HGT-RND    HGT-TP HGT-TPIDF  TP-TPIDF 
# 0.7120265 1.6838443 1.0380862 1.6763330


## generate overlapping histogram plot
myhist <- hist(
  p10_mat[,"HGT-RND"],
  breaks = seq(from = 0, to = top, by = 1),
  plot = FALSE)
counts1 <- myhist$counts

myhist <- hist(
  p10_mat[,"HGT-TPIDF"],
  breaks = seq(from = 0, to = top, by = 1),
  plot = FALSE)
counts2 <- myhist$counts

myhist <- hist(
  p10_mat[,"HGT-TP"],
  breaks = seq(from = 0, to = top, by = 1),
  plot = FALSE)
counts3 <- myhist$counts

y1 <- 0
y2 <- max(c(counts1, counts2, counts3))

pdf(paste0(plots_output_dir, "figure-1b.pdf"), width = 6, height = 6)
plot(
  x = counts1,
  main = "Document Summarization Scenario",
  xlab = "P@10 Score",
  ylab = "Number of Documents",
  log = "x",
  axes = FALSE,
  type = "n")
magaxis()

for (i in 1 : top) {
  rect(xleft = i - 1 + 0.5, ybottom = 0, xright = i + 0.5, ytop = counts1[i], col = tred)
}

for (i in 1 : top) {
  rect(xleft = i - 1 + 0.5, ybottom = 0, xright = i + 0.5, ytop = counts2[i], col = tblue)
}

for (i in 1 : top) {
  rect(xleft = i - 1 + 0.5, ybottom = 0, xright = i + 0.5, ytop = counts3[i], col = tgreen)
}

legend( 
  "topright", 
  c("Hypergeometric Test / Random", "Hypergeometric Test / TP", "Hypergeometric Test / TPIDF"), 
  fill = c(tred, tgreen, tblue),
  cex = 0.9,
  inset = 0.03,
  box.lwd = 0,
  box.col = "grey90",
  bg = "grey90", 
  seg.len = 4)
dev.off()

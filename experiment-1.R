########################################################################################################################
#
# Description: This R script compares tp, tpidf, and the hypergeometric test on NYSK documents in the one-term query
#  document retrieval scenario
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
working_dir <- getwd()
input_dir <- paste0(working_dir, "/stats/")
stats_output_dir <- paste0(working_dir, "/experiment-2-stats/")
plots_output_dir <- paste0(working_dir, "/plots/")
setwd(working_dir)
set.seed(24329) # initialize random seed


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
load(paste0(input_dir, "T-uppercase.Rdata"))
load(paste0(input_dir, "N_bar-uppercase.Rdata"))
load(paste0(input_dir, "K_bar-uppercase.Rdata"))
load(paste0(input_dir, "n_bar-lowercase.Rdata"))
load(paste0(input_dir, "doc.Rdata"))
load(paste0(input_dir, "term_to_doc_ids_lst.Rdata"))
load(paste0(input_dir, "term.Rdata"))


## calculate P@10 scores for each individual terms in the collection
p10_mat <-  mat.or.vec(nr = T, nc = 4)
colnames(p10_mat) <- c("HGT-RND", "HGT-TP", "HGT-TPIDF", "TP-TPIDF")

# tp and tp-tpidf
load(paste0(input_dir, "tp_score_mat.Rdata"))
load(paste0(input_dir, "tpidf_score_mat.Rdata"))

for (i in 1 : T) {
  doc_ids <- term_to_doc_ids_lst[[i]]
  tp_ranked_docs <- sort(rank(-tp_score_mat[i, doc_ids], ties.method = "random"), index.return = TRUE)$ix
  tpidf_ranked_docs <- sort(rank(-tpidf_score_mat[i, doc_ids], ties.method = "random"), index.return = TRUE)$ix
  
  if (K_bar[i] <= top) {
    p10_mat[i, "TP-TPIDF"] <- top 
  } else {
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }
}

rm(tpidf_score_mat)

# tp and hgt
load(paste0(input_dir, "hgt_score_mat.Rdata"))

for (i in 1 : T) {
  doc_ids <- term_to_doc_ids_lst[[i]]
  tp_ranked_docs <- sort(rank(-tp_score_mat[i, doc_ids], ties.method = "random"), index.return = TRUE)$ix
  hgt_ranked_docs <- sort(rank(-hgt_score_mat[i, doc_ids], ties.method = "random"), index.return = TRUE)$ix
  
  if (K_bar[i] <= top) {
    p10_mat[i, "HGT-TP"] <- top 
  } else {
    p10_mat[i, "HGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
  }
}

rm(tp_score_mat)

# tpidf and hgt
load(paste0(input_dir, "tpidf_score_mat.Rdata"))

for (i in 1 : T) {
  doc_ids <- term_to_doc_ids_lst[[i]]
  tpidf_ranked_docs <- sort(rank(-tpidf_score_mat[i, doc_ids], ties.method = "random"), index.return = TRUE)$ix
  hgt_ranked_docs <- sort(rank(-hgt_score_mat[i, doc_ids], ties.method = "random"), index.return = TRUE)$ix
  
  if (K_bar[i] <= top) {
    p10_mat[i, "HGT-TPIDF"] <- top 
  } else {
    p10_mat[i, "HGT-TPIDF"] <- length(intersect(tfidf_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
  }
}

rm(tpidf_score_mat)

# hgt and random
for (i in 1 : T) {
  doc_ids <- term_to_doc_ids_lst[[i]]
  hgt_ranked_docs <- sort(rank(-hgt_score_mat[i, doc_ids], ties.method = "random"), index.return = TRUE)$ix
  
  if (K_bar[i] <= top) {
    p10_mat[i, "HGT-RND"] <- top 
  } else {
    p10_mat[i, "HGT-RND"] <- top * top / K_bar[i]
  }
}

rm(hgt_score_mat)

apply(p10_mat, 2, mean)
# > apply(p10_mat, 2, mean)
#   HGT-RND    HGT-TP HGT-TPIDF  TP-TPIDF 
#  8.684827  9.509527  9.509283  9.984217 

save(p10_mat, file = paste0(stats_output_dir, "p10_mat.Rdata"))


## P@10 scores by cutoff C
M <- max(K_bar)
p10_mat_at_C <-  mat.or.vec(nr = M, nc = 4)
colnames(p10_mat_at_C) <- c("HGT-RND", "HGT-TP", "HGT-TPIDF", "TP-TPIDF")

for (C in 1 : M) {
  indices <- which(K_bar >= C)
  p10_mat_at_C[C, "HGT-RND"] <- mean(p10_mat[indices, "HGT-RND"])
  p10_mat_at_C[C, "HGT-TP"] <- mean(p10_mat[indices, "HGT-TP"])
  p10_mat_at_C[C, "HGT-TPIDF"] <- mean(p10_mat[indices, "HGT-TPIDF"])
  p10_mat_at_C[C, "TP-TPIDF"] <- mean(p10_mat[indices, "TP-TPIDF"])
}

save(p10_mat_at_C, file = paste0(stats_output_dir, "p10_mat_at_C.Rdata"))

## generate cutoff plot
x <- 1 : M
y1 <- p10_mat_at_C[, "TP-TPIDF"]
y2 <- p10_mat_at_C[, "HGT-TPIDF"]
y3 <- p10_mat_at_C[, "HGT-RND"]

pdf(paste0(plots_output_dir, "figure-1a.pdf"), width = 6, height = 6)
plot(x = x,
  y = y1,
  xlab = "Cutoff Value, C",
  ylab = "Average P@10 Score",
  main = "Document Retrieval Scenario",
  ylim = c(0, 10),
  type = "n",
  log = "x",
  axes = FALSE)
grid(lwd = 2)
magaxis()
points(x = x,
  y = y1,
  type = "l",
  lty = "dashed",
  lwd = 2,
  col = gray25)
points(x = x,
  y = y2,
  type = "l",
  lwd = 2,
  col = gray25)
points(x = x,
  y = y3,
  type = "l",
  lty = "dotted",
  lwd = 3,
  col = gray25)
legend( 
  "topright", 
  c("TP / TPIDF", "Hypergeometric Test / TPIDF", "Hypergeometric Test / Random"), 
  col = c(gray25, gray25), 
  lty = c("dashed", "solid", "dotted"), 
  lwd = c(2, 2, 3), 
  cex = 0.9,
  inset = 0.08,
  box.lwd = 0,
  box.col = "grey90",
  bg = "grey90", 
  seg.len = 4)
dev.off()

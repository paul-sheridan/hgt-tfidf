########################################################################################################################
#
# Description: This R script calculates some basic statistics about the NYSK collection documents.
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
library(rjson)
working_dir <- getwd()
output_dir <- paste0(working_dir, "/stats/")
setwd(working_dir)
set.seed(12345) # initialize random seed


## process NYSK collection documents and calculate basic statistics
infile <- "nysk-processed.json"
mydata <- fromJSON(paste(readLines(infile), collapse = ""))
raw_docs <- mydata
N_bar <- length(raw_docs)
k <- vector("list", length = N_bar)
doc <- vector("list", length = N_bar)

for (j in 1 : N_bar) {
  mydoc <- raw_docs[[j]]
  mytable <- sort(table(mydoc), decreasing = T)
  k[[j]] <- as.numeric(mytable)
  doc[[j]] <- names(mytable)
}

n <- unlist(lapply(k, sum))
n_bar <- unlist(lapply(k, length))
term <- unique(unlist(doc))
T <- length(term)
N <- sum(unlist(k))
K <- table(unlist(raw_docs))[term]
K_bar <- table(unlist(doc))[term]
save(doc, file = paste0(output_dir, "doc.Rdata"))
save(k, file = paste0(output_dir, "k-lowercase.Rdata"))
save(n, file = paste0(output_dir, "n-lowercase.Rdata"))
save(n_bar, file = paste0(output_dir, "n_bar-lowercase.Rdata"))
save(term, file = paste0(output_dir, "term.Rdata"))
save(T, file = paste0(output_dir, "T-uppercase.Rdata"))
save(N, file = paste0(output_dir, "N-uppercase.Rdata"))
save(N_bar, file = paste0(output_dir, "N_bar-uppercase.Rdata"))
save(K, file = paste0(output_dir, "K-uppercase.Rdata"))
save(K_bar, file = paste0(output_dir, "K_bar-uppercase.Rdata"))


## calculate tp/tpidf/hgt term score lists
tp_score_lst <- vector("list", length = N_bar)
tpidf_score_lst <- vector("list", length = N_bar)
hgt_score_lst <- vector("list", length = N_bar)

for (j in 1 : N_bar) {
  tp_score_lst[[j]] <- numeric(n_bar[j])
  tpidf_score_lst[[j]] <- numeric(n_bar[j])
  hgt_score_lst[[j]] <- numeric(n_bar[j])

  for (i in 1 : n_bar[j]) {
    tp_score_lst[[j]][i] <- k[[j]][i] / n[j]
    tpidf_score_lst[[j]][i] <- - k[[j]][i] / n[j] * log(K_bar[doc[[j]][i]] / N_bar)
    hgt_score_lst[[j]][i] <- - phyper(k[[j]][i] - 1, K[doc[[j]][i]], N - K[doc[[j]][i]], n[j], lower.tail = FALSE, log.p = TRUE)
  }
}

save(tp_score_lst, file = paste0(output_dir, "tp_score_lst.Rdata"))
save(tpidf_score_lst, file = paste0(output_dir, "tpidf_score_lst.Rdata"))
save(hgt_score_lst, file = paste0(output_dir, "hgt_score_lst.Rdata"))


## calculate list of which docs each term occurs in
term_to_doc_ids_lst <- vector("list", length = T)
names(term_to_doc_ids_lst) <- term

for (i in 1 : T) {
  if (i %% 10000 == 0) { print(i) }
  count <- 1
  term_to_doc_ids <- numeric(K_bar[i])
  
  for (j in 1 : N_bar) {
    if (term[i] %in% doc[[j]]) {
      term_to_doc_ids[count] <- j
      count <- count + 1
    }
  }

  term_to_doc_ids_lst[[i]] <- term_to_doc_ids
}
save(term_to_doc_ids_lst, file = paste0(output_dir, "term_to_doc_ids_lst.Rdata"))


## calculate tp term score by document matrix
tp_score_mat <- mat.or.vec(nr = T, nc = N_bar)
rownames(tp_score_mat) <- term

for (i in 1 : T) {
  for (l in 1 : K_bar[i]) {
    j <- term_to_doc_ids_lst[[i]][l]
    tp_score_mat[i, j] <- tp_score_lst[[j]][which(doc[[j]] == term[i])]
  }
}
save(tp_score_mat, file = paste0(output_dir, "tp_score_mat.Rdata"))

## calculate term burstiness scores
term_burstiness_score <- apply(tp_score_mat, 1, sum) / K_bar
save(term_burstiness_score, file = paste0(output_dir, "term_burstiness_score.Rdata"))


## calculate tpidf term score by document matrix
tpidf_score_mat <- mat.or.vec(nr = T, nc = N_bar)
rownames(tpidf_score_mat) <- term

for (i in 1 : T) {
  for (l in 1 : K_bar[i]) {
    j <- term_to_doc_ids_lst[[i]][l]
    tpidf_score_mat[i, j] <- tpidf_score_lst[[j]][which(doc[[j]] == term[i])]
  }
}
save(tpidf_score_mat, file = paste0(output_dir, "tpidf_score_mat.Rdata"))


## calculate hypergeometric test term score by document matrix
hgt_score_mat <- mat.or.vec(nr = T, nc = N_bar)
rownames(hgt_score_mat) <- term

for (i in 1 : T) {
  for (l in 1 : K_bar[i]) {
    j <- term_to_doc_ids_lst[[i]][l]
    hgt_score_mat[i, j] <- hgt_score_lst[[j]][which(doc[[j]] == term[i])]
  }
}
save(hgt_score_mat, file = paste0(output_dir, "hgt_score_mat.Rdata"))

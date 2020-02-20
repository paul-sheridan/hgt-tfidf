########################################################################################################################
#
# Description: This R script compares tp, tpidf, and the hypergeometric test on NYSK documents in the two-term query
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
working_dir <- getwd()
input_dir <- paste0(working_dir, "/stats/")
stats_output_dir <- paste0(working_dir, "/experiment-3-stats/")
setwd(working_dir)
set.seed(24329) # initialize random seed


## load some NYSK data summary statistics
load(paste0(input_dir, "T-uppercase.Rdata"))
load(paste0(input_dir, "N_bar-uppercase.Rdata"))
load(paste0(input_dir, "K_bar-uppercase.Rdata"))
load(paste0(input_dir, "n_bar-lowercase.Rdata"))
load(paste0(input_dir, "doc.Rdata"))
load(paste0(input_dir, "term_to_doc_ids_lst.Rdata"))
load(paste0(input_dir, "term_burstiness_score.Rdata"))


## calculate top 106 most bursty terms and load term scores
top <- 10
cut <- 106
no_of_queries <- cut - top
top_bursty_score <- sort(term_burstiness_score[which(K_bar >= 10)], decreasing = TRUE)[1 : cut]
top_bursty_term <- names(top_bursty_score)
load(paste0(input_dir, "tp_score_mat.Rdata"))
tp_score_submat <- tp_score_mat[top_bursty_term, ]
rm(tp_score_mat)
load(paste0(input_dir, "tpidf_score_mat.Rdata"))
tpidf_score_submat <- tpidf_score_mat[top_bursty_term, ]
rm(tpidf_score_mat)
load(paste0(input_dir, "hgt_score_mat.Rdata"))
hgt_score_submat <- hgt_score_mat[top_bursty_term, ]
rm(hgt_score_mat)


## retrieve the three bursty terms with maximum K_bar values
top_six_bursty_terms <- names(sort(K_bar[top_bursty_term], decreasing = TRUE))[1 : 6]
# > sort(K_bar[top_bursty_term], decreasing = TRUE)[1 : 6]
# strausskahn         say         new         imf     comment     lagarde 
#        9634        8891        8666        8094        4106        2376

## empty rows ids
empty_row_id <- numeric(length(top_six_bursty_terms))
for (i in 1 : length(top_six_bursty_terms)) {
  empty_row_id[i] <- which(top_bursty_term == top_six_bursty_terms[i])
}


## calculate P@10 scores for all queries q = ("strausskahn", t)
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("HGT-RND", "HGT-TP", "HGT-TPIDF", "TP-TPIDF")
t1 <- "strausskahn"

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    hgt_ranked_docs <- sort(rank(-apply(hgt_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    p10_mat[i, "HGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "HGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "HGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }
}

save(hgt_score_submat, file = paste0(stats_output_dir, "hgt_score_submat_strausskahn.Rdata"))
save(p10_mat) file = paste0(stats_output_dir, "p10_mat_strausskahn.Rdata"))

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#    HGT-RND     HGT-TP  HGT-TPIDF   TP-TPIDF 
# 0.01037412 0.29000000 2.90000000 0.63000000
apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#      HGT-RND       HGT-TP    HGT-TPIDF     TP-TPIDF 
# 1.401856e-05 1.335566e+00 3.520503e+00 1.599590e+00


## calculate P@10 scores for all queries q = ("say", t)
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("HGT-RND", "HGT-TP", "HGT-TPIDF", "TP-TPIDF")
t1 <- "say"

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    hgt_ranked_docs <- sort(rank(-apply(hgt_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    p10_mat[i, "HGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "HGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "HGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }
}

save(hgt_score_submat, file = paste0(stats_output_dir, "hgt_score_submat_say.Rdata"))
save(p10_mat, file = paste0(stats_output_dir, "p10_mat_say.Rdata"))

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#    HGT-RND     HGT-TP  HGT-TPIDF   TP-TPIDF 
# 0.01123226 1.86000000 3.91000000 0.83000000
apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#      HGT-RND       HGT-TP    HGT-TPIDF     TP-TPIDF 
# 3.114928e-05 1.477235e+00 3.648979e+00 1.896595e+00


## calculate P@10 scores for all queries q = ("new", t)
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("HGT-RND", "HGT-TP", "HGT-TPIDF", "TP-TPIDF")
t1 <- "new"

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    hgt_ranked_docs <- sort(rank(-apply(hgt_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    p10_mat[i, "HGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "HGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "HGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }
}

save(hgt_score_submat, file = paste0(stats_output_dir, "hgt_score_submat_new.Rdata"))
save(p10_mat, file = paste0(stats_output_dir, "p10_mat_new.Rdata"))

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#    HGT-RND     HGT-TP  HGT-TPIDF   TP-TPIDF 
# 0.01152406 0.25000000 2.13000000 0.79000000
apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#      HGT-RND       HGT-TP    HGT-TPIDF     TP-TPIDF 
# 3.745847e-05 1.336171e+00 2.987026e+00 1.816284e+00


## calculate P@10 scores for all queries q = ("imf", t)
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("HGT-RND", "HGT-TP", "HGT-TPIDF", "TP-TPIDF")
t1 <- "imf"

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    hgt_ranked_docs <- sort(rank(-apply(hgt_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    p10_mat[i, "HGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "HGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "HGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }
}

save(hgt_score_submat, file = paste0(stats_output_dir, "hgt_score_submat_imf.Rdata"))
save(p10_mat, file = paste0(stats_output_dir, "p10_mat_imf.Rdata"))

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#    HGT-RND     HGT-TP  HGT-TPIDF   TP-TPIDF 
# 0.01233398 0.17000000 3.38000000 1.17000000
apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#      HGT-RND       HGT-TP    HGT-TPIDF     TP-TPIDF 
# 5.185109e-05 1.128644e+00 3.778675e+00 2.261022e+00


## calculate P@10 scores for all queries q = ("comment", t)
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("HGT-RND", "HGT-TP", "HGT-TPIDF", "TP-TPIDF")
t1 <- "comment"

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    hgt_ranked_docs <- sort(rank(-apply(hgt_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    p10_mat[i, "HGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "HGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "HGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }
}

save(hgt_score_submat, file = paste0(stats_output_dir, "hgt_score_submat_comment.Rdata"))
save(p10_mat, file = paste0(stats_output_dir, "p10_mat_comment.Rdata"))

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#    HGT-RND     HGT-TP  HGT-TPIDF   TP-TPIDF 
# 0.02419001 0.20000000 1.08000000 4.81000000
apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#      HGT-RND       HGT-TP    HGT-TPIDF     TP-TPIDF 
# 0.0005321128 1.1891767800 2.2904280399 2.6655111133


## calculate P@10 scores for all queries q = ("lagarde", t)
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("HGT-RND", "HGT-TP", "HGT-TPIDF", "TP-TPIDF")
t1 <- "lagarde"

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    hgt_ranked_docs <- sort(rank(-apply(hgt_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    p10_mat[i, "HGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "HGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "HGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], hgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }
}

save(hgt_score_submat, file = paste0(stats_output_dir, "hgt_score_submat_lagarde.Rdata"))
save(p10_mat, file = paste0(stats_output_dir, "p10_mat_lagarde.Rdata"))

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#    HGT-RND     HGT-TP  HGT-TPIDF   TP-TPIDF 
# 0.04131877 0.40000000 1.34000000 4.58000000
apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#     HGT-RND      HGT-TP   HGT-TPIDF    TP-TPIDF 
# 0.001907695 1.530827993 2.161275344 3.318695131

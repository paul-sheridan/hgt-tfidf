########################################################################################################################
#
# Description: This R script compares tp, tpidf, and the multivariate hypergeometric test on NYSK documents in 
#  the two-term query document retrieval scenario
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
stats_output_dir <- paste0(working_dir, "/experiment-4-stats/")
setwd(working_dir)
library("extraDistr")
source("mvhgt-test.R")
set.seed(24329) # initialize random seed


## load some NYSK data summary statistics
load(paste0(input_dir, "T-uppercase.Rdata"))
load(paste0(input_dir, "N-uppercase.Rdata"))
load(paste0(input_dir, "K-uppercase.Rdata"))
load(paste0(input_dir, "k-lowercase.Rdata"))
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


## calculate P@10 scores for all queries q = ("strausskhan", t)
count <- 0
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("MVHGT-RND", "MVHGT-TP", "MVHGT-TPIDF", "TP-TPIDF")
t1 <- "strausskahn"
n1 <- as.numeric(K[t1])
mvhgt_score_submat <- mat.or.vec(nr = cut, nc = N_bar)

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]
  n2 <- as.numeric(K[t2])
  n3 <- N - n1 - n2
  n <- c(n1, n2, n3)
  if (i %% 1 == 0) { print(paste0("i = ", i, ", t2 = ", t2)) }

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix

    for (j in 1 : N_bar) {
      if (j %% 1000 == 0) { print(paste0("i = ", i, " & j = ", j)) }
      p <- k[[j]]
      names(p) <- doc[[j]]
      if ((t1 %in% doc[[j]]) == TRUE) {
        x1 <- as.numeric(p[t1])
      } else {
        x1 <- 0
      }
      if ((t2 %in% doc[[j]]) == TRUE) {
        x2 <- as.numeric(p[t2])
      } else {
        x2 <- 0
      }
      x3 <- sum(p) - x1 - x2
      x <- c(x1, x2, x3)

      mvhgt_score_submat[i, j] <- pmvhyper(x = x, n = n, k = sum(p), log = TRUE)

      if (mvhgt_score_submat[i, j] == Inf) {
        mvhgt_score_submat[i, j] <- pmvhyper3(x = x, n = n, k = sum(p), log = TRUE)
        count <- count + 1
        print(paste0("(i,j) = (", i, ",", j, "), c = ", count, ", -log(P) = ", mvhgt_score_submat[i, j], ", x = (", x1, ",", x2, ",", x3, ")", ", n = (", n1, ",", n2, ",", n3, ")"))
      }
    }

    mvhgt_ranked_docs <- sort(rank(-mvhgt_score_submat[i, ], ties.method = "random"), index.return = TRUE)$ix

    p10_mat[i, "MVHGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "MVHGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "MVHGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }

  save(mvhgt_score_submat, file = paste0(output_dir, "mvhgt_score_submat_strausskahn.Rdata"))
  save(p10_mat, file = paste0(output_dir, "p10_mat_strausskahn.Rdata"))
}

## number of infinite scores
length(which(mvhgt_score_submat == Inf))
# > length(which(mvhgt_score_submat == Inf))
# [1] 0

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#   MVHGT-RND    MVHGT-TP MVHGT-TPIDF    TP-TPIDF 
#  0.01037412  0.30000000  2.36000000  0.62000000
apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#    MVHGT-RND     MVHGT-TP  MVHGT-TPIDF     TP-TPIDF 
# 1.401856e-05 1.306549e+00 2.779870e+00 1.549063e+00


## calculate P@10 scores for all queries q = ("say", t)
count <- 0
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("MVHGT-RND", "MVHGT-TP", "MVHGT-TPIDF", "TP-TPIDF")
t1 <- "say"
n1 <- as.numeric(K[t1])
mvhgt_score_submat <- mat.or.vec(nr = cut, nc = N_bar)

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]
  n2 <- as.numeric(K[t2])
  n3 <- N - n1 - n2
  n <- c(n1, n2, n3)
  if (i %% 1 == 0) { print(paste0("i = ", i, ", t2 = ", t2)) }

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix

    for (j in 1 : N_bar) {
      if (j %% 1000 == 0) { print(paste0("i = ", i, " & j = ", j)) }
      p <- k[[j]]
      names(p) <- doc[[j]]
      if ((t1 %in% doc[[j]]) == TRUE) {
        x1 <- as.numeric(p[t1])
      } else {
        x1 <- 0
      }
      if ((t2 %in% doc[[j]]) == TRUE) {
        x2 <- as.numeric(p[t2])
      } else {
        x2 <- 0
      }
      x3 <- sum(p) - x1 - x2
      x <- c(x1, x2, x3)

      mvhgt_score_submat[i, j] <- pmvhyper(x = x, n = n, k = sum(p), log = TRUE)

      if (mvhgt_score_submat[i, j] == Inf) {
        mvhgt_score_submat[i, j] <- pmvhyper3(x = x, n = n, k = sum(p), log = TRUE)
        count <- count + 1
        print(paste0("(i,j) = (", i, ",", j, "), c = ", count, ", -log(P) = ", mvhgt_score_submat[i, j], ", x = (", x1, ",", x2, ",", x3, ")", ", n = (", n1, ",", n2, ",", n3, ")"))
      }
    }

    mvhgt_ranked_docs <- sort(rank(-mvhgt_score_submat[i, ], ties.method = "random"), index.return = TRUE)$ix

    p10_mat[i, "MVHGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "MVHGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "MVHGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }

  save(mvhgt_score_submat, file = paste0(output_dir, "mvhgt_score_submat_say.Rdata"))
  save(p10_mat, file = paste0(output_dir, "p10_mat_say.Rdata"))
}

## number of infinite scores
length(which(mvhgt_score_submat == Inf))
# > length(which(mvhgt_score_submat == Inf))
# [1] 0

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#   MVHGT-RND    MVHGT-TP MVHGT-TPIDF    TP-TPIDF 
#  0.01123226  1.72000000  3.20000000  0.82000000
apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#    MVHGT-RND     MVHGT-TP  MVHGT-TPIDF     TP-TPIDF 
# 3.114928e-05 1.146624e+00 2.902455e+00 1.855268e+00


## calculate p@10 scores for all queries q = ("new", t)
count <- 0
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("MVHGT-RND", "MVHGT-TP", "MVHGT-TPIDF", "TP-TPIDF")
t1 <- "new"
n1 <- as.numeric(K[t1])
mvhgt_score_submat <- mat.or.vec(nr = cut, nc = N_bar)

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]
  n2 <- as.numeric(K[t2])
  n3 <- N - n1 - n2
  n <- c(n1, n2, n3)
  if (i %% 1 == 0) { print(paste0("i = ", i, ", t2 = ", t2)) }

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix

    for (j in 1 : N_bar) {
      if (j %% 1000 == 0) { print(paste0("i = ", i, " & j = ", j)) }
      p <- k[[j]]
      names(p) <- doc[[j]]
      if ((t1 %in% doc[[j]]) == TRUE) {
        x1 <- as.numeric(p[t1])
      } else {
        x1 <- 0
      }
      if ((t2 %in% doc[[j]]) == TRUE) {
        x2 <- as.numeric(p[t2])
      } else {
        x2 <- 0
      }
      x3 <- sum(p) - x1 - x2
      x <- c(x1, x2, x3)

      mvhgt_score_submat[i, j] <- pmvhyper(x = x, n = n, k = sum(p), log = TRUE)

      if (mvhgt_score_submat[i, j] == Inf) {
        mvhgt_score_submat[i, j] <- pmvhyper3(x = x, n = n, k = sum(p), log = TRUE)
        count <- count + 1
        print(paste0("(i,j) = (", i, ",", j, "), c = ", count, ", -log(P) = ", mvhgt_score_submat[i, j], ", x = (", x1, ",", x2, ",", x3, ")", ", n = (", n1, ",", n2, ",", n3, ")"))
      }
    }

    mvhgt_ranked_docs <- sort(rank(-mvhgt_score_submat[i, ], ties.method = "random"), index.return = TRUE)$ix

    p10_mat[i, "MVHGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "MVHGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "MVHGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }

  save(mvhgt_score_submat, file = paste0(output_dir, "mvhgt_score_submat_new.Rdata"))
  save(p10_mat, file = paste0(output_dir, "p10_mat_new.Rdata"))
}

## number of infinite scores
length(which(mvhgt_score_submat == Inf))
# > length(which(mvhgt_score_submat == Inf))
# [1] 0

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#   MVHGT-RND    MVHGT-TP MVHGT-TPIDF    TP-TPIDF 
#  0.01152406  0.18000000  2.29000000  0.80000000
apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#    MVHGT-RND     MVHGT-TP  MVHGT-TPIDF     TP-TPIDF 
# 3.745847e-05 9.361430e-01 2.602621e+00 1.814643e+00



## calculate P@10 scores for all queries q = ("imf", t)
count <- 0
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("MVHGT-RND", "MVHGT-TP", "MVHGT-TPIDF", "TP-TPIDF")
t1 <- "imf"
n1 <- as.numeric(K[t1])
mvhgt_score_submat <- mat.or.vec(nr = cut, nc = N_bar)

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]
  n2 <- as.numeric(K[t2])
  n3 <- N - n1 - n2
  n <- c(n1, n2, n3)
  if (i %% 1 == 0) { print(paste0("i = ", i, ", t2 = ", t2)) }

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix

    for (j in 1 : N_bar) {
      if (j %% 1000 == 0) { print(paste0("i = ", i, " & j = ", j)) }
      p <- k[[j]]
      names(p) <- doc[[j]]
      if ((t1 %in% doc[[j]]) == TRUE) {
        x1 <- as.numeric(p[t1])
      } else {
        x1 <- 0
      }
      if ((t2 %in% doc[[j]]) == TRUE) {
        x2 <- as.numeric(p[t2])
      } else {
        x2 <- 0
      }
      x3 <- sum(p) - x1 - x2
      x <- c(x1, x2, x3)

      mvhgt_score_submat[i, j] <- pmvhyper(x = x, n = n, k = sum(p), log = TRUE)

      if (mvhgt_score_submat[i, j] == Inf) {
        mvhgt_score_submat[i, j] <- pmvhyper3(x = x, n = n, k = sum(p), log = TRUE)
        count <- count + 1
        print(paste0("(i,j) = (", i, ",", j, "), c = ", count, ", -log(P) = ", mvhgt_score_submat[i, j], ", x = (", x1, ",", x2, ",", x3, ")", ", n = (", n1, ",", n2, ",", n3, ")"))
      }
    }

    mvhgt_ranked_docs <- sort(rank(-mvhgt_score_submat[i, ], ties.method = "random"), index.return = TRUE)$ix

    p10_mat[i, "MVHGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "MVHGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "MVHGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }

  save(mvhgt_score_submat, file = paste0(output_dir, "mvhgt_score_submat_imf.Rdata"))
  save(p10_mat, file = paste0(output_dir, "p10_mat_imf.Rdata"))
}

## number of infinite scores
length(which(mvhgt_score_submat == Inf))
# > length(which(mvhgt_score_submat == Inf))
# [1] 0

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#   MVHGT-RND    MVHGT-TP MVHGT-TPIDF    TP-TPIDF 
#  0.01233398  0.13000000  2.20000000  1.18000000
apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#    MVHGT-RND     MVHGT-TP  MVHGT-TPIDF     TP-TPIDF 
# 5.185109e-05 8.487067e-01 2.666667e+00 2.293513e+00


## calculate P@10 scores for all queries q = ("comment", t)
count <- 0
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("MVHGT-RND", "MVHGT-TP", "MVHGT-TPIDF", "TP-TPIDF")
t1 <- "comment"
n1 <- as.numeric(K[t1])
mvhgt_score_submat <- mat.or.vec(nr = cut, nc = N_bar)

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]
  n2 <- as.numeric(K[t2])
  n3 <- N - n1 - n2
  n <- c(n1, n2, n3)
  if (i %% 1 == 0) { print(paste0("i = ", i, ", t2 = ", t2)) }

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix

    for (j in 1 : N_bar) {
      if (j %% 1000 == 0) { print(paste0("i = ", i, " & j = ", j)) }
      p <- k[[j]]
      names(p) <- doc[[j]]
      if ((t1 %in% doc[[j]]) == TRUE) {
        x1 <- as.numeric(p[t1])
      } else {
        x1 <- 0
      }
      if ((t2 %in% doc[[j]]) == TRUE) {
        x2 <- as.numeric(p[t2])
      } else {
        x2 <- 0
      }
      x3 <- sum(p) - x1 - x2
      x <- c(x1, x2, x3)

      mvhgt_score_submat[i, j] <- pmvhyper(x = x, n = n, k = sum(p), log = TRUE)

      if (mvhgt_score_submat[i, j] == Inf) {
        mvhgt_score_submat[i, j] <- pmvhyper3(x = x, n = n, k = sum(p), log = TRUE)
        count <- count + 1
        print(paste0("(i,j) = (", i, ",", j, "), c = ", count, ", -log(P) = ", mvhgt_score_submat[i, j], ", x = (", x1, ",", x2, ",", x3, ")", ", n = (", n1, ",", n2, ",", n3, ")"))
      }
    }

    mvhgt_ranked_docs <- sort(rank(-mvhgt_score_submat[i, ], ties.method = "random"), index.return = TRUE)$ix

    p10_mat[i, "MVHGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "MVHGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "MVHGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }

  save(mvhgt_score_submat, file = paste0(output_dir, "mvhgt_score_submat_comment.Rdata"))
  save(p10_mat, file = paste0(output_dir, "p10_mat_comment.Rdata"))
}

## number of infinite scores
length(which(mvhgt_score_submat == Inf))
# > length(which(mvhgt_score_submat == Inf))
# [1] 0

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#   MVHGT-RND    MVHGT-TP MVHGT-TPIDF    TP-TPIDF 
#  0.02419001  1.70000000  0.89000000  4.80000000
apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#    MVHGT-RND     MVHGT-TP  MVHGT-TPIDF     TP-TPIDF 
# 0.0005321128 1.0396191999 1.5820648827 2.6514718611


## calculate P@10 scores for all queries q = ("lagarde", t)
count <- 0
p10_mat <-  mat.or.vec(nr = cut, nc = 4)
colnames(p10_mat) <- c("MVHGT-RND", "MVHGT-TP", "MVHGT-TPIDF", "TP-TPIDF")
t1 <- "lagarde"
n1 <- as.numeric(K[t1])
mvhgt_score_submat <- mat.or.vec(nr = cut, nc = N_bar)

for (i in 1 : cut) {
  t2 <- top_bursty_term[i]
  n2 <- as.numeric(K[t2])
  n3 <- N - n1 - n2
  n <- c(n1, n2, n3)
  if (i %% 1 == 0) { print(paste0("i = ", i, ", t2 = ", t2)) }

  if(t2 %in% top_six_bursty_terms == FALSE) {
    query <- c(t1, t2)
    tp_ranked_docs <- sort(rank(-apply(tp_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix
    tpidf_ranked_docs <- sort(rank(-apply(tpidf_score_submat[query,], 2, sum), ties.method = "random"), index.return = TRUE)$ix

    for (j in 1 : N_bar) {
      if (j %% 1000 == 0) { print(paste0("i = ", i, " & j = ", j)) }
      p <- k[[j]]
      names(p) <- doc[[j]]
      if ((t1 %in% doc[[j]]) == TRUE) {
        x1 <- as.numeric(p[t1])
      } else {
        x1 <- 0
      }
      if ((t2 %in% doc[[j]]) == TRUE) {
        x2 <- as.numeric(p[t2])
      } else {
        x2 <- 0
      }
      x3 <- sum(p) - x1 - x2
      x <- c(x1, x2, x3)

      mvhgt_score_submat[i, j] <- pmvhyper(x = x, n = n, k = sum(p), log = TRUE)

      if (mvhgt_score_submat[i, j] == Inf) {
        mvhgt_score_submat[i, j] <- pmvhyper3(x = x, n = n, k = sum(p), log = TRUE)
        count <- count + 1
        print(paste0("(i,j) = (", i, ",", j, "), c = ", count, ", -log(P) = ", mvhgt_score_submat[i, j], ", x = (", x1, ",", x2, ",", x3, ")", ", n = (", n1, ",", n2, ",", n3, ")"))
      }
    }

    mvhgt_ranked_docs <- sort(rank(-mvhgt_score_submat[i, ], ties.method = "random"), index.return = TRUE)$ix

    p10_mat[i, "MVHGT-RND"] <- top * top / length(union(term_to_doc_ids_lst[[t1]], term_to_doc_ids_lst[[t2]]))
    p10_mat[i, "MVHGT-TP"] <- length(intersect(tp_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "MVHGT-TPIDF"] <- length(intersect(tpidf_ranked_docs[1 : top], mvhgt_ranked_docs[1 : top]))
    p10_mat[i, "TP-TPIDF"] <- length(intersect(tp_ranked_docs[1 : top], tpidf_ranked_docs[1 : top]))
  }

  save(mvhgt_score_submat, file = paste0(output_dir, "mvhgt_score_submat_lagarde.Rdata"))
  save(p10_mat, file = paste0(output_dir, "p10_mat_lagarde.Rdata"))
}

## empty rows ids
empty_row_id <- numeric(length(top_six_bursty_terms))
for (i in 1 : length(top_six_bursty_terms)) {
  empty_row_id[i] <- which(top_bursty_term == top_six_bursty_terms[i])
}

## number of infinite scores
length(which(mvhgt_score_submat == Inf))
# > length(which(mvhgt_score_submat == Inf))
# [1] 0

## summary stats
apply(p10_mat[-empty_row_id, ], 2, mean)
# > apply(p10_mat[-empty_row_id, ], 2, mean)
#   MVHGT-RND    MVHGT-TP MVHGT-TPIDF    TP-TPIDF 
#  0.04131877  0.27000000  1.21000000  4.58000000

apply(p10_mat[-empty_row_id, ], 2, sd)
# > apply(p10_mat[-empty_row_id, ], 2, sd)
#   MVHGT-RND    MVHGT-TP MVHGT-TPIDF    TP-TPIDF 
# 0.001907695 0.972916058 1.965844718 3.318695131

########################################################################################################################
#
# Description: This R script performs document proportion versus total term proportion regression analyses.
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
rm(list = ls()) # deletes all variables initialized in R session
library(magicaxis)
input_dir <- paste0(working_dir, "/stats/") # deletes any already initialized variables in R session
plots_output_dir <- paste0(working_dir, "/plots/")
setwd(working_dir)
set.seed(83725) # initialize random seed


## set plot colours 
gg_color_hue = function(n, alpha = 1) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1 : n]
}
cols <- c("grey25", gg_color_hue(n = 3))
gray25 <- cols[1]; red <- cols[2]; green <- cols[3]; blue <- cols[4]


## load some NYSK data summary statistics
load(paste0(input_dir, "N-uppercase.Rdata"))
load(paste0(input_dir, "K-uppercase.Rdata"))
load(paste0(input_dir, "N_bar-uppercase.Rdata"))
load(paste0(input_dir, "K_bar-uppercase.Rdata"))

x <- as.numeric(K_bar / N_bar)
y <- as.numeric(K / N)

## calculate correlation between K_bar / N_bar and K / N
corr <- format(round(cor(x, y), 2), nsmall = 2)
# > corr
# [1] "0.89"

## R squared
rsq <- format(round(cor(x, y)^2, 2), nsmall = 2)
# > rsq
# [1] "0.79"

## fit linear regression model
fit <- lm(y ~ x)
m <- fit$coeff[1]
b <- fit$coeff[2]

## sample points for plotting purposes
s <- sample(1 : length(x), size = 500, prob = 1 / (1 + (x / (1 - x))^(-3)))
x_sample <- x[s]
y_sample <- y[s]

pdf(paste0(plots_output_dir, "figure-2a.pdf"), width = 6, height = 6, useDingbats = FALSE)
plot(x = x,
  y = y,
  main = bquote("NYSK Dataset" ~ K / N ~ "vs." ~ bold(K) / bold(N) ~ "Plot"),
  xlab = expression(K / N),
  ylab = expression(bold(K) / bold(N)),
  axes = FALSE,
  type = "n")
magaxis()
grid(lwd = 2)
points(x = x_sample, y = y_sample, pch = 16, cex = 0.75)
abline(m, b, col = red, lwd = 3)
legend( 
  "topleft", 
  paste0("R-sq: ", rsq), 
  cex = 1,
  inset = 0.03,
  box.lwd = 0,
  box.col = "grey90",
  bg = "grey90")
dev.off()

## fit linear regression model
fit2 <- lm(log10(1 / y) ~ log10(1 / x))
m2 <- fit2$coeff[1]
b2 <- fit2$coeff[2]

rsq2 <- format(round(cor(log10(1 / x), log10(1 / y))^2, 2), nsmall = 2)

s <- sample(1 : length(x), size = 2000)
x_sample <- x[s]
y_sample <- y[s]

pdf(paste0(plots_output_dir, "figure-2b.pdf"), width = 6, height = 6, useDingbats = FALSE)
plot(x = 1 / x,
  y = 1 / y,
  main = bquote("NYSK Dataset" ~ -log(K / N) ~ "vs." ~ -log(bold(K) / bold(N)) ~ "Plot"),
  xlab = expression(-log(K / N)),
  ylab = expression(-log(bold(K) / bold(N))),
  axes = FALSE,
  log = "xy",
  type = "n")
magaxis()
grid(lwd = 2)
points(x = 1 / x_sample, y = 1 / y_sample, pch = 16, cex = 0.75)
abline(m2, b2, col = red, lwd = 3)
legend( 
  "topleft", 
  paste0("R-sq: ", rsq2), 
  cex = 1,
  inset = 0.03,
  box.lwd = 0,
  box.col = "grey90",
  bg = "grey90")
dev.off()

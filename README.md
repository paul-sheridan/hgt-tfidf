# A Hypergeometric Test Interpretation of a Common tf-idf Variant Computer Code
This repository contains computer code for reproducing the results described in the manuscript “A hypergeometric test interpretation of a common tf-idf variant”. ArXiv preprint link: https://arxiv.org/abs/2002.11844. Research Gate link: https://www.researchgate.net/publication/339550104_A_hypergeometric_test_interpretation_of_a_common_tf-idf_variant.

## Initial Setup

* Clone this repository by running the command `git clone https://github.com/paul-sheridan/hgt-tfidf.git`.
* Download the file `nysk.xml` from the UCI Machine Learning Repository at https://archive.ics.uci.edu/ml/datasets/NYSK and place it inside the `hgt-tfidf` folder.

## NYSK Dataset Preprocessing and Summary Stats Calculation

* To prepare the NYSK data for analysis, open the `nysk-preprocessing.ipynb` file as a Jypyter Notebook and run it. This script will save  a preprocessed version of the NYSK data in JSON format to the file `nysk-processed.json`.
* Run the R script `nysk-stats-calc.R` to calculate some basic statistics that will be used in the experiments.

## Reproducing our Results

* Run the R script `experiment-1.R` to reproduce the results of comparing tp-idf with the hypergeometric test in the one-term query document retrieval scenario.
* Run the R script `experiment-2.R` to reproduce the results of comparing tp-idf with the hypergeometric test in the document summarization scenario.
* Run the R script `experiment-3.R` to reproduce the results of comparing tp-idf with the hypergeometric test in the two-term query document retrieval scenario.
* Run the R script `experiment-4.R` to reproduce the results of comparing tp-idf with the multivariate hypergeometric test in the two-term query document retrieval scenario.
* Run the R script `experiment-5.R` to reproduce the regression analyses for the negated logarithmically scaled document proportion versus the negated logarithmically scaled total term proportion and the document proportion versus the total term proportion, respectively.

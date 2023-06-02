# The Hypergeometric Test Performs Comparably to TF-IDF on Standard Text Analysis Tasks Computer Code
This repository contains computer code for reproducing the results described in the manuscript “The hypergeometric test performs comparably to TF-IDF on standard text analysis tasks”. ArXiv preprint link: https://arxiv.org/abs/2002.11844.

## Initial Setup

Clone this repository by running the command

```
git clone https://github.com/paul-sheridan/hgt-tfidf.git
```

## Reproducing our NYSK Dataset Case Study Results

* Download the file `nysk.xml` from the UCI Machine Learning Repository ([download page](https://archive.ics.uci.edu/ml/datasets/NYSK)) and place it inside the `hgt-tfidf/nysk/dataset` folder.
* To prepare the NYSK data for analysis, open the `hgt-tfidf/nysk/dataset/nysk-preprocessing.ipynb` file as a Jypyter Notebook and run it. This script will save a preprocessed version of the NYSK data in JSON format to the file `nysk-processed.json`.
* Create the folders `hgt-tfidf/nysk/stats` and `hgt-tfidf/nysk/plots`.
* Run the `hgt-tfidf/nysk/nysk-generate-stats.Rmd` R Markdown file to calculate some basic statistics that will be used in the experiments.
* To reproduce the one-term query document retrieval case results of Figure 1a, run the `hgt-tfidf/nysk/nysk-doc-retrieval-one-term-query-experiments.Rmd` R Markdown file.
* To reproduce the two-term query document retrieval case of Table 1, run the `hgt-tfidf/nysk/nysk-doc-retrieval-two-term-query-experiments.Rmd` R Markdown file.
* To reproduce the document summarization case result of Figure 1b, run the `hgt-tfidf/nysk/nysk-doc-summary-experiments.Rmd` R Markdown file.

## Reproducing our Cranfield 1400 Test Collection Case Study Results

* Download the files `cran.all.1400.xml`, `cran.qry.xml`, and `cranqrel.trec.txt` from the [cranfield-trec-dataset](https://github.com/oussbenk/cranfield-trec-dataset) GitHub repo and place them inside the `hgt-tfidf/cranfield/dataset` folder. We used the files associated with commit Id `1208e6edfb6cb2527b2c44398d3d8fefd3249144` in our experiments.
* To prepare the Cranfield data for analysis, open the `hgt-tfidf/cranfield/dataset/cran-preprocessing.ipynb` file as a Jypyter Notebook and run it. This script will save a preprocessed version of the Cranfield documents in JSON format to the file `cran-docs-preprocessed.json`, and a preprocessed version of the Cranfield documents in JSON format to the file `cran-queries-preprocessed.json`.
* To reproduce the preformance evaluation results of Table 2, run the `hgt-tfidf/cranfield/cranfield-experiments.Rmd` R Markdown file.

## Reproducing our Twenty Newsgroups Dataset Case Study Results

* To reproduce the preformance evaluation results of Table 3, run the `hgt-tfidf/twenty-newsgroups/dataset/twenty-newsgroups-analysis.ipynb` Jypyter Notebook file.

# cofa

This R package provides an implementation of the methods described in:

[Categorical Co-Frequency Analysis: Clustering Diagnosis Codes to Predict Hospital Readmissions](https://arxiv.org/abs/1909.00306)

Hallee E. Wong, Brianna C. Heggeseth, Steven J. Miller

# Installation

To install and load this package in R from GitHub, run the following commands:
```
install.packages("devtools")
library(devtools) 
install_github("halleewong/cofa")
library(cofa)
```

`test/test_adult.R` shows an example script using functions from this package on the [adult data set from the UCI machine learning repository](https://archive.ics.uci.edu/ml/datasets/Adult). 
 
# Development

This package is under active development and may change substantially with each commit. 

# Reference

If you use our code, please cite our [paper](https://arxiv.org/abs/1909.00306).

```
@misc{wong2019categorical,
    title={Categorical Co-Frequency Analysis: Clustering Diagnosis Codes to Predict Hospital Readmissions},
    author={Hallee E. Wong and Brianna C. Heggeseth and Steven J. Miller},
    year={2019},
    eprint={1909.00306},
    archivePrefix={arXiv},
    primaryClass={stat.AP}
}
```

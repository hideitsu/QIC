# QIC
Quantile-based Information Theoretic Clustering algorithm


Originally coded at April, 17, 2014


REFERENCES

This code implements algorithms proposed in the following papers by R.  
- "A Non-parametric Clustering Algorithm With A Quantile-based Likelihood Estimator",
   Hideitsu Hino, Noboru Murata
   Neural Computation

- "A Non-Parametric Maximum Entropy Clustering",
   Hideitsu Hino, Noboru Murata
   International Conference on Artificial Neural Networks 2014 (ICANN2014)

When you use this source code in your study and write papers, please cite above two papers.



## You must install following libraries in advance:
"MASS", "clues", and "pracma"


## a test code using the data "wine" is written in "demo.r"

> source("demo.r")

The main clustering routine is coded in "QIC.r", which provides a function QIC.

If you want to perform clustering with gradually decreasing alpha (the method proposed in ICANN2014 paper), please set the parameter 

> anneal <- TRUE

in the 10th line of the file "demo.r"

If you want to optimize alpha using the cluster conditional entropy (the method proposed in Neural Computation paper), please set the parameter

> anneal <- FALSE



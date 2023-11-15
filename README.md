# Convexity - Data driven model discovery in F#

This library aims to provide a series of novel algorithms for data driven model discovery.

In particular sparsity promoting algorithms are becoming a powerful method to find good models which
also generalise well and prevent overfitting.

The Dynamic Mode Decomposition (Dmd) is one such algorithm. It can be used to find dominant system dynamics from
high dimensional time series data. When combined with the singular value decomposition and appropriate thresholding the algorithm
provides a powerfull method for finding the low rank structures in data, even in the presence of some noise.

A vast amount of information can be found online and a good introductory reference I have found is [Data Driven Science and Engineering](https://databookuw.com/).

Another novel algorithm is the Sparse Identification of Non Linear Dynamics. This algorithm takes the Dmd and extends it by instead solving
the regression over a projected nonlinear manifold of candidate terms which are selected by the user. The hope is that in this new space only a few dominant terms will prevail when solved using a sparsity promoting technique such as Least Squares with L1 norm regularisation.

Other algorithms for computing the sparse solution can be employed giving various benifits and drawbacks from each.
Currently the following sparse solvers exist.


| Solver Name | Description |
| ----------- | ------------:
|    stlsq    | The sequential thresholded least squares |
|    lasso    | The least absolute shrinkage and selection operator |

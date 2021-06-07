# Goal of this simulative study

Johansen test is adequate cointegration test when there are more than two tested series at the same time. This test is performed to estimate the number of cointegrated vectors (r) in the system. If there is any cointegration in the model then 0 < r < k, where k is the number of tested time-series. To perform the simulation I generate 1000 trajectories of three different data generating process (DGP): (1) bivariate cointegrated system, (2) treevariate system with only one cointegrated vector, (3) treevariate system with two cointegrated vector.

[DGPs with formula]

In this short analyses I focus on the effect of the different parameters of the DGPs: length of the simulated series (n), auto-correlation parameter (p) and the variance of the innovations (s).

The system decomposition is not unique, so we can only estimate the cointegration rank r [Kirchg¨assner and Wolters, 2007]. The method can be performed with several tests, in this paper I chose the trace test. It gives a vector of the test statistics as a result and that may be compared to critical values. The null hypothesis is that r ≤ x, where x = 0, 1, 2, ..., k − 1. The number of cointegrated vectors is the smallest x, under which the null hypothesis is not rejected at a given significance level. 

The methodology of estimation means that one has to perform tree different statistical test to identify the cointegration rank in a system, and so can commit multiple type I error (reject a true null hyphothesis, over estimation of the rank) or type II error (accept a false null hyphothesis, underestimation of the rank).


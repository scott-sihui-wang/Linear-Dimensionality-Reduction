# Linear Dimensionality Reduction

## 1. Introduction

**Topics:** _Statistical Machine Learning_, _Information Theory_, _Optimization_

**Skills:** _Matrix Calculus_, _Optimization on Manifolds_, _Information Geometry_, _Lie Algebra_, _Matlab_, _Manopt (library)_

This project corresponds to my [paper](https://ieeexplore.ieee.org/abstract/document/9518004), `Linear Discriminant Analysis under f-divergence Measures` on `IEEE International Symposium on Information Theory (ISIT) 2021`, and the [journal](https://www.mdpi.com/1099-4300/24/2/188), `Discriminant Analysis under f-Divergence Measures` on `MDPI Entropy`.

## 2. Contributions

In this project, my contributions are:

- I derived and mathematically proved the formulas for matrix calculus that is the foundation for gradient computation; [report](/reports/Maximizing_Divergence_after_Linear_Dimensionality_Reduction.pdf)

- I proposed a first-order `gradient descent algorithm` based on `Lie algebra` and its `exponential maps`; [report](/reports/Statistical_LDAs_Results_on_Various_Datasets.pdf)

- I formulated `linear dimensionality reduction` as the optimization problem on `Stiefel manifolds` and `Grassmann manifolds`; I tried first order methods (`deepest descent method`) and second order methods (`trust region method`) to solve the optimization problems; I conducted experiments on various datasets (e.g., `MNIST`, `IRIS`) and reported the results;

- I performed an `information geometric analysis` for `linear dimensionality reduction` and proposed an dimensionality reduction algorithm to maximally preserve `information geometric distance`; [report](/reports/An_Information_Geometric_Interpretation_of_Linear_Dimensionality_Reduction.pdf)

- I performed a `variational analysis of Total Variation` and proposed a gradient descent algorithm. [report](/reports/Gradient_Based_Optimization_on_Total_Variation.pdf)

## 3. Demo

Here is the result of `linear dimensionality reduction` on `MNIST` dataset. The `784` dimensional data is linearly projected to a `3` dimensional subspace with their distinctive features preserved (in an `information theoretic` view) as much as possible.

![demo](/demo/MNIST_dim3.png)

# Tight Concentrations and Confidence Sequences from the Regret of Universal Portfolio

## Overview

This repository contains the software implementation of the methods presented in [Tight Concentrations and Confidence Sequences from the Regret of Universal Portfolio](https://ieeexplore.ieee.org/document/10315047) (also available on [Arxiv](https://arxiv.org/pdf/2110.14099.pdf)). We introduce novel methods for deriving time-uniform concentration inequalities and confidence sequences using portfolio algorithms. These results have direct applications in statistics and machine learning, specifically for constructing valid and tight confidence intervals.

## Problem Statement

Estimating the mean of bounded random variables from samples is a fundamental task in statistics. Confidence sequences---confidence intervals that hold uniformly over time---are particularly valuable for adaptive decision making. However, traditional methods often produce vacuous intervals for small sample sizes, limiting their practical usability.

## Our Approach

Building upon the theory of online betting and portfolio selection, we propose algorithms that leverage the regret of universal portfolio algorithms to generate tighter and more practical confidence sequences. In particular, the obtained confidence sequences are state-of-the-art and never vacuous, even with a single sample.

## Features

- **Time-uniform Confidence Sequences**: Generate confidence intervals that hold across all time steps.
- **Never-vacuous Guarantees**: Ensures practical usability even for small sample sizes.
- **Optimal Asymptotics**: Matches the best possible performance for large datasets.
- **Efficient Algorithms**: Includes both exact and approximated methods for fast computations.

## Algorithms
   
- PRECiSE-CO96: A method leveraging regret guarantees of universal portfolio algorithms for optimal confidence sequences.
- PRECiSE-A-CO96: A computationally efficient variant with constant-time updates per sample.
- PRECiSE-R70: A portfolio algorithm inspired by Robbins' mixture method, achieving state-of-the-art asymptotic performance while satisfying the law of iterated logarithm.

The matlab folder contains the implementation of all the algorithms in Matlab, while the Python folder contains the implementation of PRECiSE-CO96 only.
Take a look at the demo scripts in both directories for an example of usage.

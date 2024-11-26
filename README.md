# PRECiSE

This software package implements the algorithms in the paper

Orabona and Jun,
"Tight Concentrations and Confidence Sequences from the Regret of Universal Portfolio",
[IEEE Trans. on Information Theory](https://ieeexplore.ieee.org/document/10315047) (also available on [Arxiv](https://arxiv.org/pdf/2110.14099.pdf))

The algorithms construct a sequence of time-uniform confidence intervals, also known as confidence sequences, using a betting approach.
In particular, the regret of universal portfolio algorithms is used to construct implicit concentration inequalities that are numerically inverted.

The matlab folder contains the implementation of all the algorithms in Matlab, while the pyhton folder contains the implementation of PRECiSE-CO96 only.
Take a look at the demo scripts in both directories for an example of usage.

# Multicomponent T2 estimation code

This repository provides the source code and examples of a robust algorithm for the estimation of multi-component T2 distributions from MRI or NMR relaxometry data. The algorithm uses the extended phase graph (EPG) algorithm to jointly estimate the flip angle and T2 distribution components. This compensates for effects such as stimulated echoes.

One application for the algorithm is the estimation of myelin water content in the human brain. A Bayesian formulation is used to incorporate prior information about the location of the component peaks. The resulting T2 estimates contain fewer outliers than standard algorithms such as least-squares. The code is written is MATLAB version 2011b.

## Key functions

* `bayesian_estimate.m` – Main function to estimate components of T2 distribution.
* `demo_bayesian.m` – Demonstration of the Bayesian estimation algorithm for multicomponent T2 estimation. See Example 1 below.
* `demo_epg_derivatives.m` – Demonstration of the EPG signal derivative. See [Example 2] below.
* `epg.m` – MATLAB implementation of the Extended Phase Graph (EPG) algorithm.
* `epg_mex.cpp` – MEX file implementation of the Extended Phase Graph (EPG) algorithm.
* `epg_derivatives.m` – MATLAB implementation of the EPG signal derivatives.
* `epg_derivatives_mex.cpp` – MEX implementation of the EPG signal derivatives
* `test_implementations.m` – Script to test the different function implementations.

## Example 1: T2 Distribution Estimation

The plot below demonstrates an example of the Bayesian estimation algorithm. Synthetic data was generated using the continuous distribution. The corresponding discrete components were estimated using the proposed Bayesian algorithm. For low SNR this algorithm is more reliable than a simple least squares fitting.

![T2 components](https://github.com/kelvinlayton/T2estimation/raw/master/images/example_bayesian.png "Estimation performance")


## Example 2: Derivative of EPG signal

A key part of the Bayesian estimation is to approximate likelihood function using the derivatives of the EPG signal. The derivates are calculated using a recursive algorithm that is computationally efficient. The derivative could also be used to speed up gradient-based optimisation algorithms. The following plot shows a simple demonstration of a local linear approximation using the derivative algorithm.

![EPG derivative](https://github.com/kelvinlayton/T2estimation/raw/master/images/example_epg_derivatives.png "Derivative of EPG algorithm")

## Reference

The detailed algorithm and additional results are presented in the following publication:

Layton, K. J., Morelande, M., Wright, D., Farrell, P. M., Moran, B., & Johnston, L. A. (2013). Modelling and Estimation of Multicomponent T2 Distributions. *IEEE Transactions on Medical Imaging*, 32(8), 1423–1434.
DOI: [10.1109/TMI.2013.2257830](http://dx.doi.org/10.1109/TMI.2013.2257830)

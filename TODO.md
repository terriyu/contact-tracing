Codetags used: TODO, FIXME, HACK, IDEA, QUESTION, NOTE

To do list
==========

+ Implement other methods for generating the initial infection.  Currently, a
  simple homogeneous Bernoulli process is used.
+ Use a better directory structure (source code in its own directory?)
+ Put parameters for run_simulations.R script in separate file?
+ (?) Make software into an R package

Where I left off the project:

+ Look into R package `coda` https://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/
+ Write tests for fit_ct_sample.R
+ Try to separate out the behavior of the three stochastic processes: Y, Z0, and W
+ Create some good initial networks for testing MCMC (can try simulating some samples and picking one)

Contact tracing
===============

Repository for development of R code, which implements the contact tracing
method in social networks.

See Krista Gile's thesis in references below for more details about contact
tracing.

Features
--------
* Simulate a network and spread an infection across it
* Perform contact tracing sampling on a simulated infected network
* [not yet implemented] Perform Markov Chain Monte Carlo Maximum Likelihood
  Estimation (MCMCMLE) for a Exponential Random Graph Model (ERGM), given a
  contact tracing sample

Dependencies
------------

The `network` package is needed to run this code. It is part of a larger
package called `statnet`.  If you want to run the tests, the `testthat` package
must be installed.

How to run the tests
--------------------
The testthat unit testing framework is used in this project.  The `testthat`
package can be downloaded and installed via the usual methods from CRAN.

Individual tests are in the `tests/testthat` subdirectory, while scripts to run
multiple tests at a time are in the tests directory.

To run all tests, there are three ways to do it:

* (Command line) run the following command in the root directory:
`Rscript tests/run_all_tests.R`

* (In a R session) run the following command in the root directory:
`source('tests/run_all_tests.R')`

* (Bash script) run the following Bash script in the root directory:
`./run_all_tests.sh`

The bash script should work on Mac OS X, since bash is the default shell.

To run an individual test:

1. Load the `testthat` library: `library(testthat)`
2. Run the individual test, for example: `test_file('test_simulate_network.R')`
   Use whatever is the correct path to the test you want to run.

Code style check
----------------
If you are using the R package lintr to check code style, there is already a
lintr config file `.lintr` which avoids running a couple linters.

References
----------

*Contact tracing*
+ Krista Gile PhD thesis, Inference from Partially-Observed Network Data,
  University of Washington (2008) -- see Chapter 4

*statnet package (which contains the network package)*
+ [Official statnet site](http://www.statnet.org/)

*testthat unit testing framework*

+ [Hadley Wickham Github repository](https://github.com/hadley/testthat)
+ [Hadley Wickham R journal article](https://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf)
+ [Hadley Wickham test workflow](http://r-pkgs.had.co.nz/tests.html)
+ [Stefan Feuerriegel slides](http://www.is.uni-freiburg.de/ressourcen/algorithm-design-and-software-engineering-oeffentlicher-zugriff/11_softwaretesting.pdf)
+ [John Cook example](http://www.johndcook.com/blog/2013/06/12/example-of-unit-testing-r-code-with-testthat/)

*R style guides*

+ [Google R style guide](https://google.github.io/styleguide/Rguide.xml)
+ [Hadley Wickham R style guide](http://adv-r.had.co.nz/Style.html)

*Markdown syntax*

+ [Daring Fireball](https://daringfireball.net/projects/markdown/syntax)

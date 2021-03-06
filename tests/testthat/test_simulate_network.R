context("Test simulate_network.R")

source('../../simulate_network.R')

test_that("Test create.edge.prob.mtx", {

 # Test for case where there is only one node class
  expect_equal(create.edge.prob.mtx(rep(1, times = 3), 0.5), matrix(c(0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0), nrow = 3, ncol = 3))

  # Test for case where there is more than one node class
  # Number of nodes in each class
  class.labels <- rep(1:3, times = 1:3)
  # Matrix containing probability of link between class i and class j
  P.ij <- matrix(c(0.5, 0.4, 0, 0.4, 0.8, 0.4, 0, 0.4, 0.5), nrow = 3, ncol = 3)
  # Expected matrix
  edge.prob.mtx <- matrix(c(0, 0.4, 0.4, 0, 0, 0, 0.4, 0, 0.8, 0.4, 0.4, 0.4, 0.4, 0.8, 0, 0.4, 0.4, 0.4, 0, 0.4, 0.4, 0, 0.5, 0.5, 0, 0.4, 0.4, 0.5, 0, 0.5, 0, 0.4, 0.4, 0.5, 0.5, 0), nrow = 6, ncol = 6)
  expect_equal(create.edge.prob.mtx(class.labels, P.ij), edge.prob.mtx)
})

test_that("Test generate.network", {
  set.seed(137)

  class.labels <- rep(1:3, times = c(25,5,25))
  # Test non-square P.ij matrix
  expect_error(generate.network(class.labels, matrix(runif(10), nrow=2, ncol=5)), regexp = "P.ij not a square matrix")
  # Test non-symmetric P.ij matrix
  expect_error(generate.network(class.labels, matrix(1:9, nrow=3, ncol=3)), regexp = "P.ij is not symmetric")
  # Test class.labels and P.ij are consistent
  expect_error(generate.network(class.labels, matrix(rep(1,4), nrow=2, ncol=2)), regexp = "class.labels and P.ij have inconsistent dimensions")

  # Number of nodes in each class
  class.labels <- rep(1:3, times = 1:3)
  num.nodes <- length(class.labels)
  # Matrix containing probability of link between class i and class j
  P.ij <- matrix(c(0.5, 0.4, 0, 0.4, 0.8, 0.4, 0, 0.4, 0.5), nrow = 3, ncol = 3)

  # Test that class.labels and class.names are consistent
  expect_error(generate.network(class.labels, P.ij, c("a", "b")), regexp = "Number of class names in class.names inconsistent with class.labels")

  # Test that class attribute set correctly
  net <- generate.network(rep(1, times = num.nodes), 0.2)
  expect_equal(get.vertex.attribute(net, "class.label"), c(1, 1, 1, 1, 1, 1))

  net <- generate.network(class.labels, P.ij)
  expect_equal(get.vertex.attribute(net, "class.label"), c(1, 2, 2, 3, 3, 3))

  class.names <- rep(c("a", "b", "c"), times = 1:3)
  net <- generate.network(class.labels, P.ij, class.names)
  expect_equal(get.vertex.attribute(net, "class.name"), c("a", "b", "b", "c", "c", "c"))

  socio.mtx <- as.sociomatrix(net)

  # Test that all values of sociomatrix are 0 or 1
  expect_true(all((socio.mtx == 0) | (socio.mtx == 1)))

  # Test that sociomatrix is symmetric
  expect_true(all(socio.mtx == t(socio.mtx)))
})

test_that("Test spread.infection", {
  set.seed(137)

  net <- network(matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 6, ncol = 6), directed = FALSE)
  W <- matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 6, ncol = 6)
  Z0 <- c(1,0,0,0,0,0)
  Z.expect <- c(1, 1, 1, 1, 0, 0) # Expected value for Z after infection has been spread over network

  # Test error checking for valid Bernoulli parameter eta
  expect_error(spread.infection(net, -9, NULL, "nodes_fixed"), regexp = "Bernoulli process parameter eta must be between 0 and 1")
  expect_error(spread.infection(net, 81, NULL, "nodes_fixed"), regexp = "Bernoulli process parameter eta must be between 0 and 1")
  # Test error checking for valid Bernoulli parameter tau
  expect_error(spread.infection(net, NULL, -2, "edges"), regexp = "Bernoulli process parameter tau must be between 0 and 1")
  expect_error(spread.infection(net, NULL, 100, "edges"), regexp = "Bernoulli process parameter tau must be between 0 and 1")

  # Test error checking for initial.infect.method variable
  expect_error(spread.infection(net, 0.2, 0.1, "invalid"), regexp = "initial.infect.method incorrectly specified")
  expect_error(spread.infection(net, 0.2, 0.1, "nodes_fixed", options = list(num.infect = 10)), regexp = "options\\$num.infect must be between 1 and the number of nodes in the network")
  expect_error(spread.infection(net, 0.2, 0.1, "nodes_fixed", options = list(num.infect = 0)), regexp = "options\\$num.infect must be between 1 and the number of nodes in the network")

  # Test that number of infected nodes in Z0 matches given number, if number of infected nodes is fixed ahead of time
  result <- spread.infection(net, NULL, 0.1, "nodes_fixed", options = list(num.infect = 3))
  expect_equal(sum(result$Z0), 3)

  # Test that spreading infection with fixed Z0 and W gives expected results
  result <- spread.infection(net, NULL, NULL, "nodes", options = list(Z0.fixed = Z0, W.fixed = W))
  expect_equal(result$Z0, Z0)
  expect_equal(result$Z, Z.expect)
  expect_true(all(as.sociomatrix(result$W.net) == W))
  expect_equal(get.vertex.attribute(result$disease.net, "initial.infection"), Z0)
  expect_equal(get.vertex.attribute(result$disease.net, "infected"), Z.expect)
  expect_equal(get.edge.attribute(result$disease.net, "spread"), c(1, 1, 1, 0))
})

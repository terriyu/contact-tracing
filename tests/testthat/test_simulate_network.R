context("Test simulate_network.R")

source('../../simulate_network.R')

test_that("Test create.edge.prob.mtx", {

  # Test non-square P.ij matrix
  expect_error(create.edge.prob.mtx(c(25,5,25), matrix(runif(10), nrow=2, ncol=5)), regexp = "P.ij not a square matrix")
  # Test non-symmetric P.ij matrix
  expect_error(create.edge.prob.mtx(c(25,5,25), matrix(1:9, nrow=3, ncol=3)), regexp = "P.ij is not symmetric")
  # Test nodes.per.class and P.ij are consistent
  expect_error(create.edge.prob.mtx(c(25,5,25), matrix(rep(1,4), nrow=2, ncol=2)), regexp = "nodes.per.class and P.ij have inconsistent dimensions")

  # Test for case where there is only one node class
  expect_equal(create.edge.prob.mtx(3, 0.5), matrix(c(0, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, 0), nrow = 3, ncol = 3))

  # Test for case where there is more than one node class
  # Number of nodes in each class
  nodes.per.class <- c(1,2,3)
  # Matrix containing probability of link between class i and class j
  P.ij <- matrix(c(0.5, 0.4, 0, 0.4, 0.8, 0.4, 0, 0.4, 0.5), nrow = 3, ncol = 3)
  # Expected matrix
  edge.prob.mtx <- matrix(c(0, 0.4, 0.4, 0, 0, 0, 0.4, 0, 0.8, 0.4, 0.4, 0.4, 0.4, 0.8, 0, 0.4, 0.4, 0.4, 0, 0.4, 0.4, 0, 0.5, 0.5, 0, 0.4, 0.4, 0.5, 0, 0.5, 0, 0.4, 0.4, 0.5, 0.5, 0), nrow = 6, ncol = 6)
  expect_equal(create.edge.prob.mtx(nodes.per.class, P.ij), edge.prob.mtx)
})

test_that("Test generate.network", {
  set.seed(137)
  # Number of nodes in each class
  nodes.per.class <- c(1,2,3)
  # Matrix containing probability of link between class i and class j
  P.ij <- matrix(c(0.5, 0.4, 0, 0.4, 0.8, 0.4, 0, 0.4, 0.5), nrow = 3, ncol = 3)
  net <- generate.network(nodes.per.class, P.ij)
  socio.mtx <- as.sociomatrix(net)

  # Test that all values of sociomatrix are 0 or 1
  expect_true(all((socio.mtx == 0) | (socio.mtx == 1)))

  # Test that sociomatrix is symmetric
  expect_true(all(socio.mtx == t(socio.mtx)))
})

test_that("Test spread.infection", {
  net <- network(matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 6, ncol = 6), directed = FALSE)
  W <- matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 6, ncol = 6)
  Z0 <- c(1,0,0,0,0,0)
  result <- spread.infection(net, NULL, NULL, FALSE, TRUE, Z0, W)
  Z <- c(1, 1, 1, 1, 0, 0)

  # Test that spreading infection with fixed Z0 and W gives expected results
  expect_equal(result$Z0, Z0)
  expect_equal(result$Z, Z)
  expect_true(all(as.sociomatrix(result$W.net) == W))
  expect_equal(get.vertex.attribute(result$disease.net, "initial.infection"), Z0)
  expect_equal(get.vertex.attribute(result$disease.net, "infected"), Z)
  expect_equal(get.edge.attribute(result$disease.net, "spread"), c(1, 1, 1, 0))
})

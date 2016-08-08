context("Test ct_sampling.R")

source('../../ct_sampling.R')

###################################################

test_that("Test .sample.next", {
  set.seed(137)
  # Infected nodes
  Z <- c(1, 1, 0, 0, 1, 0)
  # Sociomatrix Y.ij
  socio.mtx <- matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 6, ncol = 6)

  # Test error checking
  S.prev <- c(0, 0, 1, 0, 0, 0)
  S <- c(0, 0, 1, 1, 1, 0)
  expect_error(.sample.next(S.prev, S, Z, socio.mtx, "invalid_design", NULL, NULL), regexp = "Invalid value specified for ct.design")

  # Test each of the 4 sampling designs
  S.prev <- c(0, 0, 1, 0, 0, 0)
  S <- c(1, 0, 1, 0, 0, 1)
  expect_equal(.sample.next(S.prev, S, Z, socio.mtx, "infected_only", NULL, NULL), c(0, 1, 0, 0, 0, 0))

  S.prev <- c(0, 0, 1, 0, 0, 0)
  S <- c(1, 0, 1, 0, 0, 1)
  expect_equal(.sample.next(S.prev, S, Z, socio.mtx, "infected_and_edge_units", NULL, NULL), c(0, 0, 0, 0, 0, 0))

  S.prev <- c(1, 0, 0, 0, 0, 0)
  S <- c(1, 0, 0, 0, 0, 1)
  expect_equal(.sample.next(S.prev, S, Z, socio.mtx, "infected_and_edge_units", NULL, NULL), c(0, 0, 1, 0, 0, 0))

  S.prev <- c(0, 0, 1, 0, 0, 0)
  S <- c(1, 0, 1, 0, 0, 1)
  expect_equal(.sample.next(S.prev, S, Z, socio.mtx, "contacts_of_edge_units", NULL, NULL), c(0, 0, 0, 0, 0, 0))

  S.prev <- c(1, 0, 0, 0, 0, 0)
  S <- c(1, 0, 0, 0, 0, 1)
  expect_equal(.sample.next(S.prev, S, Z, socio.mtx, "contacts_of_edge_units", NULL, NULL), c(0, 0, 1, 0, 0, 0))

  S.prev <- c(0, 0, 1, 0, 0, 0)
  S <- c(1, 0, 1, 0, 0, 1)
  expect_equal(.sample.next(S.prev, S, Z, socio.mtx, "full_contact_components", NULL, NULL), c(0, 1, 0, 1, 0, 0))

  # Test p.trace.infected and p.trace.uninfected
  S.prev <- c(0, 0, 1, 0, 0, 0)
  S <- c(1, 0, 1, 0, 0, 1)
  expect_equal(.sample.next(S.prev, S, Z, socio.mtx, "full_contact_components", 1, 1), c(0, 1, 0, 1, 0, 0))

  S.prev <- c(0, 0, 1, 0, 0, 0)
  S <- c(1, 0, 1, 0, 0, 1)
  expect_equal(.sample.next(S.prev, S, Z, socio.mtx, "full_contact_components", 0, 0), c(0, 0, 0, 0, 0, 0))

  S.prev <- c(0, 0, 1, 0, 0, 0)
  S <- c(1, 0, 1, 0, 0, 1)
  expect_equal(.sample.next(S.prev, S, Z, socio.mtx, "full_contact_components", 1, 0), c(0, 1, 0, 0, 0, 0))

  S.prev <- c(0, 0, 1, 0, 0, 0)
  S <- c(1, 0, 1, 0, 0, 1)
  expect_equal(.sample.next(S.prev, S, Z, socio.mtx, "full_contact_components", 0, 1), c(0, 0, 0, 1, 0, 0))
})

test_that("Test ct.sample", {
  set.seed(137)
  # Disease network
  disease.net <- network(matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 6, ncol = 6), directed = FALSE)
  # Initial infected nodes
  Z0 <- c(1, 0, 0, 0, 0, 0)
  # Final infected nodes
  Z <- c(1, 1, 1, 1, 0, 0)
  # Set infection attributes in network
  set.vertex.attribute(disease.net, "initial.infection", Z0)
  set.vertex.attribute(disease.net, "infected", Z)

  # Test error checking
  expect_error(ct.sample(disease.net, -1, "infected_only"), regexp = "Bernoulli process parameter sigma must be between 0 and 1")
  expect_error(ct.sample(disease.net, 99, "infected_only"), regexp = "Bernoulli process parameter sigma must be between 0 and 1")
  expect_error(ct.sample(disease.net, 0.2, "infected_only", NULL, NULL, -2), regexp = "p.trace.infected must be between 0 and 1")
  expect_error(ct.sample(disease.net, 0.2, "infected_only", NULL, NULL, 256), regexp = "p.trace.infected must be between 0 and 1")
  expect_error(ct.sample(disease.net, 0.2, "infected_only", NULL, NULL, 0.5, 137), regexp = "p.trace.uninfected must be between 0 and 1")
  expect_error(ct.sample(disease.net, 0.2, "infected_only", NULL, NULL, 0.5, -5), regexp = "p.trace.uninfected must be between 0 and 1")
  expect_error(ct.sample(disease.net, 0.2, "infected_only", NULL, NULL, 0.5, NULL), regexp = "Both p.trace.infected and p.trace.uninfected need to be specified")
  expect_error(ct.sample(disease.net, 0.2, "infected_only", NULL, NULL, NULL, 0.5), regexp = "Both p.trace.infected and p.trace.uninfected need to be specified")
  expect_error(ct.sample(disease.net, 0.2, "infected_only", 0), regexp = "size.S0 must be between 1 and number of nodes")
  expect_error(ct.sample(disease.net, 0.2, "infected_only", 10), regexp = "size.S0 must be between 1 and number of nodes")
  expect_error(ct.sample(disease.net, 0.2, "infected_only", NULL, 0), regexp = "Number of waves needs to be at least 1")
  expect_error(ct.sample(disease.net, 0.2, "invalid_design"), regexp = "Invalid value specified for ct.design")

  # Test size.S0
  result <- ct.sample(disease.net, 0.2, "contacts_of_edge_units", 2)
  expect_equal(sum(result$S0), 2)

  # Test num.waves
  S0 <- c(0, 0, 0, 0, 0, 1)
  result <- ct.sample(disease.net, NULL, "full_contact_components", NULL, 100, NULL, NULL, FALSE, S0)
  expect_equal(result$S, c(0, 0, 0, 0, 0, 1))

  S0 <- c(1, 0, 0, 0, 0, 1)
  result <- ct.sample(disease.net, NULL, "full_contact_components", NULL, 1, NULL, NULL, FALSE, S0)
  expect_equal(result$S, c(1, 0, 1, 0, 0, 1))

  # Test "infected_only" design
  S0 <- c(1, 1, 1, 1, 1, 1)
  result <- ct.sample(disease.net, NULL, "infected_only", NULL, NULL, NULL, NULL, FALSE, S0)
  expect_equal(result$S0, Z)

  S0 <- c(0, 0, 0, 1, 0, 0)
  result <- ct.sample(disease.net, NULL, "infected_only", NULL, NULL, NULL, NULL, FALSE, S0)
  expect_equal(result$S, c(1, 1, 1, 1, 0, 0))

  # Test "infected_and_edge_units" design
  S0 <- c(0, 0, 0, 0, 0, 1)
  result <- ct.sample(disease.net, NULL, "infected_and_edge_units", NULL, NULL, NULL, NULL, FALSE, S0)
  expect_equal(result$S0, S0)
  expect_equal(result$S, S0)

  S0 <- c(0, 0, 0, 1, 0, 0)
  result <- ct.sample(disease.net, NULL, "infected_and_edge_units", NULL, NULL, NULL, NULL, FALSE, S0)
  expect_equal(result$S0, S0)
  expect_equal(result$S, c(1, 1, 1, 1, 1, 0))

  # Test "contacts_of_edge_units" design
  S0 <- c(0, 0, 0, 0, 0, 1)
  result <- ct.sample(disease.net, NULL, "contacts_of_edge_units", NULL, NULL, NULL, NULL, FALSE, S0)
  expect_equal(result$S0, S0)
  expect_equal(result$S, S0)

  S0 <- c(0, 0, 0, 1, 0, 0)
  result <- ct.sample(disease.net, NULL, "contacts_of_edge_units", NULL, NULL, NULL, NULL, FALSE, S0)
  expect_equal(result$S0, S0)
  expect_equal(result$S, c(1, 1, 1, 1, 1, 0))

  # Test "full_contact_components" design
  S0 <- c(0, 0, 0, 0, 0, 1)
  result <- ct.sample(disease.net, NULL, "full_contact_components", NULL, NULL, NULL, NULL, FALSE, S0)
  expect_equal(result$S0, S0)
  expect_equal(result$S, S0)

  S0 <- c(0, 0, 0, 0, 1, 0)
  result <- ct.sample(disease.net, NULL, "full_contact_components", NULL, NULL, NULL, NULL, FALSE, S0)
  expect_equal(result$S0, S0)
  expect_equal(result$S, c(1, 1, 1, 1, 1, 0))

  S0 <- c(0, 0, 0, 0, 1, 0)
  result <- ct.sample(disease.net, NULL, "full_contact_components", NULL, NULL, 1, 1, FALSE, S0)
  expect_equal(result$S0, S0)
  expect_equal(result$S, c(1, 1, 1, 1, 1, 0))

  S0 <- c(0, 0, 0, 0, 1, 0)
  result <- ct.sample(disease.net, NULL, "full_contact_components", NULL, NULL, 0, 1, FALSE, S0)
  expect_equal(result$S0, S0)
  expect_equal(result$S, S0)

  S0 <- c(0, 0, 0, 0, 1, 0)
  result <- ct.sample(disease.net, NULL, "full_contact_components", NULL, NULL, 1, 0, FALSE, S0)
  expect_equal(result$S0, S0)
  expect_equal(result$S, c(1, 1, 1, 1, 1, 0))
})

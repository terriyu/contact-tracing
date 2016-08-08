# Contact tracing sampling

# Sampling designs are in Chapter 4.6 of Krista Gile's PhD thesis

# Load network package without all the startup messages
suppressMessages(library(network))

.sample.next <- function(S.prev, S, Z, socio.mtx, ct.design, p.trace.infected, p.trace.uninfected) {
  # Return next sample in contact tracing process, given sample from previous wave and current sample
  #
  # NOTE: This is a private helper function for ct.sample()
  #
  # Args:
  #   S.prev - sample from previous wave as vector
  #   S - current sample as vector
  #   Z - infected nodes as vector
  #   socio.mtx - adjacency matrix corresponding to socio matrix
  #   ct.design - contact tracing design to use: "infected_only", "infected_and_edge_units", "contacts_of_edge_units", "full_contact_components"
  #   p.trace.infected - Bernoulli parameter for tracing infected contact
  #   p.trace.uninfected - Bernoulli parameter for tracing uninfected contact
  #
  # Returns:
  #   S.next - next sample as vector

  # Links can be followed with complete accuracy
  transmission.mtx <- socio.mtx
  # Links cannot be followed with complete accuracy
  if ((! is.null(p.trace.infected)) & (! is.null(p.trace.uninfected))) {
    num.nodes <- length(Z)
    num.infect <- sum(Z)
    # Degrade links to uninfected nodes with Bernoulli parameter p.trace.uninfected
    transmission.mtx[ , Z == 0] <- transmission.mtx[ , Z == 0] * matrix(rbinom((num.nodes - num.infect) * num.nodes, 1, p.trace.uninfected), nrow = num.nodes, ncol = num.nodes - num.infect)
    # Degrade links to infected nodes with Bernoulli parameter p.trace.infected
    transmission.mtx[ , Z == 1] <- transmission.mtx[ , Z == 1] * matrix(rbinom(num.infect * num.nodes, 1, p.trace.infected), nrow = num.nodes, ncol = num.infect)
  }

  # Compute next sample
  if (ct.design == "infected_only") {
    # Only sample infected nodes
    S.next <- (S.prev %*% transmission.mtx) * Z * (1 - S)
    S.next <- as.vector(S.next)
    S.next[S.next > 0] <- 1
  } else if ((ct.design == "infected_and_edge_units") | (ct.design == "contacts_of_edge_units")) {
    # Sample both infected and uninfected nodes but only trace infected nodes
    S.next <- ((Z * S.prev) %*% transmission.mtx) * (1 - S)
    S.next <- as.vector(S.next)
    S.next[S.next > 0] <- 1
  } else if (ct.design == "full_contact_components") {
    # Sample and trace both infected and uninfected nodes
    S.next <- (S.prev %*% transmission.mtx) * (1 - S)
    S.next <- as.vector(S.next)
    S.next[S.next > 0] <- 1
  } else {
    stop("Invalid value specified for ct.design")
  }

  # Return next sample
  return(S.next)
}

ct.sample <- function(disease.net, sigma, ct.design, size.S0 = NULL, num.waves = NULL, p.trace.infected = NULL, p.trace.uninfected = NULL, verbose = FALSE, S0.fixed = NULL) {
  # Perform contact tracing sample of an infected network
  #
  # Args:
  #   disease.net - network with infection attributes
  #   sigma - Bernoulli parameter for initial sampling process
  #   ct.design - contact tracing design to use: "infected_only", "infected_and_edge_units", "contacts_of_edge_units", "full_contact_components"
  #   size.S0 - [optional] fix size of initial sample to avoid returning empty sample
  #   num.waves - [optional] number of waves to use in sampling (NULL if you keep sampling until no new nodes are sampled)
  #   p.trace.infected - [optional] Bernoulli parameter for tracing infected contact
  #   p.trace.uninfected - [optional] Bernoulli parameter for tracing uninfected contact
  #   verbose - [optional] flag to print variables during sampling process (for debugging)
  #   S0.fixed - [optional] nodes in initial sample (only use for debugging/testing)
  #
  # Returns:
  #   S0 - initial nodes sampled as vector
  #   S - full contact tracing sample as vector

  # Error checking for sigma, p.trace.infected, and p.trace.uninfected
  if (! is.null(sigma)) {
    if ((sigma < 0) | (sigma > 1)) {
      stop("Bernoulli process parameter sigma must be between 0 and 1")
    }
  }
  if (! is.null(p.trace.infected)) {
    if ((p.trace.infected < 0) | (p.trace.infected > 1)) {
      stop("p.trace.infected must be between 0 and 1")
    }
  }
  if (! is.null(p.trace.uninfected)) {
    if ((p.trace.uninfected < 0) | (p.trace.uninfected > 1)) {
      stop("p.trace.uninfected must be between 0 and 1")
    }
  }
  if (xor(is.null(p.trace.infected), is.null(p.trace.uninfected))) {
    stop("Both p.trace.infected and p.trace.uninfected need to be specified")
  }

  # Number of nodes in network
  num.nodes <- network.size(disease.net)
  # Extract infected nodes in network
  Z <- get.vertex.attribute(disease.net, "infected")
  # Extract adjacency matrix for disease network
  socio.mtx <- as.sociomatrix(disease.net)

  # Initial sample
  if (! is.null(S0.fixed)) {
    # Initial sample is completely specified ahead of time
    S0 <- S0.fixed
  } else if (! is.null(size.S0)) {
    if ((size.S0 < 1) | (size.S0 > num.nodes)) {
      stop("size.S0 must be between 1 and number of nodes")
    }
    # Fixed number of nodes in initial sample
    S0 <- rep(0, num.nodes)
    S0[sample(1:num.nodes, size.S0)] <- 1
  } else {
    # Homogenous Bernoulli sample with parameter sigma
    S0 <- rbinom(num.nodes, 1, sigma)
  }

  # For "infected_only" ct.design, only keep infected nodes from initial sample
  if (ct.design == "infected_only") {
    # Only keep infected nodes
    S0 <- S0 * Z
  }

  # Current sample
  S <- S0
  # Previous wave of samples
  S.prev <- S0

  if (verbose) {
    cat("S0:\n")
    print(S0)
  }

  if (is.null(num.waves)) {
    # Keep sampling until no new nodes sampled
    while (sum(S.prev) > 0) {
      S.next <- .sample.next(S.prev, S, Z, socio.mtx, ct.design, p.trace.infected, p.trace.uninfected)

      if (verbose) {
        cat("S.next:\n")
        print(S.next)
      }

      # Update samples
      S <- S + S.next
      S.prev <- S.next
    }
  } else {
    if (num.waves < 1) {
      stop("Number of waves needs to be at least 1")
    }

    # Sample for num.waves
    for (i in 1:num.waves) {
      S.next <- .sample.next(S.prev, S, Z, socio.mtx, ct.design, p.trace.infected, p.trace.uninfected)

      if (verbose) {
        cat("S.next:\n")
        print(S.next)
      }

      # Update samples
      S <- S + S.next
      S.prev <- S.next
    }
  }

  return(list(S0 = S0, S = S))
}

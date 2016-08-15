# subgraphs with size 1 Functions for fitting observed networks and sampled nodes

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ! NOTE: We don't estimate sigma, because we use missing-at-random assumption, !
# !     and therefore sampling parameters can be ignored                        !
# !                                                                             !
# ! NOTE: We assume that links can be observed perfectly (no degradation), i.e. !
# !       p.trace.infected = 1, p.trace.uninfected = 1                          !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 # Need for components functions, components() and component.dist()
suppressMessages(library(sna))
 # Need for spanning tree function mst()
suppressMessages(library(ape))

transform.tree.directed <- function(root, tree) {
  # Delete edges going into root node
  tree[ , root] <- 0

  if (sum(tree[root, ] == 1) == 0) {
    # Base case for recursion
    # No children, so return NULL
    return(NULL)
  } else {
    children <- which(tree[root, ] == 1)
    dyads <- lapply(children, function(child) c(root, child))
    for (cc in children) {
      dyads <- c(dyads, transform.tree.directed(cc, tree))
    }
    return(dyads)
  }
}

estimate.initial.params <- function(Y.obs, Z.obs, S, ct.design, class.labels) {
  # Estimate initial model parameters P.ij, eta, and tau for MLE
  #
  # NOTE: We don't estimate sigma, because we use missing-at-random assumption,
  #       and therefore sampling parameters can be ignored
  #
  # Args:
  #   Y.obs - observed network, derived from sample
  #   Z.obs - observed infection, derived from sample
  #   S - contact tracing sample as vector
  #   ct.design - contact tracing design used for sample: "infected_only",
  #               "infected_and_edge_units", "contacts_of_edge_units", "full_contact_components"
  #   class.labels - class labels (string) for each node, as vector of strings
  #
  # Returns:
  #   P.ij0 - square, symmetric matrix containing probabilities; each entry
  #          is the probability of a link between node class i and node class j
  #   eta0 - Bernoulli process parameter for initial infection
  #   tau0 - Bernoulli process parameter for generating transmission matrix

  # FIXME: error checking

  # ---------- ERROR CHECKING ---------- #

  num.nodes <- length(class.labels)
  num.classes <- length(unique(class.labels))

  # Compute design matrix
  design.mtx <- compute.design.mtx(Z.obs, S, ct.design)

  # ---------- INITIAL ESTIMATE FOR P.ij ---------- #

  # Simple calculation if only one class
  if (num.classes == 1) {
    # Vector of ones
    ones <- rep(1, num.nodes)

    # Number of observed links
    num.links.obs <- t(ones) %*% Y.obs %*% ones
    # Number of possible links
    num.links <- t(ones) %*% design.mtx %*% ones
    # Set initial parameter estimate to proportion of observed links
    P.ij0 <- num.links.obs / num.links
  } else {
    # Compute unique labels and number of labels
    unique.labels <- sort(unique(class.labels))
    num.labels <- length(unique.labels)
    # Initialize P.ij matrix
    P.ij0 <- matrix(NA, nrow = num.labels, ncol = num.labels)

    # Compute estimates for entries of P.ij matrix
    for (i in 1:num.labels) {
      for (j in i:num.labels) {
        if (i == j) {
          # Number of nodes with class i
          dim.ii <- sum(class.labels == unique.labels[i])
          # Vector of ones
          ones.ii<- rep(1, dim.ii)

          # Extract sub-matrix of observed network corresponding to class i
          Y.obs.ii<- Y.obs[class.labels == unique.labels[i]]
          # Extract sub-matrix of design matrix corresponding to class i
          design.mtx.ii <- design.mtx[class.labels == unique.labels[i]]

          # Number of observed links
          num.links.obs <- t(ones.ii) %*% Y.obs.ii %*% ones.ii
          # Number of possible links
          num.links <- t(ones.ii) %*% design.mtx.ii %*% ones.ii

          # Set initial parameter estimate to proportion of observed links
          P.ij0[i, i] <- num.links.obs / num.links
        } else {
          # Number of nodes with class i
          dim.ii <- sum(class.labels == unique.labels[i])
          # Vector of ones
          ones.ii <- rep(1, dim.ii)

          # Number of nodes with class j
          dim.jj <- sum(class.labels == unique.labels[j])
          # Vector of ones
          ones.jj <- rep(1, dim.jj)

          # Number of nodes with class i or j
          dim.ij <- sum(class.labels == unique.labels[i]) + sum(class.labels == unique.labels[j])
          # Vector of ones
          ones.ij <- rep(1, dim.ij)

          # Extract sub-matrix of observed network corresponding to class i
          Y.obs.ii <-  Y.obs[class.labels == unique.labels[i]]
          # Extract sub-matrix of observed network corresponding to class j
          Y.obs.jj <-  Y.obs[class.labels == unique.labels[j]]
          # Extract sub-matrix of observed network corresponding to classes i and j
          Y.obs.ij <- Y.obs[class.labels == unique.labels[i] | class.labels == unique.labels[j]]

          # Extract sub-matrix of design matrix corresponding to class i
          design.mtx.ii <- design.mtx[class.labels == unique.labels[i]]
          # Extract sub-matrix of design matrix corresponding to class j
          design.mtx.jj <- design.mtx[class.labels == unique.labels[j]]
          # Extract sub-matrix of design matrix corresponding to classes i and j
          design.mtx.ij <- design.mtx[class.labels == unique.labels[i] | class.labels == unique.labels[j]]

          # Number of observed links
          num.links.obs <- (t(ones.ij) %*% Y.obs.ij %*% ones.ij) - (t(ones.ii) %*% Y.obs.ii %*% ones.ii) - (t(ones.jj) %*% Y.obs.jj %*% ones.jj)
          # Number of possible links
          num.links <- (t(ones.ij) %*% design.mtx.ij %*% ones.ij) - (t(ones.ii) %*% design.mtx.ii %*% ones.ii) - (t(ones.jj) %*% design.mtx.jj %*% ones.jj)
          # Set initial parameter estimate to proportion of observed links
          P.ij0[i, j] <- num.links.obs / num.links
          P.ij0[j, i] <- P.ij0[i, j]
        }
      }
    }
  }

  # Error checking on P.ij0

  if (sum(is.na(P.ij0)) != 0) {
    stop("Error calculating initial parameters: Not all entries of P.ij0 were calculated")
  }

  if (! all(P.ij0 == t(P.ij0))) {
    stop("Error calculating initial parameters: P.ij0 matrix is not symmetric")
  }

  if (! all((P.ij0 >= 0) | (P.ij <= 1))) {
    stop("Error calculating initial parameters: All entries of P.ij0 matrix must be between 0 and 1")
  }
  # ---------- INITIAL ESTIMATE FOR ETA ---------- #

  # Compute number of infected subgraphs

  # QUESTION: Best way to compute number of infected subgraphs for full contact components design?
  Y.obs.infected <- Y.obs
  # Remove all links to observed uninfected nodes
  Y.obs.infected[(Z.obs == 0) & (S == 1), ] <- 0
  Y.obs.infected[ , (Z.obs == 0) & (S == 1)] <- 0

  cd.result <- component.dist(Y.obs.infected)
  cd.cdist <- cd.result$cdist
  cd.membership <- cd.result$membership
  cd.csize <- cd.result$csize

  # Number of infected subgraphs with size greater than 1
  Q <- sum(cd.cdist[-1])
  # Add number of infected subgraphs with size 1
  Q <- Q + sum(cd.csize[cd.membership[Z.obs == 1]] == 1)

  # Initial estimate for eta
  eta0 <- Q / (t(S) %*% ones)

  # Error checking for eta
  if ((eta0 < 0) | (eta0 > 1)) {
    stop("Error calculating initial parameters: eta0 must be between 0 and 1")
  }

  # ---------- INITIAL ESTIMATE FOR TAU ---------- #

  # Initial estimate for tau
  # Number of transmission events which occurred
  A <- t(Z.obs) %*% ones - Q
  # Number of possible transmission events
  B <- t(Z.obs) %*% Y.obs %*% ones
  tau0 <- A / B

  # Error checking for tau
  if ((tau0 < 0) | (tau0 > 1)) {
    stop("Error calculating initial parameters: tau0 must be between 0 and 1")
  }

  # Return estimate of initial model parameters
  return(list(P.ij0 = P.ij0, eta0 = eta0, tau0 = tau0))
}

initial.net.sample <- function(Y.obs, Z.obs, S, class.labels, P.ij0, eta0, tao0) {
  # FIXME: Error checking

  # ---------- ERROR CHECKING ---------- #

  # ---------- PERFORM Y SAMPLE ---------- #

  # Compute number of nodes and classes
  num.nodes <- length(class.labels)
  num.classes <- length(unique(class.labels))

  # Compute design matrix
  design.mtx <- compute.design.mtx(Z.obs, S, ct.design)

  # Boolean matrix of unobserved links in upper triangle
  unobserved.links <- (design.mtx == 0) & upper.tri(design.mtx)
  num.unobserved.links <- sum(unobserved.links)
  unobserved.indices <- which(unobserved.links, arr.ind = TRUE)

  Y <- Y.obs
  # Set lower triangle to 0, will fill in later
  Y[lower.tri(Y)] <- 0

  if (num.classes == 1) {
    # Only one class, so P.ij0 is a scalar
    Y[unobserved.indices] <- rbinom(sum(unobserved.links), 1, P.ij0)
  } else {
    # Case of multiple classes
    for (i in 1:num.unobserved.links) {
      bernoulli.prob <- P.ij0(unobserved.indices[i, 1], unobserved.indices[i, 2])
      Y[unobserved.indices[i, ]] <- rbinom(1, 1, bernoulli.prob)
    }
  }

  # Fill in lower triangle
  Y <- Y + t(Y)

  if (sum(diag(Y)) != 0) {
    stop("Error calculating initial MCMC sample: Diagonal of Y must be zero")
  }

  if (! all(Y == t(Y))) {
    stop("Error calculating initial MCMC sample: Y is not a symmetric matrix")
  }

  if (! all((Y == 0) | (Y == 1))) {
    stop("Error calculating initial MCMC sample: All entries of Y matrix must be 0 or 1")
  }

  # ---------- PERFORM Z0 AND W SAMPLES ---------- #

  # QUESTION: Best way to compute number of infected subgraphs for full contact components design?
  Y.obs.infected <- Y.obs
  # Remove all links to observed uninfected nodes
  Y.obs.infected[(Z.obs == 0) & (S == 1), ] <- 0
  Y.obs.infected[ , (Z.obs == 0) & (S == 1)] <- 0

  # Compute component distribution
  cd.result <- component.dist(Y.obs.infected)
  cd.membership <- cd.result$membership
  num.subgraphs <- components(Y.obs.infected)

  # Initialize Z0
  Z0 <- rep(NA, num.nodes)

  # Initialize the transmissiblity matrix W
  W <- matrix(NA, nrow = num.nodes, ncol = num.nodes)
  # A node can't transmit an infection to itself
  diag(W) <- 0

  # Infect one node in each infected subgraph
  # QUESTION: Could be more efficient to treat one-node subgraphs and larger subgraphs separately?

  # Loop through infected subgraphs
  for (i in 1:num.subgraphs) {
    # Boolean vector for infected nodes in subgraph i
    subgraph.infected.nodes <- (cd.membership == i) & (Z.obs == 1)

    # If subgraph i is infected, set one of its infected nodes to Z0 = 1

    # Check that subgraph has at least one infected node
    if (sum(subgraph.infected.nodes) > 0) {
      # Numerical indices for infected nodes in subgraph i
      subgraph.infected.indices <- which(subgraph.infected.nodes)

      # Sample one of the infected nodes in subgraph i
      infected.idx <- sample(subgraph.infected.indices, 1)

      # Set corresponding node in Z0 to infected status
      Z0[infected.idx] <- 1

      # Case of infected subgraph of size 2 or larger
      if (sum(subgraph.infected.nodes) > 1) {
        # Compute a minimum spanning tree for infected subgraph
        # NOTE: We don't need a minimum spanning tree, any spanning tree is sufficient
        # NOTE: Multiply by -1, since we are computing minimum spanning tree
        spanning.tree <- mst(-1*Y.obs.infected[subgraph.infected.nodes, subgraph.infected.nodes])

        # Turn spanning tree into a directed spanning tree, with root at initial infected node Z0

        # Get index for root
        Z0.root <- which(subgraph.infected.indices == infected.idx)
        # Get list of directed edges corresponding to directed spanning tree
        directed.edges.list <- transform.tree.directed(Z0.root, spanning.tree)
        # Convert list to array indices
        flattened.indices <- sapply(unlist(directed.edges.list), function(idx) subgraph.infected.indices[idx])
        directed.edges.indices <- matrix(flattened.indices, nrow = length(directed.edges.list), ncol = 2, byrow = TRUE)

        # Set W = 1 for all arcs in infected spanning tree
        W[directed.edges.indices] <- 1
      }
    }
  }

  # QUESTION: PhD thesis and Krista's code are different
  # Set values for all other nodes in Z0

  # Any other observed nodes (not already set) are set to 0
  Z0[is.na(Z0) & (S == 1)] <- 0
  # Unobserved nodes randomly infected with Bernoulli process
  Z0[is.na(Z0) & (S == 0)] <- rbinom(sum(is.na(Z0) & (S == 0)), 1, eta0)

  # Set W = 0 for observed links originating from observed infected nodes
  # Only do this, if W has not already been set for links (entry not NA)

  # Matrix of links that should have W set to 0
  # If entry in matrix is 1, that link should be set
  # If entry in matrix is 0, that link should not be set
  # You can think of this matrix as a representation for a set of indicies

  # Initialize all entries to 0
  Y.obs.untransmitted <- matrix(0, nrow = num.nodes, ncol = num.nodes)
  # Links from observed infected nodes to observed uninfected nodes could potentially be set
  Y.obs.untransmitted[Z.obs == 1, (Z.obs == 0) & (S == 1)] <- 1
  # Links from observed infected nodes to other observed infected nodes could potentially be set
  Y.obs.untransmitted[Z.obs == 1, Z.obs == 1] <- 1
  # Eliminate any links that are unobserved or which have already been set in W
  # HACK: R automatically casts (W == NA) from Boolean to integer matrix
  Y.obs.untransmitted <- Y.obs.untransmitted * Y.obs * (W == NA)

  # QUESTION: Krista seems to perform the W sample in a more simplified way
  # Set W = 0 according to Y.obs.untransmitted
  W[Y.obs.untransmitted == 1] <- 0

  # Any unset entries in W are assigned by Bernoulli process
  W[is.na(W)] <- rbinom(sum(is.na(W)), 1, tau0)

  # FIXME: Compute reachabiity

  # Error checking for Z0

  if (sum(is.na(Z0)) > 0) {
    stop("Error calculating initial MCMC sample: Some elements of Z0 were not calculated")
  }

  if (! all((Z0 == 0) | (Z0 == 1))) {
    stop("Error calculating initial MCMC sample: All elements of Z0 must be 0 or 1")
  }

  # Error checking for W

  if (sum(is.na(W)) > 0) {
    stop("Error calculating initial MCMC sample: Some elements of W were not calculated")
  }

  if (sum(diag(Y)) != 0) {
    stop("Error calculating initial MCMC sample: Diagonal of W must be zero")
  }

  if (! all((Y == 0) | (Y == 1))) {
    stop("Error calculating initial MCMC sample: All entries of W matrix must be 0 or 1")
  }

  # Return sample for Y, Z0, and W
  return(list(Y = Y, Z0 = Z0, W = W))
}

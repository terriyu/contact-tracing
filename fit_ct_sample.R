# Functions for fitting contact tracing samples

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ! NOTE: We don't estimate sigma, because we use missing-at-random assumption, !
# !     and therefore sampling parameters can be ignored                        !
# !                                                                             !
# ! NOTE: We assume that links can be observed perfectly (no degradation), i.e. !
# !       p.trace.infected = 1, p.trace.uninfected = 1                          !
# !                                                                             !
# ! NOTE: For observed variables Y.obs and Z.obs, we set an unobserved entry    !
# !       to 0, even if it has not been observed.                               !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# FIXME: Change Y.obs and Z.obs to use NA for unobserved entries

# Need sna package for components functions, components() and component.dist()
# NOTE: Potential conflict with network package for operator %c%
suppressMessages(library(sna))

####################################################################################
############# FUNCTIONS FOR INITIAL PARAMETERS AND INITIAL MCMC SAMPLE #############
####################################################################################

compute.dst.bf <- function(root, subgraph) {
  # Compute directed spanning tree of unweighted connected subgraph
  # Uses breadth-first traversal
  #
  # Args:
  #   root - Index of tree root (integer)
  #   subgraph - adjacency matrix for connected subgraph
  #
  # Returns:
  #   directed.edges.list <- list of edges in directed spanning tree

  # Error checking
  if (dim(subgraph)[1] != dim(subgraph)[2]) {
    stop("Subgraph is not a square matrix")
  }

  if (sum(diag(subgraph)) != 0) {
    stop("Subgraph must have zero diagonal")
  }

  if ((root < 1) | (root > dim(subgraph)[1])) {
    stop("Invalid index for root of tree")
  }

  if (! is.connected(subgraph)) {
    stop("Subgraph is not connected")
  }

  # Compute number of nodes
  num.nodes <- dim(subgraph)[1]

  # Initialize list of directed edges in spanning tree
  directed.edges.list <- list()
  # Initialize Boolean vector indicating which nodes have been visited
  visited <- rep(FALSE, num.nodes)

  # Visit root node
  # Initialize queue of nodes to explore, enqueue root
  queue <- list(root)
  # Set root as visited
  visited[root] <- TRUE

  while (length(queue) > 0) {
    # Dequeue head and set to current node
    current.node <- queue[[1]]
    queue[[1]] <- NULL
    # Find all children of current node
    children <- which(subgraph[current.node, ] == 1)

    # Check all children to see if they have been visited
    for (child in children) {
      # If child has not been visited, add edge and add child to queue
      if (! visited[child]) {
        # Mark child as visited
        visited[child] <- TRUE
        # Enqueue child
        queue <- c(queue, child)

        # Add edge to directed spanning tree
        directed.edges.list <- c(directed.edges.list, list(c(current.node, child)))
      }
    }
  }

  return(directed.edges.list)
}

compute.dst.df <- function(root, subgraph) {
  # Compute directed spanning tree of unweighted connected subgraph
  # Uses depth-first traversal
  #
  # Args:
  #   root - Index of tree root (integer)
  #   subgraph - adjacency matrix for connected subgraph
  #
  # Returns:
  #   directed.edges.list <- list of edges in directed spanning tree

  # Error checking
  if (dim(subgraph)[1] != dim(subgraph)[2]) {
    stop("Subgraph is not a square matrix")
  }

  if (sum(diag(subgraph)) != 0) {
    stop("Infected subgraph must have zero diagonal")
  }

  if ((root < 1) | (root > dim(subgraph)[1])) {
    stop("Invalid index for root of tree")
  }

  if (! is.connected(subgraph)) {
    stop("Subgraph is not connected")
  }

  # Compute number of nodes
  num.nodes <- dim(subgraph)[1]

  # Initialize list of directed edges in spanning tree
  directed.edges.list <- list()
  # Initialize Boolean vector indicating which nodes have been visited
  visited <- rep(FALSE, num.nodes)

  # Visit root node
  # Initialize stack of nodes to explore, push root onto stack
  stack <- list(root)
  # Set root as visited
  visited[root] <- TRUE

  while (length(stack) > 0) {
    # Peek at top of stack
    current.node <- stack[[length(stack)]]
    # Find all unvisited children of current node
    unvisited.children <- which(subgraph[current.node, ] == 1 & (! visited))

    if (length(unvisited.children) != 0) {
      # Get first unvisited child
      next.node <- unvisited.children[1]
      # Mark it as visited
      visited[next.node] <- TRUE
      # Push it onto the stack
      stack <- c(stack, next.node)

      # Add edge to directed spanning tree
      directed.edges.list <- c(directed.edges.list, list(c(current.node, next.node)))
    } else {
      # All ancestors of current node visited, pop from stack
      stack[[length(stack)]] <- NULL
    }
  }

  return(directed.edges.list)
}

estimate.initial.params <- function(obs, class.labels) {
  # Estimate initial model parameters P.ij, eta, and tau for MLE
  #
  # NOTE: We don't estimate sigma, because we make the missing-at-random assumption,
  #       and therefore sampling parameters can be ignored
  #
  # Args:
  #   [obs is a list of the observations]
  #   obs$Y.obs - observed network, derived from sample
  #   obs$Z.obs - observed infection, derived from sample
  #   obs$S - contact tracing sample as vector
  #   obs$ct.design - contact tracing design used for sample: "infected_only",
  #                   "infected_and_edge_units", "contacts_of_edge_units", "full_contact_components"
  #
  #   class.labels - class label (string) for each node, as vector of strings;
  #                  each label must be in 1:num.nodes
  #
  # Returns:
  #   P.ij - square, symmetric matrix containing probabilities; each entry
  #          is the probability of a link between node class i and node class j
  #   eta - Bernoulli process parameter for initial infection
  #   tau - Bernoulli process parameter for generating transmission matrix

  # ---------- ERROR CHECKING ---------- #

  dims <- c(dim(obs$Y.obs)[1], dim(obs$Y.obs)[2], length(obs$Z.obs), length(class.labels))
  if (! all(dims == length(obs$S))) {
    stop("Number of nodes inconsistent for one of obs$Y.obs, obs$Z.obs, obs$S, class.labels")
  }

  if (! (all(obs$Y.obs == 0 | obs$Y.obs == 1))) {
    stop("obs$Y.obs needs to be a matrix of 0's and 1's")
  }

  if (! (all(obs$Z.obs == 0 | obs$Z.obs == 1))) {
    stop("obs$Z.obs needs to be a vector of 0's and 1's")
  }

  if (! (all(obs$S == 0 | obs$S == 1))) {
    stop("obs$S needs to be a vector of 0's and 1's")
  }

  if (dim(obs$Y.obs)[1] != dim(obs$Y.obs)[2]) {
    stop("obs$Y.obs is not a square matrix")
  }

  if (sum(diag(obs$Y.obs)) != 0) {
    stop("obs$Y.obs must have zero diagonal")
  }

  # ----------- CALCULATE USEFUL VARIABLES ---------- #

  # Compute number of nodes and classes
  num.nodes <- length(class.labels)
  num.classes <- length(unique(class.labels))

  # Compute design matrix
  design.mtx <- compute.design.mtx(obs$Z.obs, obs$S, obs$ct.design)

  # ---------- INITIAL ESTIMATE FOR P.ij ---------- #

  # NOTE: Factor of 1/2 in num.links.obs and num.links is due to the fact
  #       that Y is undirected

  if (num.classes == 1) {
    # Simple calculation if only one class

    # Construct a vector of ones for convenience
    ones <- rep(1, num.nodes)

    # Number of observed links
    num.links.obs <- 1/2 * (t(ones) %*% obs$Y.obs %*% ones)
    # Number of possible links
    num.links <- 1/2 * (t(ones) %*% design.mtx %*% ones)
    # Set initial parameter estimate to proportion of observed links
    P.ij0 <- num.links.obs / num.links
  } else {
    # Case of multiple classes

    # Initialize P.ij matrix
    P.ij0 <- matrix(NA, nrow = num.classes, ncol = num.classes)

    # Compute estimates for entries of P.ij matrix
    for (i in 1:num.labels) {
      for (j in i:num.labels) {
        if (i == j) {
          # Number of nodes with class i
          dim.ii <- sum(class.labels == i)
          # Vector of ones
          ones.ii <- rep(1, dim.ii)

          # Extract sub-matrix of observed network corresponding to class i
          Y.obs.ii <- obs$Y.obs[class.labels == i, class.labels == i]
          # Extract sub-matrix of design matrix corresponding to class i
          design.mtx.ii <- design.mtx[class.labels == i, class.labels == i]

          # Number of observed links
          num.links.obs <- 1/2 * (t(ones.ii) %*% Y.obs.ii %*% ones.ii)
          # Number of possible links
          num.links <- 1/2 * (t(ones.ii) %*% design.mtx.ii %*% ones.ii)

          # Set initial parameter estimate to proportion of observed links
          P.ij0[i, i] <- num.links.obs / num.links
        } else {
          # Number of nodes with class i
          dim.ii <- sum(class.labels == i)
          # Vector of ones
          ones.ii <- rep(1, dim.ii)

          # Number of nodes with class j
          dim.jj <- sum(class.labels == j)
          # Vector of ones
          ones.jj <- rep(1, dim.jj)

          # Number of nodes with class i or j
          dim.ij <- sum(class.labels == i) + sum(class.labels == j)
          # Vector of ones
          ones.ij <- rep(1, dim.ij)

          # Extract sub-matrix of observed network corresponding to class i
          Y.obs.ii <-  obs$Y.obs[class.labels == i, class.labels == i]
          # Extract sub-matrix of observed network corresponding to class j
          Y.obs.jj <-  obs$Y.obs[class.labels == j, class.labels == j]
          # Extract sub-matrix of observed network corresponding to classes i and j
          Y.obs.ij <- obs$Y.obs[class.labels == i | class.labels == j, class.labels == i | class.labels == j]

          # Extract sub-matrix of design matrix corresponding to class i
          design.mtx.ii <- design.mtx[class.labels == i, class.labels == i]
          # Extract sub-matrix of design matrix corresponding to class j
          design.mtx.jj <- design.mtx[class.labels == j, class.labels == j]
          # Extract sub-matrix of design matrix corresponding to classes i and j
          design.mtx.ij <- design.mtx[class.labels == i | class.labels == j, class.labels == i | class.labels == j]

          # Number of observed links
          num.links.obs <- (t(ones.ij) %*% Y.obs.ij %*% ones.ij) - (t(ones.ii) %*% Y.obs.ii %*% ones.ii) - (t(ones.jj) %*% Y.obs.jj %*% ones.jj)
          num.links.obs <- 1/2 * num.links.obs
          # Number of possible links
          num.links <- (t(ones.ij) %*% design.mtx.ij %*% ones.ij) - (t(ones.ii) %*% design.mtx.ii %*% ones.ii) - (t(ones.jj) %*% design.mtx.jj %*% ones.jj)
          num.links <- 1/2 * num.links
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

  if (! all((P.ij0 >= 0) | (P.ij0 <= 1))) {
    stop("Error calculating initial parameters: All entries of P.ij0 matrix must be between 0 and 1")
  }
  # ---------- INITIAL ESTIMATE FOR ETA ---------- #

  # Compute number of infected subgraphs

  Y.obs.infected <- obs$Y.obs
  # Remove all links to observed uninfected nodes
  Y.obs.infected[(obs$Z.obs == 0) & (obs$S == 1), ] <- 0
  Y.obs.infected[ , (obs$Z.obs == 0) & (obs$S == 1)] <- 0

  cd.result <- component.dist(Y.obs.infected)
  cd.cdist <- cd.result$cdist
  cd.membership <- cd.result$membership
  cd.csize <- cd.result$csize

  # Number of infected subgraphs with size greater than 1
  Q <- sum(cd.cdist[-1])
  # Add number of infected subgraphs with size 1
  Q <- Q + sum(cd.csize[cd.membership[obs$Z.obs == 1]] == 1)

  # Initial estimate for eta
  num.samples <- t(S) %*% ones
  if (num.samples < 1) {
    stop("Error: Sample size of zero")
  }
  eta0 <- Q / num.samples

  # Error checking for eta
  if ((eta0 < 0) | (eta0 > 1)) {
    stop("Error calculating initial parameters: eta0 must be between 0 and 1")
  }

  # ---------- INITIAL ESTIMATE FOR TAU ---------- #

  # Initial estimate for tau
  # Number of transmission events which occurred
  A <- t(obs$Z.obs) %*% ones - Q
  # Number of possible transmission events
  B <- t(obs$Z.obs) %*% obs$Y.obs %*% ones


  tau0 <- A / B
  if (B == 0) {
    stop("Error: zero possible transmission events, divide by zero")
  }

  # Error checking for tau
  if ((tau0 < 0) | (tau0 > 1)) {
    stop("Error calculating initial parameters: tau0 must be between 0 and 1")
  }

  # Return estimate of initial model parameters
  return(list(P.ij = P.ij0, eta = eta0, tau = tau0))
}

draw.initial.sample <- function(obs, class.labels, params) {
  # Draw initial MCMC sample based on initial model parameters
  #
  # NOTE: params are in members of a list, e.g params = list(P.ij = P.ij, eta = eta, tau = tau)
  #
  # Args:
  #   [obs is a list of observations]
  #   obs$Y.obs - observed network, derived from sample
  #   obs$Z.obs - observed infection, derived from sample
  #   obs$S - contact tracing sample as vector
  #   obs$ct.design - contact tracing design used for sample: "infected_only",
  #                   "infected_and_edge_units", "contacts_of_edge_units", "full_contact_components"
  #
  #   class.labels - class labels (string) for each node, as vector of strings;
  #                  each label must be in 1:num.nodes
  #
  #   [params is a list of model parameters]
  #   params$P.ij - square, symmetric matrix containing probabilities; each entry
  #                  is the probability of a link between node class i and node class j
  #   params$eta - Bernoulli process parameter for initial infection
  #   params$tau - Bernoulli process parameter for generating transmission matrix
  #
  # Returns:
  #   Y - network sample (matrix)
  #   Z0 - initial infection sample (vector)
  #   W - transmissibility matrix sample (matrix)

  # ---------- ERROR CHECKING ---------- #

  # Compute dimensions of input variables
  dims <- c(dim(obs$Y.obs)[1], dim(obs$Y.obs)[2], length(obs$Z.obs), length(class.labels))
  if (! all(dims == length(S))) {
    stop("Number of nodes inconsistent for one of obs$Y.obs, obs$Z.obs, obs$S, class.labels")
  }

  if (! (all(obs$Y.obs == 0 | obs$Y.obs == 1))) {
    stop("obs$Y.obs needs to be a matrix of 0's and 1's")
  }

  if (! (all(obs$Z.obs == 0 | obs$Z.obs == 1))) {
    stop("obs$Z.obs needs to be a vector of 0's and 1's")
  }

  if (! (all(obs$S == 0 | obs$S == 1))) {
    stop("obs$S needs to be a vector of 0's and 1's")
  }

  if (dim(obs$Y.obs)[1] != dim(obs$Y.obs)[2]) {
    stop("obs$Y.obs is not a square matrix")
  }

  if (sum(diag(obs$Y.obs)) != 0) {
    stop("obs$Y.obs must have zero diagonal")
  }

  if (! is.null(params$eta)) {
    if ((params$eta < 0) | (params$eta > 1)) {
      stop("Bernoulli process parameter eta must be between 0 and 1")
    }
  }

  if (! is.null(params$tau)) {
    if ((params$tau < 0) | (params$tau > 1)) {
      stop("Bernoulli process parameter tau must be between 0 and 1")
    }
  }

  # Error checking: check that params$P.ij is square matrix
  if (length(params$P.ij) > 1) {
    if (dim(params$P.ij)[1] != dim(params$P.ij)[2]) {
      stop("params$P.ij not a square matrix")
    }
  }

  # Error checking: check that params$P.ij is symmetric matrix
  # Probability of link between class i and class j should be same as
  # probability of link between class j and class i
  if (! all(params$P.ij == t(params$P.ij))) {
    stop("params$P.ij is not symmetric")
  }

  # Error checking: check that class.labels and params$P.ij are consistent
  if (length(params$P.ij) > 1) {
    # Number of classes
    num.classes <- length(unique(class.labels))
    if (num.classes != dim(params$P.ij)[1]) {
      stop("class.labels and params$P.ij are inconsistent")
    }
  }

  # ---------- PERFORM Y SAMPLE ---------- #

  # Compute number of nodes and classes
  num.nodes <- length(class.labels)
  num.classes <- length(unique(class.labels))

  # Compute design matrix
  design.mtx <- compute.design.mtx(obs$Z.obs, obs$S, obs$ct.design)

  # Boolean matrix of unobserved links in upper triangle
  unobserved.links <- (design.mtx == 0) & upper.tri(design.mtx)
  num.unobserved.links <- sum(unobserved.links)
  unobserved.indices <- which(unobserved.links, arr.ind = TRUE)

  Y <- obs$Y.obs
  # Set lower triangle to 0, will fill in later
  Y[lower.tri(Y)] <- 0

  if (num.classes == 1) {
    # Only one class, so params$P.ij is a scalar
    Y[unobserved.indices] <- rbinom(sum(unobserved.links), 1, params$P.ij)
  } else {
    # Case of multiple classes
    for (i in 1:num.unobserved.links) {
      bernoulli.prob <- params$P.ij(unobserved.indices[i, 1], unobserved.indices[i, 2])
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

  Y.obs.infected <- obs$Y.obs
  # Remove all links to observed uninfected nodes
  Y.obs.infected[(obs$Z.obs == 0) & (obs$S == 1), ] <- 0
  Y.obs.infected[ , (obs$Z.obs == 0) & (obs$S == 1)] <- 0

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
    subgraph.infected.nodes <- (cd.membership == i) & (obs$Z.obs == 1)

    # If subgraph i is infected, set one of its infected nodes to Z0 = 1

    # Check that subgraph has at least one infected node
    if (sum(subgraph.infected.nodes) > 0) {
      # Numerical indices for infected nodes in subgraph i
      subgraph.infected.indices <- which(subgraph.infected.nodes)

      # Sample one of the infected nodes in subgraph i
      infected.idx <- sample(subgraph.infected.indices, 1)

      # Set corresponding node in Z0 to infected status
      Z0[infected.idx] <- 1

      # QUESTION: Better to construct W using directed spanning tree for infected subgraph
      #           or simply construct W by setting W = 1 for all edges in infected subgraph?

      # Case of infected subgraph of size 2 or larger
      if (sum(subgraph.infected.nodes) > 1) {
        # Root at initial infected node Z0
        # Get index for root
        Z0.root <- which(subgraph.infected.indices == infected.idx)

        # Compute directed spanning tree for infected subgraph

        # NOTE: Prefer breadth-first versus depth-first spanning tree
        #       Makes more sense for practical situations like needle/drug sharing party?
        directed.edges.list <- compute.dst.bf(Z0.root, Y.obs.infected[subgraph.infected.nodes, subgraph.infected.nodes])

        # Convert list to array indices
        flattened.indices <- sapply(unlist(directed.edges.list), function(idx) subgraph.infected.indices[idx])
        directed.edges.indices <- matrix(flattened.indices, nrow = length(directed.edges.list), ncol = 2, byrow = TRUE)

        # Set W = 1 for all arcs in infected spanning tree
        W[directed.edges.indices] <- 1
      }
    }
  }

  # NOTE: PhD thesis and Krista's code are different

  # Set values for all other nodes in Z0

  # Any other observed nodes (not already set) are set to 0
  Z0[is.na(Z0) & (obs$S == 1)] <- 0
  # Unobserved nodes randomly infected with Bernoulli process
  Z0[is.na(Z0) & (obs$S == 0)] <- rbinom(sum(is.na(Z0) & (obs$S == 0)), 1, params$eta)

  ##########################################################################
  # QUESTION: alternative to the above, is this better???
  # Randomly set Z0 for observed infected nodes
  #Z0[is.na(Z0) & (obs$Z.obs == 1)] <- rbinom(sum(is.na(Z0) & (obs$Z.obs == 1)), 1, params$eta)
  # Set Z0 = 0 for observed uninfected nodes
  #Z0[is.na(Z0) & (obs$Z.obs == 0) & (obs$S == 1)] <- 0
  # Unobserved nodes randomly infected with Bernoulli process
  #Z0[is.na(Z0) & (obs$S == 0)] <- rbinom(sum(is.na(Z0) & (obs$S == 0)), 1, params$eta)
  ##########################################################################

  # NOTE: Krista's code performs the W sample in a different way

  # Set W = 0 for observed links originating from observed infected nodes
  # Only do this, if W has not already been set for links (entry not NA)

  # Matrix of links that should have W set to 0
  # If entry in matrix is 1, that link should be set
  # If entry in matrix is 0, that link should not be set
  # You can think of this matrix as a representation for a set of indicies

  # Initialize all entries to 0
  Y.obs.untransmitted <- matrix(0, nrow = num.nodes, ncol = num.nodes)
  # Links from observed infected nodes to observed uninfected nodes could potentially be set
  Y.obs.untransmitted[obs$Z.obs == 1, (obs$Z.obs == 0) & (obs$S == 1)] <- 1
  # Links from observed infected nodes to other observed infected nodes could potentially be set
  # QUESTION: Might delete line below, if we follow Krista's code?
  Y.obs.untransmitted[obs$Z.obs == 1, obs$Z.obs == 1] <- 1
  # Eliminate any links that are unobserved (whether or not they actually exist)
  # or which have already been set in W
  # HACK: R automatically casts (W == NA) from Boolean to integer matrix
  Y.obs.untransmitted <- Y.obs.untransmitted * obs$Y.obs * (W == NA)
#FIXME???
  # Set W = 0 according to Y.obs.untransmitted
  W[Y.obs.untransmitted == 1] <- 0

  # QUESTION: Should links from unsampled nodes to observed uninfected nodes
  #           have W = 0?
  # Any unset entries in W are assigned by Bernoulli process
  W[is.na(W)] <- rbinom(sum(is.na(W)), 1, params$tau)

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

  # Check that Z0 and W from MCMC sample is consistent with obs$Z.obs

  # Compute reachability matrix
  R.D <- reachability(W * Y)
  # Compute resulting infection
  Z <- as.vector(Z0 %*% R.D)
  Z[Z > 1] <- 1

  success <- TRUE
  if (! all(Z[obs$S == 1] == obs$Z.obs[obs$S == 1])) {
    success <- FALSE
  }

  # Set sample for Y, Z0, and W if success
  if (success) {
    sample <- list(Y = Y, Z0 = Z0, W = W)
  } else {
    sample <- NA
  }

  return(list(success = success, sample = sample))
}

############### LIKELIHOOD ###############

log.likelihood <- function(sample, params, class.labels) {
  # Compute log-likelihood of a sample
  #
  # NOTE: Does not check if sample is consistent with observations
  #
  # Args:
  #   [sample is a MCMC sample, a list of network variables]
  #   sample$Y - network sample (matrix)
  #   sample$Z0 - initial infection sample (vector)
  #   sample$W - transmissibility matrix sample (matrix)
  #
  #   [params is a list of model parameters]
  #   params$P.ij - square, symmetric matrix containing probabilities; each entry
  #                  is the probability of a link between node class i and node class j
  #   params$eta - Bernoulli process parameter for initial infection
  #   params$tau - Bernoulli process parameter for generating transmission matrix
  #
  #   class.labels - class labels (string) for each node, as vector of strings;
  #                  each label must be in 1:num.nodes
  #
  # Returns:
  #   loglik - log-likelihood of sample

  ################################################################
  # NOTE: Assumes that sample is consistent with observations!!! #
  ################################################################

  # Compute number of nodes
  num.nodes <- length(sample$Z0)
  # Construct a vector of ones for convenience
  ones <- rep(1, num.nodes)

  # Part of log-likelihood due to Z0
  loglik.Z0 <- (t(sample$Z0) %*% ones) * log(params$eta / (1 - params$eta)) + num.nodes * log(1 - params$eta)

  # Part of log-likelihood due to W
  loglik.W <- (t(ones) %*% sample$W %*% ones) * log(params$tau / (1 - params$tau)) + num.nodes * log(1 - params$tau)

  if (length(unique(class.labels)) == 1) {
    # Case of one class
    # params$P.ij is a scalar

    # NOTE: Factor of 1/2 is due to the fact that Y is undirected, thus
    #       we are double counting links in the matrix product below

    # Part of log-likelihood due to Y
    loglik.Y <- 1/2 * (t(ones) %*% sample$Y %*% ones) * log(params$P.ij / (1 - params$P.ij)) + num.nodes * log(1 - params$P.ij)
  } else {
    # Case of multiple classes
    edge.prob.mtx <- create.edge.prob.mtx(class.labels, params$P.ij)

    # Row and column indices for each entry in a num.nodes x num.nodes matrix
    # HACK: These indices are generated for using mapply() below
    row.indices <- rep(1:num.nodes, each = num.nodes)
    col.indices <- rep(1:num.nodes, times = num.nodes)

    # For each link in the sociomatrix, compute its contribution to the log-likelihood
    # Then take product of all of these contributions
    # NOTE: We only compute contributions for the upper triangle of the sociomatrix
    #       because the sociomatrix Y is assumed to be undirected.
    #       If the column index > row index, set the contribution = 1

    # Part of log-likelihood due to Y
    loglik.Y <- prod(mapply(function(i, j) if (j > i) (sample$Y[i, j] * log(edge.prob.mtx[i, j] / (1 - edge.prob.mtx[i, j])) + log(1 - edge.prob.mtx[i, j])) else 1, row.indices, col.indices))
  }

  # Add up contributions from all variables
  loglik <- loglik.Z0 + loglik.W + loglik.Y

  return(loglik)
}


llratio.toggle <- function(proposal, prev.sample, toggled, obs, params, class.labels) {
  # Compute likelihood ratio after a toggle
  #
  # Args:
  #   [proposal is proposed MCMC sample, a list of network variables]
  #   proposal$Y - network sample (matrix)
  #   proposal$Z0 - initial infection sample (vector)
  #   proposal$W - transmissibility matrix sample (matrix)
  #
  #   [prev.sample is previous sample in MCMC, a list of network variables]
  #   prev.sample$Y - network sample (matrix)
  #   prev.sample$Z0 - initial infection sample (vector)
  #   prev.sample$W - transmissibility matrix sample (matrix)
  #
  #   [toggled is list of information about toggled variable]
  #   toggled$var - MCMC sample variable that was toggled:"Y", "Z0", or "W"
  #   toggled$idx - index of entry toggled (scalar for Z0, array index for Y and W)
  #
  #   [obs is a list of observations]
  #   obs$Y.obs - observed network, derived from sample
  #   obs$Z.obs - observed infection, derived from sample
  #   obs$S - contact tracing sample as vector
  #   obs$ct.design - contact tracing design used for sample: "infected_only",
  #                   "infected_and_edge_units", "contacts_of_edge_units", "full_contact_components"
  #
  #   [params is a list of model parameters]
  #   params$P.ij - square, symmetric matrix containing probabilities; each entry
  #                  is the probability of a link between node class i and node class j
  #   params$eta - Bernoulli process parameter for initial infection
  #   params$tau - Bernoulli process parameter for generating transmission matrix
  #
  #   class.labels - class labels (string) for each node, as vector of strings;
  #                  each label must be in 1:num.nodes
  #
  # Returns:
  #   ratio - likelihood ratio of toggled proposal compared to previous sample

  # Initialize variable for whether reachability should be checked
  check.reach <- FALSE

  # NOTE: Krista's thesis suggests that when Y or W are toggled, reachability
  #       should only be checked if a transmission pathway is removed, but not
  #       if a transmission pathway is added.  However, this suggestion is not
  #       followed in Krista's code.  I decided to follow the example of
  #       Krista's code

  # Determine whether reachability needs to be checked
  if (toggled$var == "Z0") {
    # Z0 was toggled
    check.reach <- TRUE
  } else if (toggled$var == "Y") {
    # Y was toggled for an edge
    edge <- prev.sample$W[toggled$idx]
    # Compute reverse edge (flip row and column index)
    rev.edge <- prev.sample$W[matrix(rev(toggled$idx), nrow = 1)]
    # Check if W = 1 in either direction for link in Y
    if ((edge == 1) | (rev.edge == 1)) {
      # W = 1 for that corresponding link in either direction
      check.reach <- TRUE
    }
  } else if (toggled$var == "W") {
    # W was toggled for a directed edge
    if (prev.sample$Y[toggled$idx] == 1) {
      # Y = 1 for corresponding edge
      # NOTE: Don't need to check reverse edge, since Y should be symmetric
      check.reach <- TRUE
    }
  } else {
    stop("Error: toggle$var set to invalid value")
  }

  # Set Boolean variable consistent
  # This variable says whether the proposal is consistent with observation

  #################################################################################
  # NOTE: VERY IMPORTANT!  We only check that the sample is consistent with Z.obs #
  #       We assume that consistency of the sample with Y.obs is enforced in the  #
  #       in the toggle function toggle.Y                                         #
  #################################################################################

  consistent <- TRUE
  if (check.reach) {
    # Compute reachability matrix
    R.D <- reachability(proposal$W * proposal$Y)
    # Compute resulting infection
    Z <- as.vector(proposal$Z0 %*% R.D)
    Z[Z > 1] <- 1

    # Check if infection derived from proposal is consistent with observation
    if (! all(Z[obs$S == 1] == obs$Z.obs[obs$S == 1])) {
      consistent <- FALSE
      print("Z in MCMC sample inconsistent with Z.obs")
    }
  }

  # Compute design matrix
  design.mtx <- compute.design.mtx(obs$Z.obs, obs$S, obs$ct.design)

  if (! all(proposal$Y[design.mtx == 1] == obs$Y.obs[design.mtx == 1])) {
    consistent <- FALSE
    print("Y in MCMC sample inconsistent with Y.obs")
  }

  if (consistent) {
    # If proposal is reachable, then likelihood ratio is nonzero, and we compute it

    # NOTE: We only need to compute the part of the likelihood ratio
    #       corresponding to the variable that was toggled because only
    #       one variable is toggled at a time

    # Compute number of nodes
    num.nodes <- length(obs$S)
    # Construct a vector of ones for convenience
    ones <- rep(1, num.nodes)

    # Determine which variable was toggled and compute corresponding ratio
    if (toggled$var == "Z0") {
      # Change in number of infected nodes
      delta.Z0 <- t(proposal$Z0 - prev.sample$Z0) %*% ones
      # Factor in ratio due to Z0 toggle
      ratio.Z0 <- (params$eta / (1 - params$eta)) ^ delta.Z0

      ratio <- ratio.Z0
    } else if (toggled$var == "W") {
      # Change in number of transmissible links
      delta.W <- t(ones) %*% (proposal$W - prev.sample$W) %*% ones
      # Factor in ratio due to W toggle
      ratio.W <- (params$tau / (1 - params$tau)) ^ delta.W

      ratio <- ratio.W
    } else if (toggled$var == "Y") {
      if (length(unique(class.labels)) == 1) {
        # Case of one class
        # params$P.ij is a scalar

        # Change in number of sociomatrix links
        # NOTE: Factor of 1/2 is due to the fact that Y is undirected, thus
        #       we are double counting links in the matrix product below
        delta.Y <- 1/2 * (t(ones) %*% (proposal$Y - prev.sample$Y) %*% ones)
        ratio.Y <- (params$P.ij / (1 - params$P.ij)) ^ delta.Y
      } else {
        edge.prob.mtx <- create.edge.prob.mtx(class.labels, params$P.ij)

        # Row and column indices for each entry in a num.nodes x num.nodes matrix
        # HACK: These indices are generated for using mapply() below
        row.indices <- rep(1:num.nodes, each = num.nodes)
        col.indices <- rep(1:num.nodes, times = num.nodes)

        # For each link in the sociomatrix, compute its contribution to the likelihood ratio
        # Then take product of all of these contributions
        # NOTE: We only compute contributions for the upper triangle of the sociomatrix
        #       because the sociomatrix Y is assumed to be undirected.
        #       If the column index > row index, set the contribution = 1
        ratio.Y <- prod(mapply(function(i, j) if (j > i) (edge.prob.mtx[i, j] / (1 - edge.prob.mtx[i, j])) ^ (proposal$Y[i, j] - prev.sample$Y[i, j]) else 1, row.indices, col.indices))
      }
      ratio <- ratio.Y
    }
  } else {
    # Since proposal is inconsistent with observation, likelihood ratio must be 0
    ratio <- 0
  }

  return(ratio)
}

############### TOGGLE FUNCTIONS ###############

toggle.binary <- function(binary.int) {
  # Toggle a binary integer, such that 0 -> 1 and 1 -> 0
  #
  # Args:
  #   binary.int - binary integer that is 0 or 1
  #
  # Returns:
  #   toggled.val - toggled binary integer

  # Flip binary integer 0 -> 1 or 1 -> 0
  toggled.val <- -1 * binary.int + 1

  return(toggled.val)
}

toggle.Z0 <- function(prev.sample, obs) {
  # Toggle an initial infected node, one element of Z0
  #
  # Args:
  #   [prev.sample is previous sample in MCMC, a list of network variables]
  #   prev.sample$Y - network sample (matrix)
  #   prev.sample$Z0 - initial infection sample (vector)
  #   prev.sample$W - transmissibility matrix sample (matrix)
  #
  #   [obs is a list of observations]
  #   obs$Y.obs - observed network, derived from sample
  #   obs$Z.obs - observed infection, derived from sample
  #   obs$S - contact tracing sample as vector
  #   obs$ct.design - contact tracing design used for sample: "infected_only",
  #                   "infected_and_edge_units", "contacts_of_edge_units", "full_contact_components"
  #
  # Returns:
  #   [proposal is a list containing proposed sample in Markov chain]
  #   proposal$Y - network sample (matrix)
  #   proposal$Z0 - initial infection sample (vector)
  #   proposal$W - transmissibility matrix sample (matrix)
  #
  #   [toggled is list of information about toggled variable]
  #   toggled$var - MCMC sample variable that was toggled:"Y", "Z0", or "W"
  #   toggled$idx - index of Z0 entry toggled

  # Copy previous sample into proposal
  proposal <- prev.sample

  # Sample one of the nodes to get index
  # Don't toggle observed uninfected nodes
  # NOTE: Krista's code doesn't seem to account for observed uninfected nodes
  idx <- sample(which(! ((obs$Z0 = 0) & (obs$S == 1))), 1)

  # Toggle Z0 value for selected index
  proposal$Z0[idx] <- toggle.binary(prev.sample$Z0[idx])
  # Save toggle info
  toggled <- list(var = "Z0", idx = idx)

  return(list(proposal = proposal, toggled = toggled))
}

toggle.Y <- function(prev.sample, obs) {
  # Toggle a link in the sociomatrix Y
  #
  # NOTE: The sociomatrix Y is assumed to be undirected, so the toggle
  #       changes two entries in Y, i.e. if Y[i,j] is changed, so is Y[j, i]
  #
  # Args:
  #   [prev.sample is previous sample in MCMC, a list of network variables]
  #   prev.sample$Y - network sample (matrix)
  #   prev.sample$Z0 - initial infection sample (vector)
  #   prev.sample$W - transmissibility matrix sample (matrix)
  #
  #   [obs is a list of observations]
  #   obs$Y.obs - observed network, derived from sample
  #   obs$Z.obs - observed infection, derived from sample
  #   obs$S - contact tracing sample as vector
  #   obs$ct.design - contact tracing design used for sample: "infected_only",
  #                   "infected_and_edge_units", "contacts_of_edge_units", "full_contact_components"
  #
  # Returns:
  #   [proposal is a list containing proposed sample in Markov chain]
  #   proposal$Y - network sample (matrix)
  #   proposal$Z0 - initial infection sample (vector)
  #   proposal$W - transmissibility matrix sample (matrix)
  #
  #   [toggled is list of information about toggled variable]
  #   toggled$var - MCMC sample variable that was toggled:"Y", "Z0", or "W"
  #   toggled$idx - array index of Y matrix element toggled

  # Copy previous sample into proposal
  proposal <- prev.sample

  # Sample one of the links

  # Compute design matrix
  design.mtx <- compute.design.mtx(obs$Z.obs, obs$S, obs$ct.design)
  # Adjacency matrix whose non-zero entries corresponding to sociomatrix links which cannot be toggled
  # Don't sample/toggle observed links (given by design matrix)
  invalid.mtx <- design.mtx
  # Don't sample self links
  diag(invalid.mtx) <- 1
  # Array indices corresponding to sociomatrix links which can be toggled
  array.indices <- which(invalid.mtx == 0, arr.ind = TRUE)

  # Number of links that can be sampled (multipled by factor of 2 since Y is undirected)
  num.links.sample <- dim(array.indices)[1]

  # Error message if no links can be sampled
  if (num.links.sample == 0) {
    stop("Error: no sociomatrix links in Y can be toggled")
  }

  # Perform sample
  edge <- matrix(array.indices[sample(1:num.links.sample, 1), ], nrow = 1)
  # Toggle link in Y for selected index

  # Compute reverse directed edge
  rev.edge <- matrix(rev(edge), nrow = 1)
  # Toggle edges in both directions
  proposal$Y[edge] <- toggle.binary(prev.sample$Y[edge])
  proposal$Y[rev.edge] <- toggle.binary(prev.sample$Y[rev.edge])
  # Save toggle info
  toggled <- list(var = "Y", idx = edge)

  return(list(proposal = proposal, toggled = toggled))
}

# FIXME: Need to check links from unsampled nodes to infected subgraphs in infected only design

toggle.W <- function(prev.sample, obs) {
  # Toggle a link in the sociomatrix W
  #
  # Args:
  #   [prev.sample is previous sample in MCMC, a list of network variables]
  #   prev.sample$Y - network sample (matrix)
  #   prev.sample$Z0 - initial infection sample (vector)
  #   prev.sample$W - transmissibility matrix sample (matrix)
  #
  #   [obs is a list of observations]
  #   obs$Y.obs - observed network, derived from sample
  #   obs$Z.obs - observed infection, derived from sample
  #   obs$S - contact tracing sample as vector
  #   obs$ct.design - contact tracing design used for sample: "infected_only",
  #                   "infected_and_edge_units", "contacts_of_edge_units", "full_contact_components"
  #
  # Returns:
  #   [sample is a list containing next sample in Markov chain]
  #   sample$Y - network sample (matrix)
  #   sample$Z0 - initial infection sample (vector)
  #   sample$W - transmissibility matrix sample (matrix)
  #
  #   accept - Boolean variable indicating whether proposal was accepted
  #
  #   [toggled is list of information about toggled variable]
  #   toggled$var - MCMC sample variable that was toggled:"Y", "Z0", or "W"
  #   toggled$idx - array index of W matrix element toggled

  # Compute number of nodes in network
  num.nodes <- length(prev.sample$Z0)

  # Copy previous sample into proposal
  proposal <- prev.sample

  # Sample one of the links

  # Adjacency matrix whose non-zero entries corresponding to sociomatrix links which cannot be toggled
  invalid.mtx <- matrix(0, nrow = num.nodes, ncol = num.nodes)
  # Don't sample self links
  diag(invalid.mtx) <- 1
  # Array indices corresponding to sociomatrix links which can be toggled
  array.indices <- which(invalid.mtx == 0, arr.ind = TRUE)

  # Number of links that can be sampled (multipled by factor of 2 since Y is undirected)
  num.links.sample <- dim(array.indices)[1]
  # Perform sample
  edge <- matrix(array.indices[sample(1:num.links.sample, 1), ], nrow = 1)

  # Toggle link in W for selected index
  proposal$W[edge] <- toggle.binary(prev.sample$W[edge])
  # Save toggle info
  toggled <- list(var = "W", idx = edge)

  return(list(proposal = proposal, toggled = toggled))
}

############### MCMC FUNCTIONS ###############

propose.sample <- function(toggle.var, prev.sample, obs, params, class.labels) {
  # Propose next sample in Markov Chain
  #
  # Args:
  #   toggle.var - variable to toggle, either "Y", "Z0", or "W"
  #
  #   [prev.sample is previous sample in MCMC, a list of network variables]
  #   prev.sample$Y - network sample (matrix)
  #   prev.sample$Z0 - initial infection sample (vector)
  #   prev.sample$W - transmissibility matrix sample (matrix)
  #
  #   [obs is a list of observations]
  #   obs$Y.obs - observed network, derived from sample
  #   obs$Z.obs - observed infection, derived from sample
  #   obs$S - contact tracing sample as vector
  #   obs$ct.design - contact tracing design used for sample: "infected_only",
  #                   "infected_and_edge_units", "contacts_of_edge_units", "full_contact_components"
  #
  #   [params is a list of model parameters]
  #   params$P.ij - square, symmetric matrix containing probabilities; each entry
  #                  is the probability of a link between node class i and node class j
  #   params$eta - Bernoulli process parameter for initial infection
  #   params$tau - Bernoulli process parameter for generating transmission matrix
  #
  #   class.labels - class labels (string) for each node, as vector of strings;
  #                  each label must be in 1:num.nodes
  #
  # Returns:
  #   [proposal is a list containing proposed sample in Markov chain]
  #   proposal$Y - network sample (matrix)
  #   proposal$Z0 - initial infection sample (vector)
  #   proposal$W - transmissibility matrix sample (matrix)
  #
  #   ll.ratio - likelihood ratio of proposal compared to previous sample
  #
  #   accept - Boolean variable indicating whether proposal was accepted
  #            NOTE: value is 0 or 1, binary integer
  #
  #   [toggled is list of information about toggled variable in proposal]
  #   toggled$var - MCMC sample variable that was toggled:"Y", "Z0", or "W"
  #   toggled$idx - index of entry toggled (scalar for Z0, array index for Y and W)

  # Toggle sample according to specified toggle.var
  if (toggle.var == "Y") {
    result <- toggle.Y(prev.sample, obs)
  } else if (toggle.var == "Z0") {
    result <- toggle.Z0(prev.sample, obs)
  } else if (toggle.var == "W") {
    result <- toggle.W(prev.sample, obs)
  } else {
    stop("Error: Invalid value specified for toggle.var")
  }

  # Likelihood ratio of toggled sample compared to untoggled sample
  ratio <- llratio.toggle(result$proposal, prev.sample, result$toggled, obs, params, class.labels)
  # Boolean variable for whether or not proposal was accepted
  accept <- (ratio > runif(1))

  # NOTE: The Boolean accept variable is coerced to a binary integer
  return(list(proposal = result$proposal, ll.ratio = ratio, accept = as.integer(accept), toggled = result$toggled))
}

construct.prop.data <- function(prop.result) {
  # Construct single data case for a proposal result, for assignment in data frame
  #
  # Args:
  #
  #   prop.result is an object list(proposal = proposal, ll.ratio = ll.ratio, accept = accept, toggled = toggled)
  #
  #   [proposal is a list containing next sample in Markov chain]
  #   proposal$Y - network sample (matrix)
  #   proposal$Z0 - initial infection sample (vector)
  #   proposal$W - transmissibility matrix sample (matrix)
  #
  #   ll.ratio - likelihood ratio of sample compared to previous sample
  #
  #   accept - Boolean variable indicating whether proposal was accepted
  #            NOTE: value is 0 or 1, binary integer
  #
  #   [toggled is list of information about toggled variable]
  #   toggled$var - MCMC sample variable that was toggled:"Y", "Z0", or "W"
  #   toggled$idx - index of entry toggled (scalar for Z0, array index for Y and W)
  #
  # Returns:
  #   data.case - list of proposal data, to be recorded to data frame
  #
  #   The list consists of the following (in this order)
  #     accept - whether sample was accepted
  #     toggled.var - variable toggled in sample (if applicable)
  #     node - node toggled (if toggled.var = Z0)
  #     head - head of edge toggled (if toggled.var = Y or W)
  #     tail - tail of edge toggled (if toggled.var = Y or W)
  #     ll.ratio - likelihood ratio of sample compared to previous sample

  # Copy info from sample.result
  accept <- as.integer(prop.result$accept)
  ll.ratio <- prop.result$ll.ratio
  toggled <- prop.result$toggled

  toggled.var <- toggled$var

  # Assign node, head, and tail variables

  # node is only valid variable for toggled.var = Z0
  node <- NA
  # head and tail are only valid variables for toggled.var = Y or W
  head <- NA
  tail <- NA
  if (toggled.var == "Z0") {
    # node is Z0 index toggled
    node <- as.integer(toggled$idx)
  }
  if ((toggled.var == "Y") | (toggled.var == "W")) {
    edge <- toggled$idx
    # head of edge toggled
    head <- as.integer(edge[1, 1])
    # tail of edge toggled
    tail <- as.integer(edge[1, 2])
  }

  # Coerce character to factor
  toggled.var <- factor(toggled.var, levels = c("Z0", "Y", "W"))

  # Construct data case
  data.case <- list(accept, toggled.var, node, head, tail, ll.ratio)

  return(data.case)
}

compute.sample.stats <- function(sample, params, class.labels) {
  # Compute statistics for sample in Markov Chain and construct data case for stats data frame
  #
  # NOTE: For calculation of log-likelihood, assumes that sample is consistent with observations!
  #
  # Args:
  #
  #   [sample is a MCMC sample, a list of network variables]
  #   sample$Y - network sample (matrix)
  #   sample$Z0 - initial infection sample (vector)
  #   sample$W - transmissibility matrix sample (matrix)
  #
  #   [params is a list of model parameters]
  #   params$P.ij - square, symmetric matrix containing probabilities; each entry
  #                  is the probability of a link between node class i and node class j
  #   params$eta - Bernoulli process parameter for initial infection
  #   params$tau - Bernoulli process parameter for generating transmission matrix
  #
  #   class.labels - class labels (string) for each node, as vector of strings;
  #                  each label must be in 1:num.nodes
  #
  # Returns:
  #   data.case - list of MCMC sample statistics, to be recorded to data frame
  #
  #   The list consists of the following (in this order)
  #     loglik - log-likelihood of sample
  #     num.initial.infected - number of initial infected nodes in sample
  #     num.infected - number of infected nodes in sample
  #     num.Y.links - number of links in sociomatrix for sample
  #     num.W.links - number of links in transmissibility matrix for sample

  # Compute log-likelihood of sample
  # NOTE: This assumes that the sample is consistent with observations!
  loglik <- log.likelihood(sample, params, class.labels)

  # Compute reachability matrix
  R.D <- reachability(sample$W * sample$Y)
  # Compute resulting infection
  Z <- as.vector(sample$Z0 %*% R.D)
  Z[Z > 1] <- 1

  # Compute stats
  num.initial.infected <- as.integer(sum(sample$Z0))
  num.infected <- as.integer(sum(Z))
  num.Y.links <- as.integer(sum(sample$Y) / 2)
  num.W.links <- as.integer(sum(sample$W))

  # Construct data case as list of stats
  data.case <- list(loglik, num.initial.infected, num.infected, num.Y.links, num.W.links)

  return(data.case)
}

simulate.mcmc.chain <- function(obs, params, class.labels, options) {
  # Simulate Markov chain using Monte Carlo
  #
  # Args:
  #
  #   [obs is a list of observations]
  #   obs$Y.obs - observed network, derived from sample
  #   obs$Z.obs - observed infection, derived from sample
  #   obs$S - contact tracing sample as vector
  #   obs$ct.design - contact tracing design used for sample: "infected_only",
  #                   "infected_and_edge_units", "contacts_of_edge_units", "full_contact_components"
  #
  #   [params is a list of model parameters]
  #   params$P.ij - square, symmetric matrix containing probabilities; each entry
  #                  is the probability of a link between node class i and node class j
  #   params$eta - Bernoulli process parameter for initial infection
  #   params$tau - Bernoulli process parameter for generating transmission matrix
  #
  #   class.labels - class labels (string) for each node, as vector of strings;
  #                  each label must be in 1:num.nodes
  #
  #   [optional arguments]
  #
  #   options$initial.max.tries -
  #   options$num.Z0.toggles - number of Z0 toggles to attempt
  #   options$mult.Y.toggles - number of Y toggles to attempt per Z0 toggle
  #   options$mult.W.toggles - number of W toggles to attempt per Y toggle
  #   options$burnin - size of burnin for Markov chain
  #   options$interval - interval to sample after burnin is achieved
  #   options$collect.data - Boolean variable whether to collect data / stats about chain
  #
  # Returns:
  #   samples - list of Markov chain samples where each sample has the structure
  #
  #   [each sample is a list containing the following]
  #   sample$Y - network sample (matrix)
  #   sample$Z0 - initial infection sample (vector)
  #   sample$W - transmissibility matrix sample (matrix)
  #
  #   [[optional data returned, if collect.data = TRUE]]
  #
  #   proposals - list of proposals from Metropolis algorithm
  #
  #   [each proposal is a list containing the following]
  #   proposal$Y - network sample (matrix)
  #   proposal$Z0 - initial infection sample (vector)
  #   proposal$W - transmissibility matrix sample (matrix)
  #
  #   sample.stats
  #   prop.data

  # FIXME: Finish comments above

  # Compute number of nodes in socionetwork
  num.nodes <- length(obs$S)

  # Set defaults for options
  opts <- list(initial.max.tries = 100, num.Z0.toggles = 2, mult.Y.toggles = num.nodes, mult.W.toggles = 2, burnin = 0, collect.data = TRUE)
  # Fill in user-specified options
  opts[names(options)] <- options

  # Compute total samples in MCMC chain
  total.samples <- opts$num.Z0.toggles * opts$mult.Y.toggles * opts$mult.W.toggles + 1

  # --------------- Draw initial sample --------------- #

  success <- FALSE
  tries <- 0
  while ((! success) & (tries < opts$initial.max.tries)) {
    initial.result <- draw.initial.sample(obs, class.labels, params)
    tries <- tries + 1
    success <- initial.result$success
  }

  if (tries >= opts$initial.max.tries) {
    stop("Error: Failed to draw valid initial MCMC sample")
  }

  print("Succeeded in drawing initial sample")

  # Save initial sample
  initial.sample <- initial.result$sample

  # ---------- Initialize lists for samples, proposals, and stats ---------- #

  # Initialize list of MCMC samples
  samples <- replicate(total.samples, list())
  samples[[1]] <- initial.sample

  # Initialize optional data
  if (opts$collect.data) {
    # Initialize list of proposal samples
    proposals <- replicate(total.samples, list())
    proposals[[1]] <- initial.sample

    # Initialize sample stats data frame
    sample.stats <- data.frame(loglik = numeric(total.samples), num.initial.infected = integer(total.samples), num.infected = integer(total.samples), num.Y.links = integer(total.samples), num.W.links = integer(total.samples))
    sample.stats[1, ] <- rep(NA, length(names(sample.stats)))

    # QUESTION: Implement toggled.var in data frame as character or factor?
    #
    # Alternative implementation of toggled.var in data frame, as character instead of factor
    #
    # prop.data <- data.frame(accept = integer(total.samples), toggled.var = character(total.samples), node = integer(total.samples), head = integer(total.samples), tail = integer(total.samples), ll.ratio = numeric(total.samples), stringsAsFactors = FALSE)

    # Initialize proposal data frame
    prop.data <- data.frame(accept = integer(total.samples), toggled.var = factor(rep(NA, total.samples), levels = c("Z0", "Y", "W")), node = integer(total.samples), head = integer(total.samples), tail = integer(total.samples), ll.ratio = numeric(total.samples), stringsAsFactors = TRUE)
    prop.data[1, ] <- rep(NA, length(names(prop.data)))
  }

  # ----------------- Start sampling --------------- #

  # Save previous sample
  prev.sample <- initial.sample

  # Counter for sample number
  count <- 1

  for (i in 1:opts$num.Z0.toggles) {
    # ---------------- Toggle Z0 --------------- #
    print("Toggle Z0")
    prop.result <- propose.sample("Z0", prev.sample, obs, params, class.labels)

    # ----- Tedious boilerplate to update data ----- #

    count <- count + 1
    print(count)

    accept <- prop.result$accept

    if (accept == 1) {
      samples[[count]] <- prop.result$proposal
      prev.sample <- prop.result$proposal
    } else {
      samples[[count]] <- samples[[count - 1]]
    }

    if (opts$collect.data) {
      proposals[[count]] <- prop.result$proposal

      prop.data[count, ] <- construct.prop.data(prop.result)

      if (prop.result$accept == 1) {
        # Proposal accepted, so compute stats on proposal
        sample.stats[count, ] <- compute.sample.stats(prop.result$proposal, params, class.labels)
      } else {
        if (count == 2) {
          # Proposal not accepted, need to compute stats based on initial sample
          sample.stats[count, ] <- compute.sample.stats(initial.sample, params, class.labels)
        } else {
          # Proposal not accepted, copy stats from previous sample
          sample.stats[count, ] <- sample.stats[count - 1, ]
        }
      }
    }
    # ----- End of tedious boilerplate to update data ----- #

    for (j in 1:opts$mult.Y.toggles) {
      # ---------------- Toggle Y --------------- #
      print("Toggle Y")
      prop.result <- propose.sample("Y", prev.sample, obs, params, class.labels)

      # ----- Tedious boilerplate to update data ----- #

      count <- count + 1
      print(count)

      accept <- prop.result$accept

      if (accept == 1) {
        samples[[count]] <- prop.result$proposal
        prev.sample <- prop.result$proposal
      } else {
        samples[[count]] <- samples[[count - 1]]
      }

      if (opts$collect.data) {
        proposals[[count]] <- prop.result$proposal

        prop.data[count, ] <- construct.prop.data(prop.result)

        if (prop.result$accept == 1) {
          # Proposal accepted, so compute stats on proposal
          sample.stats[count, ] <- compute.sample.stats(prop.result$proposal, params, class.labels)
        } else {
          if (count == 2) {
            # Proposal not accepted, need to compute stats based on initial sample
            sample.stats[count, ] <- compute.sample.stats(initial.sample, params, class.labels)
          } else {
            # Proposal not accepted, copy stats from previous sample
            sample.stats[count, ] <- sample.stats[count - 1, ]
          }
        }
      }
      # ----- End of tedious boilerplate to update data ----- #

      for (k in 1:opts$mult.W.toggles) {
        # ---------------- Toggle W --------------- #
        print("Toggle W")
        prop.result <- propose.sample("W", prev.sample, obs, params, class.labels)

        # ----- Tedious boilerplate to update data ----- #

        count <- count + 1
        print(count)

        accept <- prop.result$accept

        if (accept == 1) {
          samples[[count]] <- prop.result$proposal
          prev.sample <- prop.result$proposal
        } else {
          samples[[count]] <- samples[[count - 1]]
        }

        if (opts$collect.data) {
          proposals[[count]] <- prop.result$proposal

          prop.data[count, ] <- construct.prop.data(prop.result)

          if (prop.result$accept == 1) {
            # Proposal accepted, so compute stats on proposal
            sample.stats[count, ] <- compute.sample.stats(prop.result$proposal, params, class.labels)
          } else {
            if (count == 2) {
              # Proposal not accepted, need to compute stats based on initial sample
              sample.stats[count, ] <- compute.sample.stats(initial.sample, params, class.labels)
            } else {
              # Proposal not accepted, copy stats from previous sample
              sample.stats[count, ] <- sample.stats[count - 1, ]
            }
          }
        }
        # ----- End of tedious boilerplate to update data ----- #

      }
    }
  }

  print("Initial sample")
  print(initial.sample)

  if (opts$collect.data) {
    return(list(samples = samples, proposals = proposals, sample.stats = sample.stats, prop.data = prop.data))
  } else {
    return(list(samples = samples))
  }

}

convert.samples.coda <- function(samples) {

  num.samples <- length(samples)
  Z0.data <- matrix(NA, nrow = num.samples, ncol = length(samples[[1]]$Z0))
  Y.data <- matrix(NA, nrow = num.samples, ncol = prod(dim(samples[[1]]$Y)))
  W.data <- matrix(NA, nrow = num.samples, ncol = prod(dim(samples[[1]]$W)))

  for (i in 1:num.samples) {
    Z0.data[i, ] <- samples[[i]]$Z0
    Y.data[i, ] <- as.vector(samples[[i]]$Y)
    W.data[i, ] <- as.vector(samples[[i]]$W)
  }

  return(list(Z0.mcmc = mcmc(Z0.data), Y.mcmc = mcmc(Y.data), W.mcmc = mcmc(W.data)))
}

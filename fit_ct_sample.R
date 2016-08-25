# Functions for fitting contact tracing samples

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ! NOTE: We don't estimate sigma, because we use missing-at-random assumption, !
# !     and therefore sampling parameters can be ignored                        !
# !                                                                             !
# ! NOTE: We assume that links can be observed perfectly (no degradation), i.e. !
# !       p.trace.infected = 1, p.trace.uninfected = 1                          !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Need for components functions, components() and component.dist()
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
  # NOTE: We don't estimate sigma, because we use missing-at-random assumption,
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

  num.nodes <- length(class.labels)
  num.classes <- length(unique(class.labels))

  # Compute design matrix
  design.mtx <- compute.design.mtx(obs$Z.obs, obs$S, obs$ct.design)

  # ---------- INITIAL ESTIMATE FOR P.ij ---------- #

  # NOTE: Factor of 1/2 in num.links.obs and num.links is due to the fact
  #       that Y is undirected

  # Simple calculation if only one class
  if (num.classes == 1) {
    # Construct a vector of ones for convenience
    ones <- rep(1, num.nodes)

    # Number of observed links
    num.links.obs <- 1 / 2 * (t(ones) %*% obs$Y.obs %*% ones)
    # Number of possible links
    num.links <- 1 / 2 * (t(ones) %*% design.mtx %*% ones)
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
          Y.obs.ii<- obs$Y.obs[class.labels == unique.labels[i]]
          # Extract sub-matrix of design matrix corresponding to class i
          design.mtx.ii <- design.mtx[class.labels == unique.labels[i]]

          # Number of observed links
          num.links.obs <- 1 / 2 * (t(ones.ii) %*% Y.obs.ii %*% ones.ii)
          # Number of possible links
          num.links <- 1 /2 * (t(ones.ii) %*% design.mtx.ii %*% ones.ii)

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
          Y.obs.ii <-  obs$Y.obs[class.labels == unique.labels[i]]
          # Extract sub-matrix of observed network corresponding to class j
          Y.obs.jj <-  obs$Y.obs[class.labels == unique.labels[j]]
          # Extract sub-matrix of observed network corresponding to classes i and j
          Y.obs.ij <- obs$Y.obs[class.labels == unique.labels[i] | class.labels == unique.labels[j]]

          # Extract sub-matrix of design matrix corresponding to class i
          design.mtx.ii <- design.mtx[class.labels == unique.labels[i]]
          # Extract sub-matrix of design matrix corresponding to class j
          design.mtx.jj <- design.mtx[class.labels == unique.labels[j]]
          # Extract sub-matrix of design matrix corresponding to classes i and j
          design.mtx.ij <- design.mtx[class.labels == unique.labels[i] | class.labels == unique.labels[j]]

          # Number of observed links
          num.links.obs <- (t(ones.ij) %*% Y.obs.ij %*% ones.ij) - (t(ones.ii) %*% Y.obs.ii %*% ones.ii) - (t(ones.jj) %*% Y.obs.jj %*% ones.jj)
          num.links.obs <- 1 / 2 * num.links.obs
          # Number of possible links
          num.links <- (t(ones.ij) %*% design.mtx.ij %*% ones.ij) - (t(ones.ii) %*% design.mtx.ii %*% ones.ii) - (t(ones.jj) %*% design.mtx.jj %*% ones.jj)
          num.links <- 1 / 2 * num.links
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
  eta0 <- Q / (t(S) %*% ones)

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

  # Error checking for tau
  if ((tau0 < 0) | (tau0 > 1)) {
    stop("Error calculating initial parameters: tau0 must be between 0 and 1")
  }

  # Return estimate of initial model parameters
  return(list(P.ij = P.ij0, eta = eta0, tau = tau0))
}

initial.mcmc.sample <- function(obs, class.labels, params) {
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
        directed.edges.list <- compute.dst.bf(root, Y.obs.infected[subgraph.infected.nodes, subgraph.infected.nodes])

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
  # Eliminate any links that are unobserved or which have already been set in W
  # HACK: R automatically casts (W == NA) from Boolean to integer matrix
  Y.obs.untransmitted <- Y.obs.untransmitted * obs$Y.obs * (W == NA)

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
  Z <- Z0 %*% R.D
  Z[Z > 1] <- 1

  # FIXME: Should get rid of this error message
  if (! all(Z[obs$S == 1] == obs$Z.obs[obs$S == 1])) {
    stop("Error calculating initial MCMC sample: Sample inconsistent with obs$Z.obs")
  }

  # Return sample for Y, Z0, and W
  return(list(Y = Y, Z0 = Z0, W = W))
}

############### TOGGLE FUNCTIONS ###############

toggle.binary.int <- function(binary.int) {
  return(as.integer(! binary.int))
}

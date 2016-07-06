# Functions to simulate an infected (undirected) network

# To simulate an infected network, do the following steps:
#   1. Generate an undirected network [generate.network()]
#   2. Spread the infection across the network [spread.infection()]

# TODO: Other types of networks?  Directed networks?

# Load network package without all the startup messages
suppressMessages(library(network))

# ---------------- FUNCTIONS ---------------------

create.edge.prob.mtx <- function(nodes.per.class, P.ij) {
  # Creates a matrix whose entries are probabilities that an edge will link 
  # each pair of nodes
  #
  # Args:
  #   nodes.per.class - vector containing the number of nodes in class 1,
  #                     the number of nodes in class 2, etc.
  #   P.ij - square, symmetric matrix containing probabilities; each entry
  #          is the probability of a link between node class i and node class j 
  #
  # Returns:
  #   edge.prob.mtx - Matrix whose entries are probabilities that an edge will
  #                   link each pair of nodes. The dimension of the matrix is
  #                   the total number of nodes over classes, i.e.
  #                   sum(nodes.per.class)

  # Number of classes
  num.classes <- length(nodes.per.class)
  # Number of nodes in network
  num.nodes <- sum(nodes.per.class)

  # Error checking: check that P.ij is square matrix
  if (length(P.ij) > 1) {
    if (dim(P.ij)[1] != dim(P.ij)[2]) {
      stop("P.ij not a square matrix")
    }
  }

  # Error checking: check that P.ij is symmetric matrix
  # Probability of link between class i and class j should be same as
  # probability of link between class j and class i
  if (! all(P.ij == t(P.ij))) {
    stop("P.ij is not symmetric")
  }

  # Error checking: check that nodes.per.class and P.ij are consistent
  if (length(P.ij) > 1) {
    if (num.classes != dim(P.ij)[1]) {
      stop("nodes.per.class and P.ij have inconsistent dimensions")
    }
  }

  # Case where there is only one node class
  if ((num.classes == 1) && (length(P.ij) == 1)) {
    edge.prob.mtx <- matrix(P.ij, nrow = nodes.per.class, ncol = nodes.per.class)
  # Case where there is more than one node class
  } else {
    # Vector where the value of entry i is the class for node i
    node.class.vec <- rep(1:num.classes, times = nodes.per.class)
    # Row and column indices for each entry in a num.nodes x num.nodes matrix
    # Hack for using mapply
    row.indices = rep(1:num.nodes, each = num.nodes)
    col.indices = rep(1:num.nodes, times = num.nodes)
    # For each pair of nodes i and j, compute probability of link between them
    # The probability is determined by the classes of node i and j
    # We use mapply to vectorize this calculation and avoid for loops
    edge.prob.mtx <- matrix(mapply(function(i,j) P.ij[node.class.vec[i], node.class.vec[j]], row.indices, col.indices), nrow = num.nodes, ncol = num.nodes)
  }
  # Set diagonal of matrix to zero (node can't be linked to itself)
  diag(edge.prob.mtx) <- 0

  return(edge.prob.mtx)
}

generate.network <- function(nodes.per.class, P.ij) {
  # Generates a network simulated according to nodes.per.class and the
  # probability matrix for links between nodes of different classes P.ij
  #
  # Args:
  #   nodes.per.class - vector containing the number of nodes in class 1,
  #                     the number of nodes in class 2, etc.
  #   P.ij - square, symmetric matrix containing probabilities; each entry
  #          is the probability of a link between node class i and node class j 
  #
  # Returns:
  #   socio.net - network object containing the simulated network

  # Number of nodes in network
  num.nodes <- sum(nodes.per.class)
  # Create matrix containing prob each node is linked to another node
  prob.mtx <- create.edge.prob.mtx(nodes.per.class, P.ij)
  # Set lower triangle to zero, so we only generate a link between i and j once
  prob.mtx[lower.tri(prob.mtx)] <- 0
  # Create network
  # Only generate links for upper triangle to ensure consistency
  # We only want to generate a link between i and j once
  upper.socio.mtx <- matrix(rbinom(num.nodes^2, 1, prob.mtx), nrow = num.nodes, ncol = num.nodes)
  # After network has been generated, fill in lower triangle
  full.socio.mtx <- upper.socio.mtx + t(upper.socio.mtx)
  socio.net <- network(full.socio.mtx, directed = FALSE, hyper = FALSE, loops = FALSE, multiple = FALSE)

  if (! all(full.socio.mtx == as.sociomatrix(socio.net))) {
    stop("Network initialization inconsistent with input data")
  }

  return(socio.net)
}

spread.infection <- function(socio.net, eta, tau, verbose = FALSE, nonstochastic = FALSE, Z0 = NULL, W = NULL) {
  # Spread infection across given network and return infected network
  #
  # Args:
  #   socio.net - uninfected, undirected network (needs to be network object)
  #   eta - Bernoulli process parameter for initial infection 
  #   tau - Bernoulli process parameter for generating transmission matrix
  #   verbose - flag to print variables during infection process (for debugging)
  #   nonstochastic - flag to use pre-determined values for Z0 and W (for debugging)
  #   Z0 - (to be specified if nonstochastic = TRUE) vector of initial infected nodes
  #   W - (to be specified if nonstochastic = TRUE) infection trasmission matrix
  #
  # Returns:
  #   Z0 - vector of initial infected nodes (1 for infected, 0 for not infected) 
  #   Z - vector of final infected nodes (1 for infected, 0 for not infected)
  #   W.net - directed network object corresponding to infection transmission matrix
  #   infected.net - infected network (network object) with edge and vertex attributes describing infection

  # Error checking: make sure Z0 and W are specified if in nonstochastic mode
  if ((nonstochastic) & is.null(Z0) & is.null(W)) {
    stop("If nonstochastic flag is TRUE, need to specify values for Z0 and W")
  }

  # Number of nodes in network
  num.nodes <- network.size(socio.net)
  # Extract sociomatrix from network
  socio.mtx <- as.sociomatrix(socio.net)
  # Generate initial infection Z0 (1 for infected, 0 for not infected)
  # Multiple ways to do this: see tosource.R
  # Assume independent homogenous Bernoulli processes on all nodes with parameter eta
  # TODO: Implement other ways of generating initial infection [see tosource.r spread()]
  if (! nonstochastic) {
    Z0 <- rbinom(num.nodes, 1, eta)
  }
  # Current infected nodes
  Z <- Z0
  # Transmissibility matrix W (1 for edge that can transmit infection, 0 otherwise)
  # Assume independent homogenous Bernoulli processes on all nodes wih parameter tau
  if (! nonstochastic) {
    W <- matrix(rbinom(num.nodes^2, 1, tau), nrow = num.nodes, ncol = num.nodes)
  }
  diag(W) <- 0 # A node can't transmit infection to itself
  # Multiply element-wise by sociomatrix to get all possible transmissions in
  # this particular network instance
  # If we didn't do this multiplication, W would contain transmissions that
  # aren't possible in the network
  W <- W * socio.mtx
  W.net <- as.network(W, directed = TRUE, hyper = FALSE, loops = FALSE, multiple = FALSE)
  # Initialize spread.edges to zero
  # Elements of the vector which are 1 represented edges that spread infection
  spread.edges <- rep(0, network.edgecount(socio.net))

  if (verbose) {
    cat("Y:\n")
    print(socio.mtx)
    cat("\nZ0: ", Z0, "\n")
    cat("\nW:\n")
    print(W)
    cat("\n")
  }

  # Keep transmitting infection across network until no new infections
  repeat {
    # Perform next wave of infections
    Z.next <- Z %*% W
    Z.next <- as.vector(Z.next)
    Z.next[which(Z.next & Z)] <- 0 # Only keep new infections
    Z.next[which(Z.next > 0)] <- 1 # Truncate any positve values to 1

    if (verbose) {
      cat("Z.next: ", Z.next, "\n")
    }

    # Mark edges that spread infection in this wave
    for (j in which(Z.next == 1)) {
      # Nodes that could have spread the infection to node j
      # Multiply by Z, only previously infected nodes can spread infection
      in.nodes <- which(Z * W[ , j] > 0)
      for (k in in.nodes) {
        # Mark edge corresponding to nodes j and k as spreading infection
        spread.edges[get.edgeIDs(socio.net, j, alter = k)] <- 1

        if (verbose) {
          cat("Infection spread on edge", "(", j, ",", k, ") with ID ")
          print(get.edgeIDs(socio.net, j, alter = k))
        }
      }
    }

    # Update current infected nodes
    Z <- Z + Z.next

    if (verbose) {
      cat("Z: ", Z, "\n")
    }

    # If no new infections, break out of loop
    if (sum(Z.next) < 1) {
      break
    }
  }

  # Add infection attributes to network
  set.vertex.attribute(socio.net, "initial.infection", Z0)
  set.vertex.attribute(socio.net, "infected", Z)
  set.edge.attribute(socio.net, "spread", spread.edges)

  # Return data from infected network
  # FIXME: Seems redundant to return Z0 and Z if they are already attributes in infected.net 
  return(list(Z0 = Z0, Z = Z, W.net = W.net, infected.net = socio.net))
}

# Number of nodes in each class
#nodes.per.class <- c(25,5,25)
# Matrix containing probability of link between class i and class j
#P.ij <- matrix(c(0.05, 0.04, 0, 0.04, 0.08, 0.04, 0, 0.04, 0.05), nrow = 3, ncol = 3)
#tau = 0.5
#eta = 0.05
#zeta = 0.015

#edge.prob.mtx <- create.edge.prob.mtx(nodes.per.class, P.ij)

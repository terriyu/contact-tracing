# Functions to simulate a network

create.edge.prob.mtx <- function(nodes.per.class, P_ij) {
  # Number of classes
  num.classes <- length(nodes.per.class)
  # Number of nodes in network
  num.nodes <- sum(nodes.per.class)

  # Check that matrix P.ij is symmetric
  # Probability of link between class i and class j should be same as
  # probability of link between class j and class i
  if (! all(P.ij == t(P.ij))) {
    print("Error: P_ij is not symmetric")
    return()
  } 

  if (dim(P.ij)[1] != dim(P.ij)[2]) {
    print("Error: P_ij not a square matrix")
    return()
  }

  # Check that nodes.per.class and P.ij are consistent
  if (num.classes != dim(P.ij)[1]) {
    print("Error: nodes.per.class and P.ij have inconsistent dimensions")
    return()  
  }

  # Case where there is only one node class
  if ((num.classes == 1) && (length(P.ij) == 1)) {
    edge.prob.mtx <- matrix(P.ij, nrow = nodes.per.class, ncol = nodes.per.class)	
    # Set diagonal of matrix to zero (node can't be linked to itself)
    diag(edge.prob.mtx) <- 0
  } else {
    # Vector where the value of entry i is the class for node i
    node.class.vec <- rep(1:num.classes, times = nodes.per.class)
    # Row and column indices for each entry in a num.nodes x num.nodes matrix 
    # Hack for using mapply
    row.indices = rep(1:num.nodes, each = num.nodes)
    col.indices = rep(1:num.nodes, times = num.nodes)
    # For each pair of nodes i and j, compute probability of link between them
    # The probability is determined by the classes of node i and j
    edge.prob.mtx <- matrix(mapply(function(i,j) P.ij[node.class.vec[i], node.class.vec[j]], row.indices, col.indices), nrow = num.nodes, ncol = num.nodes) 
  } 

  return(edge.prob.mtx)
}

#generate.network <- function(nodes.per.class, P.ij, tau= 0.5, eta = 0.05, zeta = 0.015) {
#}


# Number of nodes in each class
nodes.per.class <- c(25,5,25)
# Matrix containing probability of link between class i and class j
P.ij <- matrix(c(0.05, 0.04, 0, 0.04, 0.08, 0.04, 0, 0.04, 0.05), nrow = 3, ncol = 3)
tau = 0.5
eta = 0.05
zeta = 0.015

edge.prob.mtx <- create.edge.prob.mtx(nodes.per.class, P_ij)


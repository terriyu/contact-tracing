# Script to simulate networks and contact tracing samples

source("simulate_network.R")
source("ct_sampling.R")

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !! More general formulation is multiple node classes, but here   !!
# !! we only consider one class, so nodes.per.class and P.ij       !!
# !! are scalars num.nodes and p, instead of a vector and a matrix !!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

###############################################################
# ------------------ PARAMETERS TO SPECIFY ------------------ #
###############################################################

# ---------- GENERAL SIMULATION PARAMETERS ---------- #

# Fix random seed
set.seed(1492)
#set.seed(137)

# File name prefix for writing simulation data
file.prefix <- "sim_data"
# If use.rds is set to TRUE, uses saveRDS() instead of save()
use.rds <- TRUE

# Number of networks and samples to simulate
num.sims <- 100
#num.sims <- 1

# ---------- MODEL PARAMETERS ---------- #

# Parameters for network, disease, and sampling models
num.nodes <- 200  # Number of nodes in network (nodes.per.class for multiple classes)
p <- 0.01         # Bernoulli parameter for link between two nodes (P.ij for multiple classes)
eta <- 0.065      # Bernoulli parameter for initial infection
tau <- 0.55       # Bernoulli parameter for transmission
sigma <- 0.11     # Bernoulli parameter for initial sample
#num.nodes <- 6
#p <- 0.5
#eta <- 0.2
#tau <- 0.2
#sigma <- 0.2

# Save model parameters in list
model.params <- list(nodes.per.class = num.nodes, P.ij = p, eta = eta, tau = tau, sigma = sigma)

# ---------- INFECTION DESIGN PARAMETERS ---------- #

# Fix number of nodes infected
initial.infect.method <- "nodes_fixed"
num.infect <- round(eta * num.nodes)

# Save infection parameters in list
infect.params <- list(initial.infect.method = initial.infect.method, num.infect = num.infect)

# ---------- CONTACT TRACING PARAMETERS ---------- #

ct.design <- "contacts_of_edge_units"
size.S0 <- round(sigma * num.nodes)
num.waves <- NULL
p.trace.infected <- NULL
p.trace.uninfected <- NULL

# Save contact tracing parameters in list
ct.params <- list(ct.design = ct.design, size.S0 = size.S0, num.waves = num.waves, p.trace.infected = p.trace.infected, p.trace.uninfected = p.trace.uninfected)

###############################################################
# ----------------------- SIMULATION ------------------------ #
###############################################################

# Initialize list of empty lists
sim.data.list <- replicate(num.sims, list())

for (i in 1:num.sims) {
  # Generate network and spread infection
  socio.net <- generate.network(num.nodes, p)
  spread.result <- spread.infection(socio.net, eta, tau, initial.infect.method, num.infect)
  # Unpack network and disease results
  Z0 <- spread.result$Z0                    # Initial infected nodes
  Z <- spread.result$Z                      # Infected nodes
  W.net <- spread.result$W.net              # Transmission network
  disease.net <- spread.result$disease.net  # Network with disease attributes

  # Perform contact tracing sample
  samp.result <- ct.sample(disease.net, sigma, ct.design, size.S0, num.waves, p.trace.infected, p.trace.uninfected)
  # Unpack sampling results
  S0 <- samp.result$S0  # Initial sampled nodes
  S <- samp.result$S    # Sampled nodes

  # Compute observed network and infection from network and contact tracing sample
  obs.result <- compute.observed(disease.net, S, ct.design)
  # Unpack observed results
  Y.obs <- obs.result$Y.obs  # Observed network
  Z.obs <- obs.result$Z.obs  # Observed infection

  # Save simulation data to list
  sim.data.list[[i]] <- list(disease.net = disease.net, Z0 = Z0, Z = Z, W.net = W.net, S0 = S0, S = S, Y.obs = Y.obs, Z.obs = Z.obs)
}

################################################################
# ------------------ WRITE RESULTS TO FILE ------------------- #
################################################################

# Write all parameters and simulation data to file
if (use.rds) {
  filename <- paste0(file.prefix, ".rds")
  saveRDS(list(model.params = model.params, infect.params = infect.params, ct.params = ct.params, sim.data.list = sim.data.list), file = filename)
} else {
  filename <- paste0(file.prefix, ".RData")
  save(model.params, infect.params, ct.params, sim.data.list, file = filename)
}

# Script to check and diagnose how the MCMC fit is working

source("simulate_network.R")
source("ct_sampling.R")
source("fit_ct_sample.R")

library(coda)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !! More general formulation is multiple node classes, but here   !!
# !! we only consider one class, so nodes.per.class and P.ij       !!
# !! are scalars num.nodes and p, instead of a vector and a matrix !!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

###############################################################
# ------------------ PARAMETERS TO SPECIFY ------------------ #
###############################################################

# ---------- GENERAL SIMULATION PARAMETERS ---------- #

# Number of networks and samples to simulate
#num.sims <- 100
num.sims <- 1

# ---------- MODEL PARAMETERS ---------- #

# Parameters for network, disease, and sampling models
#set.seed(1)
#num.nodes <- 50  # Number of nodes in network (nodes.per.class for multiple classes)
#p <- 0.01         # Bernoulli parameter for link between two nodes (P.ij for multiple classes)
#eta <- 0.065      # Bernoulli parameter for initial infection
#tau <- 0.55       # Bernoulli parameter for transmission
#sigma <- 0.11     # Bernoulli parameter for initial sample

set.seed(1492)
num.nodes <- 10
p <- 0.3
eta <- 0.5
tau <- 0.3
sigma <- 0.2

# ---------- INFECTION DESIGN PARAMETERS ---------- #

# Fix number of nodes infected
initial.infect.method <- "nodes_fixed"
num.infect <- round(eta * num.nodes)

# ---------- CONTACT TRACING PARAMETERS ---------- #

ct.design <- "contacts_of_edge_units"
size.S0 <- round(sigma * num.nodes)
num.waves <- NULL
p.trace.infected <- NULL
p.trace.uninfected <- NULL

###############################################################
# ----------------------- SIMULATION ------------------------ #
###############################################################

# Initialize list of empty lists
sim.data.list <- replicate(num.sims, list())

for (i in 1:num.sims) {
  # Generate network and spread infection
  class.labels <- rep(1, times = num.nodes)
  socio.net <- generate.network(class.labels, p)
  spread.result <- spread.infection(socio.net, eta, tau, initial.infect.method, options = list(num.infect = num.infect))
  # Unpack network and disease results
  Z0 <- spread.result$Z0                    # Initial infected nodes
  Z <- spread.result$Z                      # Infected nodes
  W.net <- spread.result$W.net              # Transmission network
  disease.net <- spread.result$disease.net  # Network with disease attributes

  # Perform contact tracing sample
  samp.result <- ct.sample(disease.net, sigma, ct.design, options = list(size.S0 = size.S0, num.waves = num.waves, p.trace.infected = p.trace.infected, p.trace.uninfected = p.trace.uninfected))
  # Unpack sampling results
  S0 <- samp.result$S0  # Initial sampled nodes
  S <- samp.result$S    # Sampled nodes

  # Compute observed network and infection from network and contact tracing sample
  obs.result <- compute.observed(disease.net, S, ct.design)
  # Unpack observed results
  Y.obs <- obs.result$Y.obs  # Observed network
  Z.obs <- obs.result$Z.obs  # Observed infection
  D <- obs.result$D # Design matrix

  # Save simulation data to list
  sim.data.list[[i]] <- list(disease.net = disease.net, Z0 = Z0, Z = Z, W.net = W.net, S0 = S0, S = S, Y.obs = Y.obs, Z.obs = Z.obs, D = D)
}

# Assuming one observation
obs <- list(Y.obs = Y.obs, Z.obs = Z.obs, S = S, ct.design = ct.design)

initial.params <- estimate.initial.params(obs, class.labels)
sim.result <- simulate.mcmc.chain(obs, initial.params, class.labels, options = list(num.Z0.toggles = 100, initial.max.tries = 1000, collect.data = TRUE))
samples <- sim.result$samples
proposals <- sim.result$proposals
prop.data <- sim.result$prop.data
sample.stats <- sim.result$sample.stats

# Plot network
# To view network plotting options, run "help(plot.network)"

# Vector of infection status for network nodes
infection.status <- disease.net %v% "infected"

# Vector of colors for infection status
# To view color choices, run "palette()"
infection.col <- infection.status
# Use red color for infected nodes
infection.col[infection.status == 1] <- which(palette() == "red")
# Use magenta color for observed infected nodes
infection.col[infection.status == 1 & S == 1] <- which(palette() == "magenta")
# Use blue color for uninfected nodes
infection.col[infection.status == 0] <- which(palette() == "blue")
# Use cyan color for observed uninfected nodes
infection.col[infection.status == 0 & S == 1] <- which(palette() == "cyan")

par(mfrow = c(1, 2))
socio.plot <- plot(disease.net, vertex.col = infection.col, displaylabels = TRUE, label = 1:num.nodes, main = "Sociomatrix")
plot(W.net, vertex.col = infection.col, displaylabels = TRUE, label = 1:num.nodes, main = "Transmissibility matrix", coord = socio.plot)

mcmc.result <- convert.samples.coda(samples)
Z0.mcmc <- mcmc.result$Z0.mcmc
Y.mcmc <- mcmc.result$Y.mcmc
W.mcmc <- mcmc.result$W.mcmc

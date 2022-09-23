using Printf
println("Starting compilation of libraries.")
using Turing
using DifferentialEquations

# Load StatsPlots for visualizations and diagnostics.
using StatsPlots

using LinearAlgebra

# Set a seed for reproducibility.
using Random


println("Setting the seed for Random to be reproducible.")
Random.seed!(14);

println("Defining LoktaVolterra Model.")
# Define Lotka-Volterra model.
function lotka_volterra(du, u, p, t)
    # Model parameters.
    α, β, γ, δ = p
    # Current state.
    x, y = u

    # Evaluate differential equations.
    du[1] = (α - β * y) * x # prey
    du[2] = (δ * x - γ) * y # predator

    return nothing
end

println("Defining initial-value problem.")
# Define initial-value problem.
u0 = [1.0, 1.0]
p = [1.5, 1.0, 3.0, 1.0]
tspan = (0.0, 10.0)
prob = ODEProblem(lotka_volterra, u0, tspan, p)

println("Solving problem.")
sol = solve(prob, Tsit5(); saveat=0.1)

println("Generating noisy data.")
odedata = Array(sol) + 0.8 * randn(size(Array(sol)))

println("Plotting simulation and noisy data.")
# Plot simulation and noisy observations.
plot(sol; alpha=0.3)
scatter!(sol.t, odedata'; color=[1 2], label="")
charts_file = "./charts/loktaVolterra.png"
savefig(charts_file)
println("The chart was saved to $charts_file")
# println("Fin")

# Direct handling of Bayesian estimation with Turing
println("Defining lokta-volterra Turing model.")
@model function fitlv(data, prob)
    # Prior distributions.
    σ ~ InverseGamma(2, 3)
    α ~ truncated(Normal(1.5, 0.5), 0.5, 2.5)
    β ~ truncated(Normal(1.2, 0.5), 0, 2)
    γ ~ truncated(Normal(3.0, 0.5), 1, 4)
    δ ~ truncated(Normal(1.0, 0.5), 0, 2)

    # Simulate Lotka-Volterra model. 
    p = [α, β, γ, δ]
    predicted = solve(prob, Tsit5(); p=p, saveat=0.1)

    # Observations.
    for i in 1:length(predicted)
        data[:, i] ~ MvNormal(predicted[i], σ^2 * I)
    end

    return nothing
end

model2 = fitlv(odedata, prob)
println("Sampling from chain.")
# Sample 3 independent chains with forward-mode automatic differentiation (the default).
chain = sample(model2, NUTS(0.65), MCMCSerial(), 1000, 3; progress=false)

println("Plotting chain samples.")
plot_chain = plot(chain)
plot_chain_path = "./charts/loktaVolterraChain.png"
savefig(plot_chain, plot_chain_path)
println("The chart was saved to $plot_chain_path")
println("Fin")

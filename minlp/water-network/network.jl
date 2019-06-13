#' # Optimizing water networks with Knitro
#' This notebook explains how to solve effectively NLP and MINLP
#' problem with the Knitro solver.
#'
#' The problem studied hereafter was originally designed by
#' Pierre Carpentier, for educational purposes. All credits
#'
#' We start by importing JuMP and Knitro

using JuMP, KNITRO

#' Import data.
include("data.jl");
# Plotting utilities
include("utils.jl");

#' ## First case study: the non-linear problem
#'
#' Define Knitro Optimizer with JuMP.
model = Model(with_optimizer(KNITRO.Optimizer, presolve=0, outlev=3, algorithm=4))

#' Write non-linear optimization problem first.
# Parameters.
α1 = Ar' * pr
# Dimension of the problem.
nx = n - md
@variable(model, qc[1:nx])
# Add dummy variable in the model.
@variable(model, q[1:n])
@constraint(model, q .== q0 + B*qc)
@NLobjective(model, Min,
             sum(r[i] * abs(q[i]) * q[i]^2 / 3 + α1[i] * q[i] for i in 1:n))
optimize!(model)

#' Display results
optimal_flow = JuMP.value.(q)
fig = figure()
plot_network(flow=optimal_flow)
display(fig)


#' ## Extension to Mixed-integer non-linear programming

#' We define hereafter the MINLP version of the problem.
# Add a maximum flows through the pipes
const QMAX = 10.

function load_mip_model!(model::JuMP.Model; nremovals=3)
    @variable(model, qc[1:nx])
    @variable(model, z[1:n], Bin)
    @variable(model, q[1:n])

    @constraint(model, q .== q0 + B*qc)
    @constraint(model, sum(z) == n - nremovals)
    @constraint(model,  q .<= QMAX * z)
    @constraint(model, -q .<= QMAX * z)

    @NLobjective(model, Min,
                sum(r[i] * abs(q[i]) * q[i]^2 / 3 + α1[i] * q[i] for i in 1:n))
    return
end

#' ### Default Knitro

#' Solve with default Knitro.
# Build non-linear solver.
model = Model(with_optimizer(KNITRO.Optimizer, outlev=3))
load_mip_model!(model)
@time JuMP.optimize!(model)

#' In the resolving, Knitro uses a Branch & Bound algorithm to find
#' the optimal solution corresponding to the (convex) MINLP problem.
#' The approach used is different than in BONMIN (outer approximation).

#' Observe that Knitro default computes 165 nodes before finding the solution.
#' Is it possible to find a better tuning for Knitro?

#' ### MIP-Tuner

#' Since Knitro 12.0, a MIP tuner was added to compute the optimal
#' parameterization to solve a given MIP problem with Knitro. The MIP
#' tuner uses an exhaustive search to find the optimal setting.
model = Model(with_optimizer(KNITRO.Optimizer, outlev=3, tuner=1))
load_mip_model!(model)
@time JuMP.optimize!(model)

# Plot!
optimal_flow = JuMP.value.(model[:q])
fig = figure()
plot_network(flow=optimal_flow)
display(fig)

#' Setting `mip_branchrule=2` and `mip_selectrule=3` seems to give
#' better results, according to the tuner.
#' Note that setting `outlev` to 0 switch off Knitro's output.
mip_solver = with_optimizer(KNITRO.Optimizer, mip_knapsack=0, mip_branchrule=2,
                            mip_selectrule=3, outlev=0)
model = Model(mip_solver)
load_mip_model!(model)
@time JuMP.optimize!(model)


#' ### Quantify impact of arcs' removals

#' The fewer arcs, the harder the problem.
#' If `nremovals >= 10`, the problem becomes infeasible.
max_removals = 9
# Save results in some arrays.
solve_time  = Float64[]
cost_values = Float64[]

for nremove in 1:max_removals
    println("Remove ", nremove, " arcs from graph.")
    model = Model(mip_solver)
    load_mip_model!(model, nremovals=nremove)
    @time JuMP.optimize!(model)

    push!(cost_values, JuMP.objective_value(model))

    # Plot!
    optimal_flow = JuMP.value.(model[:q])
    fig = figure()
    plot_network(flow=optimal_flow)
    title("Optimal solution with $nremove removals")
    display(fig)
end

#' Plot evolution of costs w.r.t. number of removals.
fig = figure()
plot(1:max_removals, cost_values, lw=3, c="k")
xlabel("#removals")
ylabel("Objective value")
display(fig)


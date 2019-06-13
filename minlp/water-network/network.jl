#' # Optimizing water networks with Knitro
#'
#' This notebook explains how to solve effectively NLP and MINLP
#' problem with the Knitro solver.
#'
#' The problem studied hereafter was originally designed by
#' Pierre Carpentier, for educational purposes.
#' The original code is available [here](http://perso.ensta-paristech.fr/~pcarpent/TP_Reseau/ENSMP/)
#' The author credits him for the following tutorial.
#'
#' # Description of the problem
#'
#' The tutorial studies the optimization of water flows in a water
#' distribution network (in a steady state). The optimizer aims at
#' minimizing the network's energy while statisfying some linear constraints
#' corresponding to the Kirchhoff first law.
#'
#' Let $G = (N, A)$ be an oriented graph. We note $n = |A|$ the number of
#' arcs and $m = |N|$ the number of nodes in this graph.
#'
#' We suppose that the network has $r$ tanks storing some waters to fulfill
#' the demands in $d$ nodes, distinct from the tanks. We split the set of nodes $N$ accordingly:
#' $$
#' N = N_r \subset N_d.
#' $$
#'
#' We suppose further that the graph is *connected*, implying
#' $$
#' n \geq m - 1 .
#' $$
#' We note $A \in R^{m \times n}$ the incidence matrix of the graph.
#'
#' Lets introduce physical variables to describe further the water network.
#' - $f$ (given), the vector of flows at nodes (that is, the demands)
#' - $p$ (given), the vector of pressures at nodes
#' - $r$ (given), the vector of resistances in all arcs
#' - $q$, the vector of flows across arcs
#'
#' The first Kirchhoff law states that:
#' $$
#' A q -f = 0,
#' $$
#' as we suppose that no accumulation occurs in the nodes.



#' ---
#' # Numerical resolution
#' We start by importing JuMP and Knitro.

using JuMP, KNITRO

#' Import data.
include("data.jl");

#' Plotting utilities
using Pkg

# Turned off by default.
PLOT_GRAPH = false
if haskey(Pkg.installed(), "PyPlot")
    PLOT_GRAPH = true
    include("utils.jl");
end


#' ---
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
# Add dummy variable in the model...
@variable(model, q[1:n])
# ... as we restrict the problem inside a 9-dimensional manifold.
@constraint(model, q .== q0 + B*qc)
@NLobjective(model, Min,
             sum(r[i] * abs(q[i]) * q[i]^2 / 3 + α1[i] * q[i] for i in 1:n))
optimize!(model)

#' Display results
if PLOT_GRAPH
    optimal_flow = JuMP.value.(q)
    fig = figure()
    plot_network(flow=optimal_flow)
    display(fig)
end

#' Observe that some arcs seem almost useless here...


#' ---
#' ## Extension to Mixed-integer non-linear programming

#' We were able to quickly formule the NLP problem with JuMP, and Knitro
#' finds a solution in few iterations.
#' We now modify slightly the problem. Consider the $n$ arcs inside the graph.
#' The questions are:
# - *how many arcs in $E$ can we remove before the problem becomes infeasible?*
# - *how much is the cost impacted by arc removals?*

#' We define hereafter the MINLP version of the problem.
# Add a maximum flows through the pipes
const QMAX = 10.

function load_mip_model!(model::JuMP.Model; nremovals=3)
    # Reconsider the NLP problem introduced previously.
    @variable(model, qc[1:nx])
    @variable(model, q[1:n])
    @constraint(model, q .== q0 + B*qc)

    # Introduce switch: $z$ is set to 0 if the corresponding arc is removed
    # from the graph.
    @variable(model, z[1:n], Bin)

    # Bounds the abs flows through the arcs by QMAX
    # The constraints write:
    # | q | <= q_{max} * z
    # Note that z_a = 0 implies that q_a = 0.
    @constraint(model,  q .<= QMAX * z)
    @constraint(model, -q .<= QMAX * z)
    # Ensure that we remove exactly nremovals arcs.
    @constraint(model, sum(z) == n - nremovals)

    # Same cost as previously.
    @NLobjective(model, Min,
                 sum(r[i] * abs(q[i]) * q[i]^2 / 3 + α1[i] * q[i] for i in 1:n))
    return
end

#' **Remark:** Note that the constraint
#' $$
#' | q | \leq q_{max} \; z
#' $$
#' could be reformulated as a complementarity constraint.


#' ### Default Knitro

#' Solve with default Knitro.
# Build non-linear solver.
model = Model(with_optimizer(KNITRO.Optimizer, outlev=3))
load_mip_model!(model)
@time JuMP.optimize!(model)

#' In the resolving, Knitro uses the Branch & Bound algorithm to find
#' the optimal solution corresponding to the (convex) MINLP problem.
#' We refer to the [documentation](https://www.artelys.com/docs/knitro/2_userGuide/minlp.html)
#' for further details.
#' The approach used here is different than in BONMIN (outer approximation).

#' Observe that Knitro default computes 165 nodes before finding the solution.
#' Is it possible to find a better tuning for Knitro?

#' ### MIP-Tuner

#' Since Knitro 12.0, a MINLP tuner was added to compute the optimal
#' parameterization to solve a given MINLP problem with Knitro. The MINLP
#' tuner uses an exhaustive search to find the optimal setting.
#'
#' By default, the tuner tests 36 combinations.
model = Model(with_optimizer(KNITRO.Optimizer, outlev=3, tuner=1))
load_mip_model!(model)
@time JuMP.optimize!(model)

# Plot!
if PLOT_GRAPH
    optimal_flow = JuMP.value.(model[:q])
    fig = figure()
    plot_network(flow=optimal_flow)
    display(fig)
end

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
    if PLOT_GRAPH
        optimal_flow = JuMP.value.(model[:q])
        fig = figure()
        plot_network(flow=optimal_flow)
        title("Optimal solution with $nremove removals")
        display(fig)
    end
end

#' Plot evolution of costs w.r.t. number of removals.
if PLOT_GRAPH
    fig = figure()
    plot(1:max_removals, cost_values, lw=3, c="k")
    xlabel("#removals")
    ylabel("Objective value")
    display(fig)
end

#' We observe that removing up to 6 arcs does not impact significantly
#' the cost.

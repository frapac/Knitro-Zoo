#' # Optimizing water networks with Knitro
#'
#' This notebook explains how to solve effectively NLP and MINLP
#' problem with the Knitro solver.
#'
#' The problem studied hereafter was originally designed by
#' Pierre Carpentier, for educational purposes.
#' The original code is available [here](http://perso.ensta-paristech.fr/~pcarpent/TP_Reseau/ENSMP/).
#' The author credits him for the following tutorial.
#'
#' # Description of the problem
#'
#' The tutorial studies the optimization of water flows in a water
#' distribution network (at steady state). The optimizer aims at
#' minimizing the network's energy while statisfying some linear constraints
#' corresponding to Kirchhoff's first law.
#'
#' Let $G = (N, E)$ be an oriented graph. We note $n = |E|$ the number of
#' arcs and $m = |N|$ the number of nodes in this graph.
#'
#' We suppose that the network has $r$ tanks storing some waters to fulfill
#' the demands in $d$ nodes (distinct from the tanks).
#' We split the set of nodes $N$ accordingly:
#' $$
#' N = N_r \cup N_d.
#' $$
#'
#' We suppose further that the graph is *connected*, implying
#' $$
#' n \geq m - 1 .
#' $$
#' We note $A \in R^{m \times n}$ the incidence matrix of the graph.
#'
#' Lets introduce physical variables to describe further the water network.
#' - the vector of resulting flows at nodes, denoted $f = (f_d, f_r)$. $f_d$ is given for demands' nodes.
#' - the vector of pressures at nodes, denoted $p = (p_d, p_r)$. Pressures $p_r$ are given for the reservoirs.
#' - the vector of resistances in all arcs, denoted $r$ (parameters of the problem)
#' - the vector of flows across arcs, denoted $q$.
#'
#' The decision variable is the vector of flows $q$.
#'
#' ### Constraints
#'
#' The first Kirchhoff law states that:
#' $$
#' A q -f = 0,
#' $$
#' as we suppose that no accumulation occurs in the nodes.
#' The second Kirchhoff law takes into account the losses in the pipes,
#' which is given on each arc by a function $\phi_{alpha}$ (corresponding
#' to the Colebrooks law):
#' $$
#' \phi_{\alpha}(q_\alpha) = r_\alpha q_\alpha | q_\alpha |
#' $$
#' The second Kirchhoff law writes, in a vectorial manner,
#' $$
#' A^\top p + r \circ q \circ | q | = 0
#' $$
#'
#' ### Objective
#' On each arc $\alpha$, we define the energy function $\Phi_\alpha$
#' as
#' $$
#' \Phi_\alpha(q_\alpha) = \dfrac{1}{3} r_\alpha q_\alpha^2 | q_\alpha |
#' $$
#' On the graph, the overall energy equates
#' $$
#' J(q, f_r) = \dfrac{1}{3} q^\top (r \circ q \circ | q |) + p_r^\top f_r
#' $$
#'
#' The global problem writes
#' $$
#' \min_{q} \dfrac{1}{3} q^\top (r \circ q \circ | q |) + p_r^\top f_r \qquad
#' s.t. \quad Aq  - f = 0
#' $$
#'
#' ### Reformulation
#' By applying some mathematical tricks, we are able to reformulate
#' the problem in the following manner.
#'
#' By considering some properties of the overall graph, we split the vector
#' of flows in two $q = (q_T, q_C)$, where $q_T$ depends linearly on $q_C$.
#' Then, we are able to prove the existence of a matrix $B$ and a fixed
#' vector $q_0$ such that the vector of flows on arcs writes
#' $$
#' q = q_0 + B \; q_C
#' $$
#' By using this formulation, we can reduce the dimension of the search
#' space by selecting as decision variable the subvector $q_C$ instead of $q$.
#'
#' Note that the problem becomes also constraint-free, as the vector $q$
#' statisfying the previous equation sastifies also the first Kirchhoff law.
#'
#' We reformulate the optimization problem as
#' $$
#' \min_{q_c} \dfrac{1}{3} q^\top (r \circ q \circ | q |) + p_r^\top f_r \qquad
#' s.t. \quad q = q_0 + B q_C
#' $$



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
model = Model(with_optimizer(KNITRO.Optimizer, outlev=3))

#' Write non-linear optimization problem first.
# Parameters.
α1 = Ar' * pr
# Dimension of the problem.
nx = n - md
@variable(model, qc[1:nx])
# Add dummy variable in the model...
@variable(model, q[1:n])
# ... as we restrict the problem inside a 9-dimensional manifold:
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

#' **Optimization related remark:** Note that the constraint
#' $$
#' | q | \leq q_{max} \; z
#' $$
#' could be reformulated as a complementarity constraint (also supported by Knitro
#' but not natively by JuMP).

#' **Physics related remark:** We are a bit loosy with the physics of the
#' problem if we remove the arcs in the previous manner. Imagine that for
#' a given arc $a \in E$, $z_a =0$ thus implying $q_a =0$. Then, if we
#' note $i$ and $j$ the two adjacent nodes, the physics tells us that
#' $$
#' p_i = 0 \, \quad p_j = 0 ,
#' $$
#' which is not the case if we solve the previous optimization problem.
#' However, the goal of this tutorial is purely pedagogical and as a
#' consequence we allow us to play a bit with the physics.


#' ### Default Knitro

#' Solve with default Knitro.
# Build non-linear solver.
model = Model(with_optimizer(KNITRO.Optimizer, outlev=3))
load_mip_model!(model)
@time JuMP.optimize!(model)

#' We plot the solution by displaying in red the arcs removed from the graph.
# Plot!
if PLOT_GRAPH
    optimal_flow = JuMP.value.(model[:q])
    fig = figure()
    plot_network(flow=optimal_flow)
    display(fig)
end

#' In the resolving, Knitro uses the Branch & Bound algorithm to find
#' the optimal solution corresponding to the (convex) MINLP problem.
#' We refer to the [documentation](https://www.artelys.com/docs/knitro/2_userGuide/minlp.html)
#' for further details.
#' The approach used here is different than in BONMIN (outer approximation).

#' Observe that Knitro with default options computes 165 nodes before
#' finding the solution.
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
#' If `nremovals >= 10`, Knitro takes too long to solve the problem.
#' However, the problem remains feasible till n_removals = 13.
#' Compute for instance the solution of the following problem:
feas_model = Model(with_optimizer(KNITRO.Optimizer))
@variable(feas_model, qc[1:nx])
@variable(feas_model, q[1:n])
@constraint(feas_model, q .== q0 + B*qc)
@variable(feas_model, z[1:n], Bin)
@constraint(feas_model,  q .<= QMAX * z)
@constraint(feas_model, -q .<= QMAX * z)
@objective(feas_model, Min, sum(z))


#' We know study the evolution of the objective cost w.r.t. the
#' number of removals.
max_removals = 9

# Save results in some arrays.
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

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


#' ## Extension to Mixed-integer non-linear programming

#' We define hereafter the MINLP version of the problem.
# Add a maximum flows through the pipes
QMAX = 10.

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

#' Setting `mip_branchrule=2` and `mip_selectrule=2` seems to give
#' better results.
model = JuMP.direct_model(KNITRO.Optimizer(mip_zerohalf=3))
load_mip_model!(model)
@time JuMP.optimize!(model)

#' The fewer arcs, the harder the problem.
# With 5 arcs removed.
model = Model(with_optimizer(KNITRO.Optimizer, outlev=0))
load_mip_model!(model, nremovals=5)
@time JuMP.optimize!(model)

# With 7 arcs removed.
model = Model(with_optimizer(KNITRO.Optimizer, outlev=0))
load_mip_model!(model, nremovals=7)
@time JuMP.optimize!(model)

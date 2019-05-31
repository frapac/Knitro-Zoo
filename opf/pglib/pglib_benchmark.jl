################################################################################
# BENCHMARKING PGLIB WITH POWERMODELS
################################################################################

"""
Usage:

```
include("opf.jl")

solver = with_optimizer(KNITRO.Optimizer, linsolver=4, presolve=0))
benchopf("/path/to/pglib", solver, relaxation=ACPPowerModel)
```

"""

using JuMP
using Memento
using PowerModels
using KNITRO

# Switch off warning in PowerModels.
setlevel!(getlogger(PowerModels), "error")


"Get instance's name."
function getname(s::String)
    sprime = split(s, ".")[1]
    case, orig = split(sprime, "_")[3:4]
    return case * orig
end

# Check that is of Matpower type
function select_mfile(x::Vector{String})
    mfiles = String[]
    for j in x
        if endswith(j,".m")
            push!(mfiles, j)
        end
    end
    return mfiles
end

"Check whether log directory exists."
function checklogdir()
    if !isdir("log")
        mkdir("log")
    end
end

"Load Knitro solver."
function knitro()
    return JuMP.with_optimizer(KNITRO.Optimizer,
                               linsolver=4,         # use MA27
                               presolve=0,          # switch-off presolve
                               outlev=3,            # print all iterations
                               outmode=2,           # output log in file
                               maxtime_cpu=15*60,   # 15mn maxtime
                               outname="knitro.log")
end

"""
Benchmark PGLIB with specified solver, upon
given relaxation (must belong to RELAXATIONS dict).

Return an array storing the following information:
```julia
    colnames = ["opf", "slv_status", "slv_cputime", "slv_obj"]
```

# Usage

```julia

solver = with_optimizer(KNITRO.Optimizer, linsolver=4)
benchopf("/path/to/pglib", solver, relaxation=ACPPowerModel)

```

"""
function benchopf(pglib_path, solver;
                  relaxation=ACPPowerModel)
    # Build log directory.
    checklogdir()
    # Import matpower files.
    opfnames = select_mfile(readdir(pglib_path))
    n_opf = length(opfnames)

    # Get solver's name.
    solver_name = "knitro"

    # Benchmark!
    status = 0
    results = zeros(n_opf, 3)
    for (idn ,name) in enumerate(opfnames)
        println("Benchmark $name...")
        res = run_opf("$(pglib_path)/$(name)", relaxation, solver)
        if res["termination_status"] == PowerModels.LOCALLY_SOLVED
            status = 0
        elseif res["termination_status"] == PowerModels.LOCALLY_INFEASIBLE
            status = 1
        else
            status = 3
        end
        results[idn, 1:3] = [status, res["solve_time"], res["objective"]]

        if isfile("$(solver_name).log")
            cp("$(solver_name).log",
               "log/$(solver_name)_$(getname(name)).log",
               force=true)
        end
    end

    return results
end

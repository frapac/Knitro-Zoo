using PyPlot

function plot_network(; flow=nothing)
    for i in 1:n
        i0, i1 = orig[i], dest[i]
        plot([absn[i0], absn[i1]], [ordn[i0], ordn[i1]],
             lw= isa(flow, Array{Float64}) ? 20*flow[i] : 5,
             zorder=1, color="darkblue")
        # If the flow is null, plot a red edge.
        if(isa(flow, Array{Float64}) && flow[i] == 0)
            plot([absn[i0], absn[i1]], [ordn[i0], ordn[i1]],
             lw= 5, zorder=1, color="red")
        end
    end
    scatter(absn, ordn, c="k", s=200, zorder=2)
end

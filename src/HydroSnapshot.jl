include("PrimitiveVars.jl")
import Base.copy

mutable struct HydroSnapshot
    grid::Array{Float64, 1}
    cells::Array{PrimitiveVars, 1}
end

function midGrid(grid)::Array{Float64, 1}
    return [0.5*(grid[n]+grid[n+1])
            for n in 1:length(grid)-1]
end

function initHydroSnapshot(grid, d_func, p_func, v_func)
    cells = [PrimitiveVars(d_func(x),
                           p_func(x),
                           v_func(x))
             for x in midGrid(grid)]
    return HydroSnapshot(grid, cells)
end

function copy(source::HydroSnapshot)::HydroSnapshot
    return HydroSnapshot(source.grid, source.cells)
end

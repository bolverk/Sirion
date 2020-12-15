struct Eulerian
end

function calcEdgeVelocities(vm::Eulerian,
                            state::HydroSnapshot)::Array{Float64, 1}
    return [0.0 for x in 1:length(state.grid)]
end

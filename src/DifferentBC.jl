struct DifferentBC
    left
    right
end

function calcLeftFlux(bc::DifferentBC,
                      state::HydroSnapshot,
                      edge_velocities::Array{Float64,1},
                      eos,
                      cached)::Conserved
    calcLeftFlux(bc.left,
                 state,
                 edge_velocities,
                 eos,
                 cached)
end

function calcRightFlux(bc::DifferentBC,
                       state::HydroSnapshot,
                       edge_velocities::Array{Float64,1},
                       eos,
                       cached)::Conserved
    calcRightFlux(bc.right,
                  state,
                  edge_velocities,
                  eos,
                  cached)
end

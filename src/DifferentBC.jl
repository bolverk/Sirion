struct DifferentBC
    left
    right
end

function calcLeftFlux(bc::DifferentBC,
                      state::HydroSnapshot,
                      edge_velocities::Array{Float64,1},
                      eos)::Conserved
    calcLeftFlux(bc.left,
                 state,
                 edge_velocities,
                 eos)
end

function calcRightFlux(bc::DifferentBC,
                       state::HydroSnapshot,
                       edge_velocities::Array{Float64,1},
                       eos)::Conserved
    calcRightFlux(bc.right,
                  state,
                  edge_velocities,
                  eos)
end

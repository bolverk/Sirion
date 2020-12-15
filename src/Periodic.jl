struct Periodic
    rs
end

function calcLeftFlux(bc::Periodic,
                      state::HydroSnapshot,
                      edge_velocities::Array{Float64,1},
                      eos)::Conserved
    return calcFlux(bc.rs,
                    state.cells[end],
                    state.cells[1],
                    edge_velocities[1],
                    eos)
end

function calcRightFlux(bc::Periodic,
                       state::HydroSnapshot,
                       edge_velocities::Array{Float64,1},
                       eos)::Conserved
    return calcFlux(bc.rs,
                    state.cells[end],
                    state.cells[1],
                    edge_velocities[end],
                    eos)
end

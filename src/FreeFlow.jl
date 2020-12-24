struct FreeFlow
    rs # Riemann solver
end

function calcLeftFlux(bc::FreeFlow,
                      state::HydroSnapshot,
                      edge_velocities::Array{Float64,1},
                      eos)::Conserved
    return calcFlux(bc.rs,
                    state.cells[1],
                    state.cells[1],
                    edge_velocities[1],
                    eos)
end

function calcRightFlux(bc::FreeFlow,
                       state::HydroSnapshot,
                       edge_velocities::Array{Float64,1},
                       eos)::Conserved
    return calcFlux(bc.rs,
                    state.cells[end],
                    state.cells[end],
                    edge_velocities[end],
                    eos)
end

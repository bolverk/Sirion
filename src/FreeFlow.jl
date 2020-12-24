struct FreeFlow
    rs # Riemann solver
end

function calcLeftFlux(bc::FreeFlow,
                      state::HydroSnapshot,
                      edge_velocities::Array{Float64,1},
                      eos,
                      cached)::Conserved
    rpcond = RPCond(state.cells[1].density,
                    state.cells[1].pressure,
                    state.cells[1].velocity,
                    cached["sound_speeds"][1],
                    cached["energies"][1])
    return calcFlux(bc.rs,
                    rpcond,
                    rpcond,
                    edge_velocities[1],
                    eos)
end

function calcRightFlux(bc::FreeFlow,
                       state::HydroSnapshot,
                       edge_velocities::Array{Float64,1},
                       eos,
                       cached)::Conserved
    rpcond = RPCond(state.cells[end].density,
                    state.cells[end].pressure,
                    state.cells[end].velocity,
                    cached["sound_speeds"][end],
                    cached["energies"][end])
    return calcFlux(bc.rs,
                    rpcond,
                    rpcond,
                    edge_velocities[end],
                    eos)
end

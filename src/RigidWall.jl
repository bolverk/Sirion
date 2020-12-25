struct RigidWall
    rs # Riemann solver
end

function reflect(p::PrimitiveVars,
                 velocity::Float64)::PrimitiveVars
    return PrimitiveVars(p.density,
                         p.pressure,
                         2*velocity-p.velocity)
end

function reflect(p::RPCond,
                 velocity::Float64)::RPCond
    return RPCond(p.density,
                  p.pressure,
                  2*velocity-p.velocity,
                  p.sound_speed,
                  p.energy)
end

function calcLeftFlux(bc::RigidWall,
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
                    reflect(rpcond,edge_velocities[1]),
                    rpcond,
                    edge_velocities[1],
                    eos)
end

function calcRightFlux(bc::RigidWall,
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
                    reflect(rpcond,edge_velocities[end]),
                    edge_velocities[end],
                    eos)
end

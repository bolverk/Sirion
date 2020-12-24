struct RigidWall
    rs # Riemann solver
end

function reflect(p::PrimitiveVars,
                 velocity::Float64)::PrimitiveVars
    return PrimitiveVars(p.density,
                         p.pressure,
                         2*velocity-p.velocity)
end

function calcLeftFlux(bc::RigidWall,
                      state::HydroSnapshot,
                      edge_velocities::Array{Float64,1},
                      eos,
                      cached)::Conserved
    return calcFlux(bc.rs,
                    reflect(state.cells[1],edge_velocities[1]),
                    state.cells[1],
                    edge_velocities[1],
                    eos)
end

function calcRightFlux(bc::RigidWall,
                       state::HydroSnapshot,
                       edge_velocities::Array{Float64,1},
                       eos,
                       cached)::Conserved
    return calcFlux(bc.rs,
                    state.cells[end],
                    reflect(state.cells[end],edge_velocities[end]),
                    edge_velocities[end],
                    eos)
end

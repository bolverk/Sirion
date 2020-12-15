struct SimpleFluxCalculator
    rs
    bc
end

function calcFluxes(fc::SimpleFluxCalculator,
                    state::HydroSnapshot,
                    edge_velocities,
                    eos)
    bulk_fluxes = [calcFlux(fc.rs,
                            state.cells[n],
                            state.cells[n+1],
                            edge_velocities[n+1],
                            eos)
                    for n in 1:length(state.cells)-1]
    left_bc = calcLeftFlux(fc.bc,
                           state,
                           edge_velocities,
                           eos)
    right_bc = calcRightFlux(fc.bc,
                             state,
                             edge_velocities,
                             eos)
    return [[left_bc];
            bulk_fluxes;
            [right_bc]]
end

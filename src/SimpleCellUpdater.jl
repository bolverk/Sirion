struct SimpleCellUpdater
end

function extensive2primitive(c::Conserved,
                             volume::Float64,
                             eos)
    density = c.mass/volume
    velocity = c.momentum/c.mass
    total_specific_energy = c.energy/c.mass
    kinetic_specific_energy = 0.5*velocity^2
    thermal_specific_energy = total_specific_energy-kinetic_specific_energy
    pressure = calcPressure(eos, density, thermal_specific_energy)
    return PrimitiveVars(density, pressure, velocity)
end

function updateCells!(cu::SimpleCellUpdater,
                      state::HydroSnapshot,
                      extensives::Array{Conserved, 1},
                      pg,
                      eos)
    cumulative_volumes = [calcVolume(pg, r) for r in state.grid]
    volumes = diff(cumulative_volumes)
    for n in 1:length(state.cells)
        state.cells[n] = extensive2primitive(extensives[n],
                                             volumes[n],
                                             eos)
    end
end

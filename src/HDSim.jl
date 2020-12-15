include("Conserved.jl")
include("HydroSnapshot.jl")

mutable struct HDSim
    time::Float64 # Time
    cycle::Int64 # Cycle number
    state::HydroSnapshot # Primitive variables and grid
    extensives::Array{Conserved, 1} # Conserved quantities
    eos # Equation of state
    tsf # Time step function
    vm # Vertex motion
    fc # Flux calculator
    pg # Physical geometry
    st # Source term
    cu # Cell updater
end

function primitive2extensive(p::PrimitiveVars,
                             volume::Float64,
                             eos)::Conserved
    mass = p.density*volume
    momentum = mass*p.velocity
    kinetic_specific_energy = 0.5*p.velocity^2
    thermal_specific_energy =
     calcSpecificThermalEnergy(eos,
                        p.density,
                        p.pressure)
    specific_energy =
     kinetic_specific_energy+
     thermal_specific_energy
    energy = specific_energy*mass
    return Conserved(mass,
                     momentum,
                     energy)
end

function initConserved(state::HydroSnapshot, eos)
    volumes = diff(state.grid)
    conserved_list = [
    primitive2extensive(p,V,eos)
    for (p,V) in zip(state.cells, volumes)]
    return conserved_list
end

function initHDSim(init_cond::HydroSnapshot,
                   eos,
                   tsf,
                   vm,
                   fc,
                   pg,
                   st,
                   cu)
    conserved_list = initConserved(init_cond, eos)
    return HDSim(0,
                 0,
                 init_cond,
                 conserved_list,
                 eos,
                 tsf,
                 vm,
                 fc,
                 pg,
                 st,
                 cu)
end

function updateConserved!(conserved_list::Array{Conserved,1},
                          fluxes::Array{Conserved,1},
                          dt::Float64)
    for n in 1:length(conserved_list)
        conserved_list[n] += dt*fluxes[n]
        conserved_list[n] -= dt*fluxes[n+1]
    end
end

function updateSourceContribution!(st,
                                  extensives::Array{Conserved,1},
                                  state::HydroSnapshot,
                                  t::Float64,
                                  dt::Float64)
    for (e,sc) in zip(extensives, calcSource(st,state,t,dt))
        e += dt*sc
    end
end

function updatePositions!(grid::Array{Float64, 1},
                          velocities::Array{Float64, 1},
                          dt::Float64)
    for n in 1:length(grid)
        grid[n] += dt*velocities[n]
    end
end

function timeAdvance(hdsim::HDSim)
    dt = calcTimeStep(hdsim.tsf, hdsim.state, hdsim.eos)
    edge_velocities = calcEdgeVelocities(hdsim.vm, hdsim.state)
    fluxes = calcFluxes(hdsim.fc,
                        hdsim.state,
                        edge_velocities,
                        hdsim.eos)
    updateConserved!(hdsim.extensives,
                     fluxes,
                     dt)
    updateSourceContribution!(hdsim.st,
                              hdsim.extensives,
                              hdsim.state,
                              hdsim.time,
                              dt)
    updatePositions!(hdsim.state.grid,
                     edge_velocities,
                     dt)

    updateCells!(hdsim.cu,
                 hdsim.state,
                 hdsim.extensives,
                 hdsim.pg,
                 hdsim.eos)

    hdsim.time += dt
    hdsim.cycle += 1
end

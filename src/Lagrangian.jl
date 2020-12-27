struct Lagrangian
    edge_behaviour
end

struct KeepEdgesFixed
end

function calcLeftVelocity(eb::KeepEdgesFixed, state::HydroSnapshot)
    return 0
end

function calcRightVelocity(eb::KeepEdgesFixed, state::HydroSnapshot)
    return 0
end

function calcEdgeVelocities(vm::Lagrangian,
                            state::HydroSnapshot)::Array{Float64, 1}
    left_velocity = calcLeftVelocity(vm.edge_behaviour, state)
    right_velocity = calcRightVelocity(vm.edge_behaviour, state)
    bulk_velocity = [0.5*(state.cells[n].velocity+
                          state.cells[n+1].velocity)
                     for n in 1:length(state.cells)-1]
    return [[left_velocity];
            bulk_velocity;
            [right_velocity]]
end

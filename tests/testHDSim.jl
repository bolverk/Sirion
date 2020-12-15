using Test

function calcRelativeDiff(x::Float64, y::Float64)::Float64
    return abs(x-y)/(abs(x)+abs(y))
end

function runRiemannProblem1()

    # Initialise simulation
    grid = range(0.0,1.0,length=100)
    d_func = x->1
    p_func = x->10-9*(sign(x-0.5)+1)/2
    v_func = x->0
    init_cond = Sirion.initHydroSnapshot(grid, d_func, p_func, v_func)
    eos = Sirion.IdealGas(5.0/3.0)
    tsf = Sirion.CFL(0.3)
    vm = Sirion.Eulerian()
    rs = Sirion.HLLC()
    bc = Sirion.RigidWall(rs)
    fc = Sirion.SimpleFluxCalculator(rs, bc)
    pg = Sirion.Planar()
    st = Sirion.ZeroSource()
    cu = Sirion.SimpleCellUpdater()
    sim = Sirion.initHDSim(init_cond,
                          eos,
                          tsf,
                          vm,
                          fc,
                          pg,
                          st,
                          cu)

    # Run simulation
    while sim.time<0.067
        Sirion.timeAdvance(sim)
    end

    # Perform diagnostics
    pressure_probe = sim.state.cells[length(sim.state.cells)÷2].pressure
    velocity_probe = sim.state.cells[length(sim.state.cells)÷2].velocity
    reference_pressure = 5.1123
    reference_velocity = 1.53795

    cond1 = calcRelativeDiff(reference_pressure, pressure_probe)<0.01
    cond2 = calcRelativeDiff(reference_velocity, velocity_probe)<0.01

    return cond1 && cond2
end

function runRiemannProblem2()

    # Initialise simulation
    grid = range(0.0,1.0,length=100)
    d_func = x->1
    p_func = x->1
    v_func = x->-5*sign(x-0.5)
    init_cond = Sirion.initHydroSnapshot(grid, d_func, p_func, v_func)
    eos = Sirion.IdealGas(5.0/3.0)
    tsf = Sirion.CFL(0.3)
    vm = Sirion.Eulerian()
    rs = Sirion.HLLC()
    bc = Sirion.RigidWall(rs)
    fc = Sirion.SimpleFluxCalculator(rs, bc)
    pg = Sirion.Planar()
    st = Sirion.ZeroSource()
    cu = Sirion.SimpleCellUpdater()
    sim = Sirion.initHDSim(init_cond,
                          eos,
                          tsf,
                          vm,
                          fc,
                          pg,
                          st,
                          cu)

    # Run simulation
    while sim.time<0.08
        Sirion.timeAdvance(sim)
    end

    # Perform diagnostics
    pressure_probe = sim.state.cells[1+length(sim.state.cells)÷2].pressure
    velocity_probe = sim.state.cells[1+length(sim.state.cells)÷2].velocity
    reference_pressure = 35.5397
    reference_velocity = 0.0

    cond1 = calcRelativeDiff(reference_pressure, pressure_probe)<0.1
    cond2 = abs(velocity_probe)<0.01

    return cond1 && cond2
end

@test runRiemannProblem1()
@test runRiemannProblem2()

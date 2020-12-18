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

function runRiemannProblem3()

    # Initialise simulation
    grid = range(0.0,1.0,length=100)
    d_func = x->1
    p_func = x->1
    v_func = x->1.0*sign(x-0.5)
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
    while sim.time<0.22
        Sirion.timeAdvance(sim)
    end

    # Perform diagnostics
    pressure_probe = sim.state.cells[1+length(sim.state.cells)÷2].pressure
    velocity_probe = sim.state.cells[1+length(sim.state.cells)÷2].velocity
    reference_pressure = 0.224614
    reference_velocity = 0.0

    cond1 = calcRelativeDiff(reference_pressure, pressure_probe)<0.01
    cond2 = abs(velocity_probe)<0.01

    return cond1 && cond2
end

function l1_fit(a1, a2)
    abs_dif = [abs(x-y) for (x,y) in zip(a1, a2)]
    return sum(abs_dif)/length(abs_dif)
end

function runPlanarNoh()

    # Initialise simulation
    grid = range(0.0,1.0,length=100)
    d_func = x->1
    p_func = x->1e-6
    v_func = x->-1.0
    init_cond = Sirion.initHydroSnapshot(grid, d_func, p_func, v_func)
    g = 5.0/3.0
    eos = Sirion.IdealGas(g)
    tsf = Sirion.CFL(0.3)
    vm = Sirion.Eulerian()
    rs = Sirion.HLLC()
    bc = Sirion.DifferentBC(Sirion.RigidWall(rs),
                            Sirion.FreeFlow(rs))
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
    while sim.time<2.0
        Sirion.timeAdvance(sim)
    end

    # Perform diagnostics
    upstream_speed = 1.0
    shock_speed = (g-1)*upstream_speed/2
    shock_position = shock_speed*sim.time
    upstream_density = 1
    downstream_density = (g+1)*upstream_density/(g-1)
    upstream_velocity = -1
    downstream_velocity = 0
    upstream_pressure = 1e-6
    downstream_pressure = ((g+1)/2)*upstream_density*
        (downstream_velocity-upstream_velocity)^2
    d_func2(x) = x<shock_position ? downstream_density : upstream_density
    p_func2(x) = x<shock_position ? downstream_pressure : upstream_pressure
    v_func2(x) = x<shock_position ? downstream_velocity : upstream_velocity
    analytic_d = [d_func2(x) for x in Sirion.midGrid(sim.state.grid)]
    analytic_p = [p_func2(x) for x in Sirion.midGrid(sim.state.grid)]
    analytic_v = [v_func2(x) for x in Sirion.midGrid(sim.state.grid)]
    cond1 = l1_fit(analytic_d, [itm.density for itm in sim.state.cells])<0.03
    cond2 = l1_fit(analytic_d, [itm.density for itm in sim.state.cells])<0.03
    cond3 = l1_fit(analytic_v, [itm.velocity for itm in sim.state.cells])<0.01
    return (cond1 && cond2 && cond3)

end

@test runRiemannProblem1()
@test runRiemannProblem2()
@test runRiemannProblem3()
@test runPlanarNoh()

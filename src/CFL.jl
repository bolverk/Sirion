struct CFL
    prefactor::Float64
end

function calcSingleTimeStep(cell::PrimitiveVars,
                            dx::Float64,
                            eos)::Float64
    sound_speed = calcSoundSpeed(eos,
                                 cell.density,
                                 cell.pressure)
    return dx/(sound_speed+abs(cell.velocity))
end

function calcTimeStep(tsf::CFL,
                      state::HydroSnapshot,
                      eos)
    dx_list = diff(state.grid)
    dt_list = [calcSingleTimeStep(cell, dx, eos)
               for (cell, dx)
               in zip(state.cells, dx_list)]
    return minimum(dt_list)*tsf.prefactor
end

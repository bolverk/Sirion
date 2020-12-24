struct CFL
    prefactor::Float64
end

function calcSingleTimeStep(cell::PrimitiveVars,
                            dx::Float64,
                            ss::Float64)::Float64
    return dx/(ss+abs(cell.velocity))
end

function calcTimeStep(tsf::CFL,
                      state::HydroSnapshot,
                      eos,
                      cached)
    dx_list = diff(state.grid)
    dt_list = [calcSingleTimeStep(cell, dx, ss)
               for (cell, dx, ss)
               in zip(state.cells, dx_list, cached["sound_speeds"])]
    return minimum(dt_list)*tsf.prefactor
end

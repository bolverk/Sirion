struct ZeroSource
end

function calcSource(source::ZeroSource,
                    state::HydroSnapshot,
                    t::Float64,
                    dt::Float64)::Array{Conserved,1}
    return [Conserved(0,0,0) for x in state.cells]
end

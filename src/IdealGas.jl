struct IdealGas
    AdiabaticIndex
end

function calcSoundSpeed(eos::IdealGas,
                        d::Float64,
                        p::Float64)
    g = eos.AdiabaticIndex
    return sqrt(g*p/d)
end

function calcSpecificThermalEnergy(eos::IdealGas,
                                   d::Float64,
                                   p::Float64)
    g = eos.AdiabaticIndex
    return p/d/(g-1)
end

function calcPressure(eos::IdealGas,
                      d::Float64,
                      e::Float64)
    g = eos.AdiabaticIndex
    return (g-1)*e*d                  
end

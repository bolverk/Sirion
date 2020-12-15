using Test

function testBackAndForthPressureEnergy()
    d_list = [10.0^n for n in -5:5]
    p_list = [10.0^n for n in -5:5]
    eos = Sirion.IdealGas(5.0/3.0)
    for d in d_list
        for p in p_list
            e = Sirion.calcSpecificThermalEnergy(eos, d, p)
            pr = Sirion.calcPressure(eos, d, e)
            if !(pr ≈ p)
                return false
            end
        end
    end
    return true
end

@test testBackAndForthPressureEnergy()

import Base.-
import Base.*
import Base.+

mutable struct Conserved
    mass::Float64
    momentum::Float64
    energy::Float64
end

function +(argl::Conserved, argr::Conserved)
    return Conserved(argl.mass+argr.mass,
                     argl.momentum+argr.momentum,
                     argl.energy+argr.energy)
end

function *(s::Float64, c::Conserved)
    return Conserved(s*c.mass,
                     s*c.momentum,
                     s*c.energy)
end

function *(s::Int64, c::Conserved)
    return convert(Float64, s)*c
end

function -(argl::Conserved, argr::Conserved)
    return argl+(-1)*argr
end

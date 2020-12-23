struct Planar
end

function calcArea(pg::Planar, r::Float64)::Float64
    return 1
end

function calcVolume(pg::Planar, r::Float64)::Float64
    return r
end

struct Cylindrical
end

function calcArea(pg::Cylindrical, r::Float64)::Float64
    return 2*pi*r
end

function calcVolume(pg::Cylindrical, r::Float64)::Float64
    return pi*r^2
end

struct Spherical
end

function calcArea(pg::Spherical, r::Float64)::Float64
    return 4*pi*r^2
end

function calcVolume(pg::Spherical, r::Float64)::Float64
    return 4*pi*r^3/3
end

struct SphericalComplamentary
end

function calcSource(source::SphericalComplamentary,
                    state::HydroSnapshot,
                    t::Float64,
                    dt::Float64)::Array{Conserved,1}
    cumulative_volumes = [4*pi*r^3/3 for r in state.grid]
    volumes = diff(cumulative_volumes)
    r_list = midGrid(state.grid)
    return [Conserved(0,2*volume*c.pressure/r,0)
            for (c,r,volume) in zip(state.cells,r_list,volumes)]
end

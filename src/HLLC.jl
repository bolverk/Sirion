include("Conserved.jl")
include("RPCond.jl")

struct HLLC
end

function calcWaveSpeeds(left::RPCond,
	 				    right::RPCond)::Array{Float64,1}
	dl = left.density;
	pl = left.pressure;
	vl = left.velocity;
	cl = left.sound_speed
    dr = right.density;
	pr = right.pressure;
	vr = right.velocity;
	cr = right.sound_speed
	sl = min(vl - cl, vr - cr);
	sr = max(vl + cl, vr + cr);
	ss = (pr - pl + dl*vl*(sl - vl) - dr*vr*(sr - vr)) /
			(dl*(sl - vl) - dr*(sr - vr));
	return [sl, ss, sr];
end

function calcTotalEnergyDensity(state::RPCond)::Float64
	d = state.density
	p = state.pressure
	v = state.velocity
	e = state.energy
	return d*(e+v^2/2)
end

function calcStarredState(state::RPCond,
	 					  sk::Float64,
						  ss::Float64)::Conserved
	dk = state.density
	pk = state.pressure
	vk = state.velocity
	ds = dk*(sk - vk) / (sk - ss)
	ek = calcTotalEnergyDensity(state)
	mass = ds
	momentum = ds*ss
	energy = ek*ds / dk +
			ds*(ss - vk)*(ss + pk / dk / (sk - vk));
	return Conserved(mass, momentum, energy)
end

function calcMassDensity(state::RPCond)::Float64
	return state.density
end

function calcMomentumDensity(state::RPCond)::Float64
	return state.density*state.velocity
end

function primitive2conserved(state::RPCond)::Conserved
	return Conserved(calcMassDensity(state),
					 calcMomentumDensity(state),
					 calcTotalEnergyDensity(state))
end

function calcMassFlux(state::RPCond)::Float64
	return state.density*state.velocity
end

function calcMomentumFlux(state::RPCond)::Float64
	return state.density*state.velocity^2+state.pressure
end

function calcEnergyFlux(state::RPCond)::Float64
	energy_density = calcTotalEnergyDensity(state)
	enthalpy = energy_density+state.pressure
	return enthalpy*state.velocity
end

function primitive2flux(state::RPCond)::Conserved
	return Conserved(calcMassFlux(state),
					 calcMomentumFlux(state),
					 calcEnergyFlux(state))
end

function calcRestFrameFluxes(left::RPCond,
							 right::RPCond)::Conserved

	ws = calcWaveSpeeds(left, right)

	ul = primitive2conserved(left)
	ur = primitive2conserved(right)

	fl = primitive2flux(left)
	fr = primitive2flux(right)

	usl = calcStarredState(left, ws[1], ws[2])
	usr = calcStarredState(right, ws[3], ws[2])

	if ws[1]>0
		return fl
	elseif ws[2]>=0
		return fl + ws[1]*(usl-ul)
	elseif ws[3]>=0
		return fr + ws[3]*(usr-ur)
	else
		return fr
	end
end

function boostFlux(flux::Conserved,
				   velocity::Float64)::Conserved
	return Conserved(flux.mass,
					 flux.momentum+
					 velocity*flux.mass,
					 flux.energy+
					 0.5*flux.mass*velocity^2+
					 flux.momentum*velocity)
end

function boostCell(p::RPCond,
				   velocity::Float64)::RPCond
	return RPCond(p.density,
				  p.pressure,
				  p.velocity-velocity,
				  p.sound_speed,
				  p.energy)
end

function calcFlux(rs::HLLC,
				  left::RPCond,
				  right::RPCond,
				  velocity::Float64,
				  eos)::Conserved

	return boostFlux(calcRestFrameFluxes(boostCell(left,-velocity),
										 boostCell(right,-velocity)),
					 velocity)
end

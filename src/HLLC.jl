include("Conserved.jl")
include("RPCond.jl")

struct HLLC
end

function calcWaveSpeeds(left::PrimitiveVars,
	 				    right::PrimitiveVars,
						eos)::Array{Float64,1}
	dl = left.density;
	pl = left.pressure;
	vl = left.velocity;
	cl = calcSoundSpeed(eos,dl,pl)
    dr = right.density;
	pr = right.pressure;
	vr = right.velocity;
	cr = calcSoundSpeed(eos, dr, pr)
	sl = min(vl - cl, vr - cr);
	sr = max(vl + cl, vr + cr);
	ss = (pr - pl + dl*vl*(sl - vl) - dr*vr*(sr - vr)) /
			(dl*(sl - vl) - dr*(sr - vr));
	return [sl, ss, sr];
end

function calcWaveSpeeds(left::RPCond,
	 				    right::RPCond,
						eos)::Array{Float64,1}
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

function calcTotalEnergyDensity(state::PrimitiveVars, eos)::Float64
	d = state.density
	p = state.pressure
	v = state.velocity
	e = calcSpecificThermalEnergy(eos, d, p)
	return d*(e+v^2/2)
end

function calcTotalEnergyDensity(state::RPCond, eos)::Float64
	d = state.density
	p = state.pressure
	v = state.velocity
	e = state.energy
	return d*(e+v^2/2)
end

function calcStarredState(state,
	 					  sk::Float64,
						  ss::Float64,
						  eos)::Conserved
	dk = state.density
	pk = state.pressure
	vk = state.velocity
	ds = dk*(sk - vk) / (sk - ss)
	ek = calcTotalEnergyDensity(state, eos)
	mass = ds
	momentum = ds*ss
	energy = ek*ds / dk +
			ds*(ss - vk)*(ss + pk / dk / (sk - vk));
	return Conserved(mass, momentum, energy)
end

function calcMassDensity(state::PrimitiveVars,
						 eos)::Float64
	return state.density
end

function calcMassDensity(state::RPCond,
						 eos)::Float64
	return state.density
end

function calcMomentumDensity(state::PrimitiveVars,
	 						 eos)::Float64
	return state.density*state.velocity
end

function calcMomentumDensity(state::RPCond,
	 						 eos)::Float64
	return state.density*state.velocity
end

function primitive2conserved(state,
							 eos)::Conserved
	return Conserved(calcMassDensity(state, eos),
					 calcMomentumDensity(state, eos),
					 calcTotalEnergyDensity(state, eos))
end

function calcMassFlux(state,
	 				  eos)::Float64
	return state.density*state.velocity
end

function calcMomentumFlux(state,
	 					  eos)::Float64
	return state.density*state.velocity^2+state.pressure
end

function calcEnergyFlux(state,
					    eos)::Float64
	energy_density = calcTotalEnergyDensity(state,eos)
	enthalpy = energy_density+state.pressure
	return enthalpy*state.velocity
end

function primitive2flux(state,
						eos)::Conserved
	return Conserved(calcMassFlux(state, eos),
					 calcMomentumFlux(state, eos),
					 calcEnergyFlux(state, eos))
end

function calcRestFrameFluxes(left,
							 right,
							 eos)::Conserved

	ws = calcWaveSpeeds(left, right, eos)

	ul = primitive2conserved(left, eos)
	ur = primitive2conserved(right, eos)

	fl = primitive2flux(left, eos)
	fr = primitive2flux(right, eos)

	usl = calcStarredState(left, ws[1], ws[2], eos)
	usr = calcStarredState(right, ws[3], ws[2], eos)

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

function boostCell(p::PrimitiveVars,
				   velocity::Float64)::PrimitiveVars
	return PrimitiveVars(p.density,
					 	 p.pressure,
					 	 p.velocity-velocity)
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
				  left,
				  right,
				  velocity::Float64,
				  eos)::Conserved

	return boostFlux(calcRestFrameFluxes(boostCell(left,-velocity),
										 boostCell(right,-velocity),
										 eos),
					 velocity)
end

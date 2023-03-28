
import OptimPackNextGen.Powell.Newuoa.newuoa

using LinearAlgebra
using OptimPackNextGen
using StatsBase

@enum Side FT=0 SC=16
@enum Diode D1=1 D2=2 D3=3 D4=4 FC

const ȷ=im

#include("Faint.jl")


@enum MetState OFF=0 LOW=1 NORMAL=2 HIGH=3 TRANSIENT=-1
 
struct FaintStates{T<:AbstractFloat,A<:AbstractVector{T}}
	timer1::A
	timer2::A
	voltage1::T
	voltage2::T
	state1::MetState
	state2::MetState
end
 
function  FaintStates(state1::AbstractVector{T},state2::AbstractVector{T},voltage1,voltage2) where {T<:AbstractFloat}
	A = typeof(state1)
	if  voltage1 > voltage2  # LOW > HIGH
	 	return FaintStates{T,A}(state2,state1,voltage2,voltage1,HIGH,LOW)
	end
	return FaintStates{T,A}(state1,state2,voltage1,voltage2,HIGH, LOW)

end

function buildstates(faintstates::FaintStates{T,A},timestamp::AbstractVector; lag::Integer=0,preswitchdelay =0,postwitchdelay =0) where {T<:AbstractFloat,A<:AbstractVector{T}}
	N = length(timestamp)
	states = Vector{MetState}(undef,N)
	timestep = timestamp[2]-timestamp[1]
	t1 = collect(faintstates.timer1 .+ lag*timestep)
	t2 = collect(faintstates.timer2 .+ lag*timestep)

	bounds(x) = floor(Int,x / (2*timestep)),ceil( Int, x / (2*timestep))
	(premin,premax) = bounds(preswitchdelay)
	(postmin,postmax) = bounds(postwitchdelay)

	currentstate = NORMAL
	first1 = popfirst!(t1)
	first2 = popfirst!(t2)

	@inbounds @simd for index ∈ 1:N
		time = timestamp[index]
		forget = 0
		# HIGH
		if time >= first1  
			currentstate = faintstates.state1
			states[index-premin:index] .= TRANSIENT
			forget = premax
			if isempty(t1) 
				first1 = last(timestamp)
				if first2 == last(timestamp)
					currentstate = NORMAL
				end
			else
				first1 = popfirst!(t1)
			end
		end
		# LOW
		if time >= first2
			currentstate = faintstates.state2
			states[index-postmin:index] .= TRANSIENT
			forget = postmax
			if isempty(t2) 
				first2 = last(timestamp)
				if first1 == last(timestamp)
					currentstate = NORMAL
				end
			else
			first2 = popfirst!(t2)
			end
		end
		if( forget>0)
			states[index] = TRANSIENT
			forget -=1
		else
			states[index] = currentstate
		end
	end	
	return states
end

function estimatelag(states::Vector{MetState} ,data::AbstractVector{Complex{T}}; range::AbstractVector{Int64}=-10:10) where {T<:AbstractFloat}
	m = [mean(abs2,data[circshift(states,i) .== HIGH]) for i ∈ range];
	return range[argmax(m)]
end

function compute_mean_power(states::Vector{MetState} ,data::AbstractVector{Complex{T}}) where {T<:AbstractFloat}
	pow = zeros(T,length(data))
	 for st ∈ instances(MetState)
		idx = states.==st
		pow[idx] .= mean(abs2,data[idx]);
	 end
	 return pow
end

function idx(side_::Side,telescope_::Integer, diode_::Diode) 
	if diode_==FC
		return 32 + Integer(side_)÷4  + (telescope_-1) + 1
	end
	return Integer(side_) + (Integer(diode_)-1) + (telescope_-1)*4 + 1 
end

mutable struct Modulation{T<:AbstractFloat}
	c::Complex{T}
	a::Complex{T}
	b::T
	ϕ::T
	ω::T
end

function Modulation(;T=Float32,c=0.,a=1.,b=1.,ϕ=0.,ω=2π)
	return Modulation{T}(c,a,b,ϕ, ω)
end

function Modulation{T}(mod::Modulation) where{T<:AbstractFloat}
	return Modulation{T}(convert.(Complex{T},(mod.c,mod.a))...,convert.(T,(mod.b,mod.ϕ, mod.ω))...)
end

function (self::Modulation{T})(timestamp::AbstractVector) where{T<:AbstractFloat,}
	timestamp = T.(timestamp)
	return self.c .+ self.a .* exp.(ȷ .* self.b .* sin.( self.ω .* timestamp .+ self.ϕ ))
end

function getphase(self::Modulation{T},timestamp::AbstractVector) where{T<:AbstractFloat}
	timestamp = T.(timestamp)
	return self.b .* sin.( self.ω .* timestamp .+ self.ϕ ) .+ angle(self.a)
end

function updatemodulation(self::Modulation{T}, timestamp::AbstractVector, data::AbstractVector{Complex{T}},  b::T, ϕ::T) where{T<:AbstractFloat}
	self.b = b
	self.ϕ = ϕ
	if b==0.
		self.c = 0
		self.a = mean(data)
		return self.a .* ones(Complex{T},size(data))
	end
	model = exp.(ȷ .* (b .* sin.( self.ω .* T.(timestamp) .+ ϕ )))
	(self.c, self.a) = linearregression( model, data)
	return 	self.c .+ self.a .* model
end

function updatemodulation!(self::Modulation{T}, model::Vector{Complex{T}},  timestamp::AbstractVector, data::AbstractVector{Complex{T}},  b::T, ϕ::T) where{T<:AbstractFloat}
	self.b = b
	self.ϕ = ϕ
	if b==0.
		self.c = 0
		self.a = mean(data)
		fill(self.a, model)
	else
		@. model = exp(ȷ * (b * sin( self.ω * timestamp + ϕ )))
		(self.c, self.a) = linearregression( model, data)
		@. model = 	self.c + self.a * model
	end
	return model
end

function linearregression( model::Vector{Complex{T}}, data::AbstractVector{Complex{T}}) where{T<:AbstractFloat}
	N = length(model)
	Sv1 = sum(data)
	Sv2 = model ⋅ data
	H12 = mean(model)
	detH = 1/(N*(1 - abs2(H12)))
	c = (Sv1 .- Sv2 .* H12) * detH
	a = (-Sv1 .* conj(H12) + Sv2) * detH
	return (c, a)
end


struct Chi2CostFunction{T<:AbstractFloat,P<:Union{Vector{T}, T}}
    N::Int64
    mod::Modulation{T}
    timestamp::Vector{T}
    data::Vector{Complex{T}}
    power::P
    function Chi2CostFunction{T,P}(mod::Modulation{T},
        							timestamp::Vector{T},
        							data::Vector{Complex{T}},
									power::P) where {T<:AbstractFloat,P<:Union{Vector{T}, T}}
        N =length(timestamp);
        @assert N == size(data,1) "voltage and time must have the same number of lines"
		if P <: Vector
			@assert N == size(power,1) "power and time must have the same number of lines"
			return new{T,Vector{T}}(N,mod,timestamp,data,power)
		end
        return new{T,T}(N,mod,timestamp,data,power)
    end
end

function Chi2CostFunction(timestamp::AbstractVector,data::AbstractVector{Complex{T}}; kwd...) where {T<:AbstractFloat}
	return Chi2CostFunction{T,T}(Modulation(;T=T,kwd...),T.(timestamp,)data,T.(1.0))
end

function Chi2CostFunction(timestamp::AbstractVector,data::AbstractVector{Complex{T}},power::Vector; kwd...) where {T<:AbstractFloat}
	return Chi2CostFunction{T,Vector{T}}(Modulation(;T=T,kwd...),T.(timestamp),data,T.(power))
end

function Chi2CostFunction(timestamp::AbstractVector,data::AbstractVector{Complex{T}},power::Number; kwd...) where {T<:AbstractFloat}
	return Chi2CostFunction{T,T}(Modulation(;T=T,kwd...),T.(timestamp),data,T.(power))
end

function (self::Chi2CostFunction{T})(b::T,ϕ::T) where{T<:AbstractFloat}
	pupilmodulation = updatemodulation(self.mod, self.timestamp, self.data, b, ϕ)
	return sum(abs2,( pupilmodulation .- self.data))
end

function (self::Chi2CostFunction{T})(pupilmodulation::AbstractVector{Complex{T}},b::T,ϕ::T) where{T<:AbstractFloat}
	updatemodulation!(self.mod, pupilmodulation, self.timestamp, self.data, b, ϕ)
	return sum(abs2,( self.power .* pupilmodulation .-  self.data))
end

(self::Chi2CostFunction{T})() where{T<:AbstractFloat}  = self(self.mod.b,self.mod.ϕ)
(self::Chi2CostFunction{T})(x::Vector{T}) where{T<:AbstractFloat}  = self(x[1],x[2])
(self::Chi2CostFunction{T})(scratch::Vector{Complex{T}},x::AbstractVector{T}) where{T<:AbstractFloat}  = self(scratch,x[1],x[2])

function minimize!(self::Chi2CostFunction{T}; xinit=[2,0]) where {T<:AbstractFloat}
	scratch =Vector{Complex{T}}(undef,self.N)
	xinit = T.(xinit)
	(status, x, χ2) =  newuoa(x ->self(scratch,x) , xinit,1,1e-3,maxeval=1500)
	return x

	# using Nelder Mead simplex method from Optim.jl is 3x slower than newoa
	# res = optimize(x ->self(scratch,x) , xinit, NelderMead())
	# return Optim.minimizer(res)
	
end

function demodulateall( timestamp::AbstractVector,data::AbstractMatrix{Complex{T}}; xinit=[0.01,0],recenter=true,faintparam::Union{Nothing,FaintStates} = nothing,onlyhigh=false,preswitchdelay=0,postwitchdelay=0)   where{T<:AbstractFloat}

	output = copy(data)
	param = Vector{Modulation{T}}(undef,32) 
	likelihood =  Vector{T}(undef,32) 
	ϕrange= range(0,π,5)

	if !isnothing(faintparam)
		state= buildstates(faintparam, timestamp)
		lag = estimatelag(state,data[:,idx(SC,1,FC)])
		@info "lag = $lag"
		state= buildstates(faintparam, timestamp; lag=lag, preswitchdelay=preswitchdelay,postwitchdelay=postwitchdelay)
	end

	Threads.@threads for (j,k) ∈ collect(Iterators.product(1:4,(FT,SC)))
		FCphase = angle.(data[:,idx(k,j,FC)])
		for i ∈ (D1,D2,D3,D4)
			d = view(data,:,idx(k,j,i)) .*  exp.(-1im.*FCphase)
			valid = (:)
			if !isnothing(faintparam) 
				if onlyhigh
					valid =  (state.== HIGH)
				end
				if any(x-> x == TRANSIENT,state) 
					valid =  (state .!= TRANSIENT)
				end
				power = compute_mean_power(state,view(data,:,idx(k,j,i)))[valid]
			else
				power = T.(1.)
			end
			lkl = Chi2CostFunction(timestamp[valid],d[valid],power)

			if xinit==:auto
				binit= initialguess(d)
				ϕinit = ϕrange[argmin(map(ϕ -> lkl(binit,ϕ),ϕrange ))]
				xinit=[binit, ϕinit]
			end
			x = minimize!(lkl,xinit=xinit)
			likelihood[idx(k,j,i)] = lkl(x)
			if recenter
				@. output[:,idx(k,j,i)] = (d  - lkl.mod.c) * exp(-1im*( $(getphase(lkl.mod, timestamp)) - FCphase- angle(lkl.mod.a)))
			else
				@. output[:,idx(k,j,i)] = data[:,idx(k,j,i)] * exp(-1im*( angle($(lkl.mod(timestamp)))))
			end
			param[idx(k,j,i)] = lkl.mod
		end
	end
	return (output, param,likelihood)
end

function initialguess(data::AbstractVector{Complex{T}}) where{T<:AbstractFloat}
	std(angle.(data .* exp.(-1ȷ .*angle(mean(data)))))
end

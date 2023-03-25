
import OptimPackNextGen.Powell.Newuoa.newuoa

using LinearAlgebra
using OptimPackNextGen
using StatsBase

@enum Side FT=0 SC=16
@enum Diode D1=1 D2=2 D3=3 D4=4 FC

const ȷ=im


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

function (self::Modulation{T})(timestamp::Vector{T2}) where{T<:AbstractFloat,T2}
	timestamp = T.(timestamp)
	return self.c .+ self.a .* exp.(ȷ .* self.b .* sin.( self.ω .* timestamp .+ self.ϕ ))
end

function getphase(self::Modulation{T},timestamp::Vector{T2}) where{T<:AbstractFloat,T2}
	timestamp = T.(timestamp)
	return self.b .* sin.( self.ω .* timestamp .+ self.ϕ ) .+ angle(self.a)
end

function updatemodulation(self::Modulation{T}, timestamp::Vector{T}, data::Vector{Complex{T}},  b::T, ϕ::T) where{T<:AbstractFloat}
	self.b = b
	self.ϕ = ϕ
	if b==0.
		self.c = 0
		self.a = mean(data)
		return self.a .* ones(Complex{T},size(data))
	end
	model = exp.(ȷ .* (b .* sin.( self.ω .* timestamp .+ ϕ )))
	(self.c, self.a) = linearregression( model, data)
	return 	self.c .+ self.a .* model
end

function updatemodulation!( model::Vector{Complex{T}}, self::Modulation{T}, timestamp::Vector{T}, data::Vector{Complex{T}},  b::T, ϕ::T) where{T<:AbstractFloat}
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

function linearregression( model::Vector{Complex{T}}, data::Vector{Complex{T}}) where{T<:AbstractFloat}
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
        							timestamp::Vector,
        							data::Vector{Complex{T}},
									power::P) where {T<:AbstractFloat,P<:Union{Vector, T}}
        N =length(timestamp);
        @assert N == size(data,1) "voltage and time must have the same number of lines"
		if typeof(P)===Vector{T}
			@assert N == size(power,1) "power and time must have the same number of lines"
		end
        return new{T,P}(N,mod,convert.(T,timestamp),data,convert.(T,power))
    end
end

function Chi2CostFunction(timestamp::Vector,data::Vector{Complex{T}}; kwd...) where {T<:AbstractFloat}
	return Chi2CostFunction{T,T}(Modulation(;T=T,kwd...),timestamp,data,1.0)
end

function Chi2CostFunction(timestamp::Vector,data::Vector{Complex{T}},power::Vector; kwd...) where {T<:AbstractFloat}
	return Chi2CostFunction{T,Vector{T}}(Modulation(;T=T,kwd...),timestamp,data,power)
end

function (self::Chi2CostFunction{T})(b::T,ϕ::T) where{T<:AbstractFloat}
	pupilmodulation = updatemodulation(self.mod, self.timestamp, self.data, b, ϕ)
	return sum(abs2,( pupilmodulation .- self.data))
end

function Chi2CostFunction(pupilmodulation::Vector{Complex{T}},timestamp::Vector,data::Vector{Complex{T}}; kwd...) where {T<:AbstractFloat}
	return Chi2CostFunction{T}(pupilmodulation,Modulation(;T=T,kwd...),timestamp,data)
end

function (self::Chi2CostFunction{T})(pupilmodulation::Vector{Complex{T}},b::T,ϕ::T) where{T<:AbstractFloat}
	updatemodulation!(pupilmodulation, self.mod, self.timestamp, self.data, b, ϕ)
	return sum(abs2,( self.power .* pupilmodulation .-  self.data))
end

(self::Chi2CostFunction{T})() where{T<:AbstractFloat}  = self(self.mod.b,self.mod.ϕ)
(self::Chi2CostFunction{T})(x::Vector{T}) where{T<:AbstractFloat}  = self(x[1],x[2])
(self::Chi2CostFunction{T})(scratch::Vector{Complex{T}},x::Vector{T}) where{T<:AbstractFloat}  = self(scratch,x[1],x[2])

function minimize!(self::Chi2CostFunction{T}; xinit=[2,0]) where {T<:AbstractFloat}
	scratch =Vector{Complex{T}}(undef,self.N)
	xinit = T.(xinit)
	(status, x, χ2) =  newuoa(x ->self(scratch,x) , xinit,1,1e-3,maxeval=1500)
	return x

	# using Nelder Mead simplex method from Optim.jl is 3x slower than newoa
	# res = optimize(x ->self(scratch,x) , xinit, NelderMead())
	# return Optim.minimizer(res)
	
end

function demodulateall( timestamp::Vector{T},data::Matrix{Complex{T}};faintparam::Union{Nothing,FaintStates} = nothing, xinit=[0.01,0],recenter=true,onlyhigh=false)   where{T<:AbstractFloat}

	output = copy(data)
	param = Vector{Modulation{T}}(undef,32) 
	likelihood =  Vector{T}(undef,32) 
	ϕrange= range(0,π,5)

	if !isnothing(faintparam)
		state= buildstates(faintparam, timestamp)
		lag = estimatelag(state,data[:,33])
		@show lag
		state= buildstates(faintparam, timestamp; lag=lag)
	end

	Threads.@threads for (j,k) ∈ collect(Iterators.product(1:4,(FT,SC)))
		FCphase = angle.(data[:,idx(k,j,FC)])
		for i ∈ (D1,D2,D3,D4)
			d = view(data,:,idx(k,j,i)) .*  exp.(-1im.*FCphase)
			if !isnothing(faintparam)
				power = compute_mean_power(state,view(data,:,idx(k,j,i)))
			else
				power = 1.
			end

			if !isnothing(faintparam) && onlyhigh
				hgh =  (state.== HIGH)
				lkl = Chi2CostFunction(timestamp[hgh],d[hgh],power[hgh])
			else
				lkl = Chi2CostFunction(timestamp,d,power)
			end

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

function initialguess(data::Vector{Complex{T}}) where{T<:AbstractFloat}
	std(angle.(data .* exp.(-1ȷ .*angle(mean(data)))))
end

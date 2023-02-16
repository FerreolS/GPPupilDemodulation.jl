
using StatsBase
using OptimPackNextGen
using LinearAlgebra
import OptimPackNextGen.Powell.Newuoa.newuoa

function idx(side_::side,telescope_::Integer, diode_::diode) 
	if diode_==FC
		return 32 + Integer(side_)÷4  + (telescope_-1) + 1
	end
	return Integer(side_) + (Integer(diode_)-1) + (telescope_-1)*4 + 1 
end

mutable struct modulation{T<:AbstractFloat}
	c::Complex{T}
	a::Complex{T}
	b::T
	ϕ1::T
	ω::T
end

function modulation(;T=Float32,c=0.,a=1.,b=1.,ϕ1=0.,ω=2π)
	return modulation{T}(c,a,b,ϕ1, ω)
end

function modulation{T}(mod::modulation) where{T<:AbstractFloat}
	return modulation{T}(convert.(Complex{T},(mod.c,mod.a))...,convert.(T,(mod.b,mod.ϕ1, mod.ω))...)
end

function (self::modulation{T})(timestamp::Vector{T2}) where{T<:AbstractFloat,T2}
	timestamp = T.(timestamp)
	return self.c .+ self.a .* exp.(ȷ .* self.b .* sin.( self.ω .* timestamp .+ self.ϕ1 ))
end

function updatemodulation(self::modulation{T}, timestamp::Vector{T}, data::Vector{Complex{T}},  b::T, ϕ1::T) where{T<:AbstractFloat}
	self.b = b
	self.ϕ1 = ϕ1
	model = exp.(ȷ .* (b .* sin.( self.ω .* timestamp .+ ϕ1 )))
	(self.c, self.a) = linearregression( model, data)
	return 	self.c .+ self.a .* model
end

function updatemodulation!( model::Vector{Complex{T}}, self::modulation{T}, timestamp::Vector{T}, data::Vector{Complex{T}},  b::T, ϕ1::T) where{T<:AbstractFloat}
	self.b = b
	self.ϕ1 = ϕ1
	@. model = exp(ȷ * (b * sin( self.ω * timestamp + ϕ1 )))
	# @inbounds @simd for i in eachindex(model,timestamp)
	# 	model[i] =  exp(ȷ * (b * sin( self.ω * timestamp[i] + ϕ1 )))
	# end
	(self.c, self.a) = linearregression( model, data)
	@. model = 	self.c + self.a * model
	return model
end

function linearregression( model::Vector{Complex{T}}, data::Vector{Complex{T}}) where{T<:AbstractFloat}
	N = length(model)
	Sv1 = sum(data)
	Sv2 = data ⋅ model

	H12 = mean(model)
	detH = 1/(N*(1 - abs2(H12)))
	c = (Sv1 .- Sv2 .* H12) * detH
	a = (-Sv1 .* conj(H12) + Sv2) * detH
	return (c, a)
end


struct Chi2CostFunction{T<:AbstractFloat}
    N::Int
    mod::modulation{T}
    timestamp::Vector{T}
    data::Vector{Complex{T}}
    function Chi2CostFunction{T}(mod::modulation{T},
        timestamp::Vector,
        data::Vector{Complex{T}}) where {T<:AbstractFloat}
        N =length(timestamp);
        @assert N == size(data,1) "voltage and time must have the same number of lines"
        return new{T}(N,mod,convert.(T,timestamp),data)
    end
end

function Chi2CostFunction(timestamp::Vector,data::Vector{Complex{T}}; kwd...) where {T<:AbstractFloat}
	return Chi2CostFunction{T}(modulation(;T=T,kwd...),timestamp,data)
end

function (self::Chi2CostFunction{T})(b::T,ϕ1::T) where{T<:AbstractFloat}
	pupilmodulation = updatemodulation(self.mod, self.timestamp, self.data, b, ϕ1)
	return sum(abs2,( pupilmodulation .- self.data))
end

function Chi2CostFunction(pupilmodulation::Vector{Complex{T}},timestamp::Vector,data::Vector{Complex{T}}; kwd...) where {T<:AbstractFloat}
	return Chi2CostFunction{T}(pupilmodulation,modulation(;T=T,kwd...),timestamp,data)
end

function (self::Chi2CostFunction{T})(pupilmodulation::Vector{Complex{T}},b::T,ϕ1::T) where{T<:AbstractFloat}
	 updatemodulation!(pupilmodulation, self.mod, self.timestamp, self.data, b, ϕ1)
	return sum(abs2,( pupilmodulation .- self.data))
end

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

function demodulateall( time::Vector{T},data::Matrix{Complex{T}})   where{T<:AbstractFloat}

	output = copy(data)
	param = Vector{modulation{T}}(undef,40) 
	Threads.@threads for (j,k) ∈ collect(Iterators.product(1:4,(FT,SC)))
		FCphase = angle.(data[:,idx(k,j,FC)])
		for i ∈ (D1,D2,D3,D4)
			lkl = Chi2CostFunction(time, view(data,:,idx(k,j,i)) .* exp.(-1im.*FCphase))
			minimize!(lkl,xinit=[π/2,0.])
			@. output[:,idx(k,j,i)] = data[:,idx(k,j,i)] * exp(-1im*( angle($(lkl.mod(time)))))
			
			param[idx(k,j,i)] = lkl.mod
		end
	end
	return (output, param)
end

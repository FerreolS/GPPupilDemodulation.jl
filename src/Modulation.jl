
import OptimPackNextGen.Powell.Newuoa.newuoa

using LinearAlgebra
using OptimPackNextGen
using StatsBase
using StaticArrays

@enum Side FT=0 SC=16
@enum Diode D1=1 D2=2 D3=3 D4=4 FC
const M_2PI = 6.283185
const ȷ=im

include("Faint.jl")


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

function (self::Modulation{T})(timestamp::D) where{T<:AbstractFloat,D<:AbstractArray{T}}
	timestamp = T.(timestamp)
	return self.c .+ self.a .* exp.(ȷ .* self.b .* sin.( self.ω .* timestamp .+ self.ϕ ))
end

function getphase(self::Modulation{T},timestamp::D) where{T<:AbstractFloat,D<:AbstractArray{T}}
	timestamp = T.(timestamp)
	return self.b .* sin.( self.ω .* timestamp .+ self.ϕ ) .+ angle(self.a)
end

function updatemodulation(self::Modulation{T}, timestamp::TT, data::D, power::P,  b::T, ϕ::T) where{T<:AbstractFloat,P<:Union{Vector{T}, T},D<:AbstractVector{Complex{T}},TT<:AbstractVector{T}}
	model =Vector{Complex{T}}(undef,length(timestamp))

	return updatemodulation!(self, model,  timestamp, data, power, b, ϕ)
end

function updatemodulation!(self::Modulation{T}, model::Vector{Complex{T}},  timestamp::TT, data::D, power::T, b::T, ϕ::T) where{T<:AbstractFloat,D<:AbstractVector{Complex{T}},TT<:AbstractVector{T}}
	self.b = b
	self.ϕ = ϕ
	if b==0.
		self.c = 0
		self.a = mean(data)
		fill!(model, self.a)
	else
		@. model = exp(ȷ * (b * sin( self.ω * timestamp + ϕ )))
		if power==1.
			(self.c, self.a) = simplelinearregression( model, data)
			@. model = 	self.c + self.a * model
		else
			(self.c, self.a) = simplelinearregression( model, power .\ data)
			(self.c, self.a) = (self.c, self.a).*power
			@. model = 	self.c + self.a * model
		end
	end
	return model
end

function updatemodulation!(self::Modulation{T}, model::Vector{Complex{T}},  timestamp::TT, data::D, power::Vector{T}, b::T, ϕ::T) where{T<:AbstractFloat,D<:AbstractVector{Complex{T}},TT<:AbstractVector{T}}
	self.b = b
	self.ϕ = ϕ
	if b==0.
		self.c = 0
		self.a = sum(map((x,y) -> x^2 * y,power,data)) / sum(x->x^4,power)
		fill!(model,self.a)
	else
		@. model = exp(ȷ * (b * sin( self.ω * timestamp + ϕ )))
		(self.c, self.a) = linearregression( model, data, power)
		@. model = 	self.c + self.a * model
	end
	return model
end


function linearregression( model::Vector{Complex{T}}, data::D, power::Vector{T}) where{T<:AbstractFloat, D<:AbstractVector{Complex{T}}}

	a11 = zero(T)
	a12 = zero(Complex{T})
	a22 = zero(T)
	b1 = zero(Complex{T})
	b2 = zero(Complex{T})
	@inbounds @simd for i in eachindex(model,data,power)
		p2 = power[i]^2
		a11 += p2
		a12 += p2*model[i]
		a22 += p2*abs2(model[i])
		b1 += power[i]*data[i]
		b2 += power[i]*conj(model[i])*data[i]
	end

	A = SMatrix{2,2}(a11, a12, conj(a12), a22)
	b = @SVector [b1,  b2]

	output = A \ b
	
	return tuple(output...) # (c,a)
end


function simplelinearregression( model::Vector{Complex{T}}, data::D) where{T<:AbstractFloat, D<:AbstractVector{Complex{T}}}
	N = length(model)
	Sv1 = sum(data)
	Sv2 = model ⋅ data
	H12 = mean(model)
	detH = 1/(N*(1 - abs2(H12)))
	c = (Sv1 .- Sv2 .* H12) * detH
	a = (-Sv1 .* conj(H12) + Sv2) * detH
	return (c, a)
end


function linearregression( model::Vector{Complex{T}}, data::D) where{T<:AbstractFloat, D<:AbstractVector{Complex{T}}}

	N = length(model)

	# A = @MMatrix zeros(Complex{T},2,2)
    # b = @MVector zeros(Complex{T},2)
	# A[1,1] = N
	# A[2,2] = sum(abs2, model)
	# A[1,2] = sum(model)
	# A[2,1] = conj(A[1,2])

	# b[1] = sum(data)
	# b[2] = model ⋅ data

	a12 = sum(model)
	A = SMatrix{2,2}(N, a12, conj(a12), sum(abs2, model))
	b = @SVector [sum(data) ,  model ⋅ data]
	
	output = A \ b
	
	return tuple(output...) # (c,a)
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
									power::P) where {T<:AbstractFloat,P<:Union{Vector{T}, Real}}
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

myeltype(::Complex{T}) where T = T
myeltype(::AbstractArray{Complex{T}}) where T = T
myeltype(x) = eltype(x)

function weighted_norm2(A::AbstractVector,W::AbstractVector)
	s = zero(promote_type(myeltype(A),myeltype(W)))
	@inbounds @simd for i in eachindex(A,W)
		s += W[i]*abs2(A[i])
	end
	return s
end

weighted_norm2(A::AbstractVector,w::Real) =  w * norm2(A)
function norm2(A::AbstractVector)
	s = zero(myeltype(A))
	@inbounds @simd for i in eachindex(A)
		s += abs2(A[i])
	end
	return s
end


function (self::Chi2CostFunction{T})(b::T,ϕ::T) where{T<:AbstractFloat}
	pupilmodulation = updatemodulation(self.mod, self.timestamp, self.data,self.power, b, ϕ)
	return weighted_norm2( pupilmodulation .- self.data./self.power,self.power.^2)
end

function (self::Chi2CostFunction{T})(pupilmodulation::AbstractVector{Complex{T}},b::T,ϕ::T) where{T<:AbstractFloat}
	updatemodulation!(self.mod, pupilmodulation, self.timestamp, self.data,self.power, b, ϕ)
	return weighted_norm2( pupilmodulation .-  self.data./self.power,self.power.^2)
end

(self::Chi2CostFunction{T})() where{T<:AbstractFloat}  = self(self.mod.b,self.mod.ϕ)
(self::Chi2CostFunction{T})(x::Vector{T}) where{T<:AbstractFloat}  = self(x[1],x[2])
(self::Chi2CostFunction{T})(scratch::Vector{Complex{T}},x::D) where{T<:AbstractFloat,D<:AbstractVector{T}}  = self(scratch,x[1],x[2])

function minimize!(self::Chi2CostFunction{T}; xinit=[2,0]) where {T<:AbstractFloat}
	scratch =Vector{Complex{T}}(undef,self.N)
	xinit = T.(xinit)
	(status, x, χ2) =  newuoa(x ->self(scratch,x) , xinit,1,1e-3; check=false)
	return x

	# using Nelder Mead simplex method from Optim.jl is 3x slower than newoa
	# res = optimize(x ->self(scratch,x) , xinit, NelderMead())
	# return Optim.minimizer(res)
	
end

function demodulateall( timestamp::AbstractVector,data::AbstractMatrix{Complex{T}}; 
							init::Union{Symbol,Vector{T}}=:auto,
							recenter::Bool=true,
							faintparam::Union{Nothing,FaintStates,S} = nothing,
							onlyhigh=false,
							preswitchdelay=0.01,
							postwitchdelay=0.3)  where{T<:AbstractFloat,S<:AbstractVector{MetState}}

	output = copy(data)
	param = Vector{Modulation{T}}(undef,32) 
	likelihood =  Vector{T}(undef,32) 
	ϕrange= range(-π,π,8)

	if !isa(init,Symbol)
		xinit = init
	end

	if isa(faintparam,FaintStates)
		state= buildstates(faintparam, timestamp)
		lag = estimatelag(state,data[:,idx(SC,1,FC)])
		@info "lag = $lag"
		state= buildstates(faintparam, timestamp; lag=lag, preswitchdelay=preswitchdelay,postwitchdelay=postwitchdelay)
	elseif isa(faintparam,AbstractVector{MetState})
		state = faintparam
	end


	Threads.@threads for (j,k) ∈ collect(Iterators.product(1:4,(FT,SC)))
		FCphase = angle.(data[:,idx(k,j,FC)])
		for i ∈ (D1,D2,D3,D4)
			valid = (:)
			if !isnothing(faintparam) 
				if onlyhigh
					valid =  (state.== HIGH)
				else
					valid =  trues(size(FCphase))
				end
				if any(x-> x == TRANSIENT,state) 
					valid .&=  (state .!= TRANSIENT)
				end
				power = compute_mean_power(state,view(data,:,idx(k,j,i)))#[valid]
			else
				power = T.(1.)
			end

			d = view(data,:,idx(k,j,i)) ./ power.*  exp.(-1im.*FCphase)
			

			lkl = Chi2CostFunction(timestamp[valid],d[valid],1,ω=M_2PI)

			if init==:auto
				binit= 0.1#initialguess(d)
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
			modparam = lkl.mod
			if modparam.b<0
				modparam.b *= -1
				modparam.ϕ +=ifelse(modparam.ϕ<0,+π,-π)
			end
			param[idx(k,j,i)] = modparam
		end
	end
	return (output, param,likelihood)
end

function initialguess(data::D) where{T<:AbstractFloat,D<:AbstractVector{Complex{T}}}
	std(angle.(data .* exp.(-1ȷ .*angle(mean(data)))))
end

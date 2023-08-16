
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

abstract type Modulation{T}  end

mutable struct ModulationWithOffsets{T<:AbstractFloat} <: Modulation{T} 
	c::Complex{T}
	a::Complex{T}
	b::T
	ϕ::T
	ω::T
end

mutable struct ModulationNoOffsets{T<:AbstractFloat} <: Modulation{T} 
	a::Complex{T}
	b::T
	ϕ::T
	ω::T
end

function ModulationWithOffsets(;T=Float32,c=0.,a=1.,b=1.,ϕ=0.,ω=2π)
	return ModulationWithOffsets{T}(c,a,b,ϕ, ω)
end

function ModulationWithOffsets{T}(mod::ModulationWithOffsets) where{T<:AbstractFloat}
	return ModulationWithOffsets{T}(convert.(Complex{T},(mod.c,mod.a))...,convert.(T,(mod.b,mod.ϕ, mod.ω))...)
end

function ModulationNoOffsets(;T=Float32,a=1.,b=1.,ϕ=0.,ω=2π)
	return ModulationNoOffsets{T}(a,b,ϕ, ω)
end

function ModulationNoOffsets{T}(mod::ModulationWithOffsets) where{T<:AbstractFloat}
	return ModulationNoOffsets{T}(convert.(Complex{T},mod.a),convert.(T,(mod.b,mod.ϕ, mod.ω))...)
end

function (self::Modulation{T})(timestamp::D) where{T<:AbstractFloat,D<:AbstractArray{T}}
	timestamp = T.(timestamp)
	if isa(self,ModulationWithOffsets{T})
		return self.c .+ self.a .* exp.(ȷ .* self.b .* sin.( self.ω .* timestamp .+ self.ϕ ))
	else
		return self.a .* exp.(ȷ .* self.b .* sin.( self.ω .* timestamp .+ self.ϕ ))
	end
end

function getphase(self::Modulation{T},timestamp::D) where{T<:AbstractFloat,D<:AbstractArray{T}}
	timestamp = T.(timestamp)
	return self.b .* sin.( self.ω .* timestamp .+ self.ϕ ) .+ angle(self.a)
end

function updatemodulation(self::Modulation{T}, timestamp::TT, data::D, weight::P,  b::T, ϕ::T) where{T<:AbstractFloat,P<:Union{Vector{T}, T},D<:AbstractVector{Complex{T}},TT<:AbstractVector{T}}
	model =Vector{Complex{T}}(undef,length(timestamp))

	return updatemodulation!(self, model,  timestamp, data, weight, b, ϕ)
end

function updatemodulation(	self::Modulation{T}, 
							timestamp::TT, 
							data::D, 
							weight::P,
							power::Q,  
							b::T, 
							ϕ::T) where{T<:AbstractFloat,
										P<:Union{AbstractVector{T}, T},
										D<:AbstractVector{Complex{T}},
										Q<:AbstractVector{Complex{T}},
										TT<:AbstractVector{T}}
	model =Vector{Complex{T}}(undef,length(timestamp))

	return updatemodulation!(self, model,  timestamp, data, weight,power, b, ϕ)
end

function updatemodulation!(	self::M,
							model::Vector{Complex{T}},  
							timestamp::TT, 
							data::D, 
							::Real, 
							b::T,
							ϕ::T) where{T<:AbstractFloat,D<:AbstractVector{Complex{T}},TT<:AbstractVector{T},M<:Modulation{T}}
	self.b = b
	self.ϕ = ϕ
	if b==0.
		if M == ModulationWithOffsets{T}
			self.c = 0
		end
		self.a = mean(data)
		fill!(model, self.a)
	else
		@. model = exp(ȷ * (b * sin( self.ω * timestamp + ϕ )))

		if M == ModulationWithOffsets{T}
			(self.c, self.a) = simplelinearregression( model, data)
			@. model = 	self.c + self.a * model
		else
			self.a = (model ⋅ data) / sum(abs2,model )
			@. model = 	self.a * model
		end
	end
	return model
end

function updatemodulation!(	self::M,
							model::Vector{Complex{T}},  
							timestamp::TT, 
							data::D, 
							weight::P, 
							power::Q,  
							b::T,
							ϕ::T) where{T<:AbstractFloat,
										D<:AbstractVector{Complex{T}},
										TT<:AbstractVector{T},
										P<:AbstractVector{T},
										Q<:AbstractVector{Complex{T}},
										M<:Modulation{T}}
	self.b = b
	self.ϕ = ϕ
	@. model = power * exp(ȷ * b * sin( self.ω * timestamp + ϕ ))
	
	if M == ModulationWithOffsets{T}
		(self.c, self.a) = linearregression( model, data, weight)
		@. model = 	self.c + self.a * model
	else
		mw = model .* weight
		self.a = (mw ⋅ data) / (mw ⋅ model)
		@. model = 	self.a * model
	end
	return model
end

# function updatemodulation!( self::ModulationWithOffsets{T}, 
# 							model::Vector{Complex{T}},  
# 							timestamp::TT, 
# 							data::D, 
# 							weight::Vector{T}, 
# 							b::T, 
# 							ϕ::T) where{T<:AbstractFloat,
# 										D<:AbstractVector{Complex{T}},
# 										TT<:AbstractVector{T}}
# 	self.b = b
# 	self.ϕ = ϕ
# 	if b==0.
# 		self.c = 0
# 		self.a = sum( weight .* data) / sum(weight)
# 		fill!(model,self.a)
# 	else
# 		@. model = exp(ȷ * (b * sin( self.ω * timestamp + ϕ )))
# 		(self.c, self.a) = linearregression( model, data, weight)
# 		@. model = 	self.c + self.a * model
# 	end
# 	return model
# end


function linearregression( model::Vector{Complex{T}}, data::D, weight::Vector{T}) where{T<:AbstractFloat, D<:AbstractVector{Complex{T}}}

	a11 = zero(T)
	a12 = zero(Complex{T})
	a22 = zero(T)
	b1 = zero(Complex{T})
	b2 = zero(Complex{T})
	@inbounds @simd for i in eachindex(model,data,weight)
		a11 += weight[i]
		a12 += weight[i]*model[i]
		a22 += weight[i]*abs2(model[i])
		b1 	+= weight[i]*data[i]
		b2	+= weight[i]*conj(model[i])*data[i]
	end

	A = SMatrix{2,2}([ a11 a12; conj(a12) a22])
	b = @SVector [b1,  b2]

	output = A \ b
	
	return tuple(output...) # (c,a)
end

function linearregression( model::Vector{Complex{T}}, data::D, ::Real) where{T<:AbstractFloat, D<:AbstractVector{Complex{T}}}
	N = T(length(model))
	a12 = zero(Complex{T})
	a22 = zero(Complex{T})
	b1 = zero(Complex{T})
	b2 = zero(Complex{T})
	@inbounds @simd for i in eachindex(model,data)
		a12 += model[i]
		a22 += abs2(model[i])
		b1 	+= data[i]
		b2	+= conj(model[i])*data[i]
	end
	A = SMatrix{2,2}([ N a12; conj(a12) a22])
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

	a12 = sum(model)
	A = SMatrix{2,2}(N, a12, conj(a12), sum(abs2, model))
	b = @SVector [sum(data) ,  model ⋅ data]
	
	output = A \ b
	
	return tuple(output...) # (c,a)
end

Base.collect(A::Vector{Float64})=A

struct Chi2CostFunction{T<:AbstractFloat,P<:Union{Vector{T}, T},M<:Modulation}
    N::Int64
    mod::M
    timestamp::Vector{T}
    data::Vector{Complex{T}}
    weight::P
	power::Vector{Complex{T}}
    function Chi2CostFunction{T,P,M}(mod::M,
        							timestamp::Vector{T},
        							data::AbstractVector{Complex{T}},
									weight::P,
									power::AbstractVector{Complex{T}}) where {T<:AbstractFloat,P<:Union{AbstractVector{T}, Real},M<:Modulation}
        N =length(timestamp);
        @assert N == size(data,1) "voltage and time must have the same number of lines"
		if P <: Vector
			@assert N == size(weight,1) "weight and time must have the same number of lines"
			return new{T,Vector{T},M}(N,mod,timestamp,collect(data),collect(weight),collect(power))
		end
		@assert N == size(power,1) "power and time must have the same number of lines"
        return new{T,T,M}(N,mod,timestamp,collect(data),weight,collect(power))
    end
end

function Chi2CostFunction(timestamp::AbstractVector,data::AbstractVector{Complex{T}},power::AbstractVector; offsets=false, kwd...) where {T<:AbstractFloat}
	if offsets
		mod = ModulationWithOffsets(;T=T,kwd...)
	else
		mod = ModulationNoOffsets(;T=T,kwd...)
	end
	return Chi2CostFunction{T,T,typeof(mod)}(mod,T.(timestamp),data,T.(1.0),power)
end

function Chi2CostFunction(timestamp::AbstractVector,data::AbstractVector{Complex{T}},weight::AbstractVector,power::AbstractVector;offsets=false, kwd...) where {T<:AbstractFloat}
	if offsets
		mod = ModulationWithOffsets(;T=T,kwd...)
	else
		mod = ModulationNoOffsets(;T=T,kwd...)
	end
	return Chi2CostFunction{T,Vector{T},typeof(mod)}(mod,T.(timestamp),data,T.(weight),Complex{T}.(power))
end

function Chi2CostFunction(timestamp::AbstractVector,data::AbstractVector{Complex{T}},weight::Number,power::AbstractVector; offsets=false, kwd...) where {T<:AbstractFloat}
	if offsets
		mod = ModulationWithOffsets(;T=T,kwd...)
	else
		mod = ModulationNoOffsets(;T=T,kwd...)
	end
	return Chi2CostFunction{T,T,typeof(mod)}(mod,T.(timestamp),data,T.(weight),Complex{T}.(power))
end

myeltype(::Complex{T}) where T = T
myeltype(::AbstractArray{Complex{T}}) where T = T
myeltype(x) = eltype(x)

function weighted_norm2(A::AbstractVector,weight::AbstractVector)
	s = zero(promote_type(myeltype(A),myeltype(weight)))
	@inbounds @simd for i in eachindex(A,weight)
		s += weight[i]*abs2(A[i])
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
	pupilmodulation = updatemodulation(self.mod, self.timestamp, self.data,self.weight,self.power, b, ϕ)
	return weighted_norm2( pupilmodulation .- self.data,self.weight)./length(pupilmodulation)
end

function (self::Chi2CostFunction{T})(pupilmodulation::AbstractVector{Complex{T}},b::T,ϕ::T) where{T<:AbstractFloat}
	updatemodulation!(self.mod, pupilmodulation, self.timestamp, self.data,self.weight,self.power, b, ϕ)
	return weighted_norm2( pupilmodulation .-  self.data,self.weight)./length(pupilmodulation)
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
							fitoffsets=false,
							preswitchdelay=0.01,
							postwitchdelay=0.3)  where{T<:AbstractFloat,S<:AbstractVector{MetState}}

	output = copy(data)
	if fitoffsets
		param = Vector{ModulationWithOffsets{T}}(undef,32) 
	else
		param = Vector{ModulationNoOffsets{T}}(undef,32) 
	end
	likelihood =  Vector{T}(undef,32) 
	ϕrange= range(-π,π,8)

	if !isa(init,Symbol)
		xinit = init
	end

	if isa(faintparam,FaintStates)
		state= buildstates(faintparam, timestamp; preswitchdelay=preswitchdelay,postwitchdelay=postwitchdelay)
	elseif isa(faintparam,AbstractVector{MetState})
		state = faintparam
	end


	valid = (:)
	if !isnothing(faintparam) 
		if onlyhigh
			valid = (state.== HIGH) .|| (state.== NORMAL)
		else
			valid =  trues(size(state))
		end
		if any(x-> x == TRANSIENT,state) 
			valid .&=  (state .!= TRANSIENT)
		end
	else
		weight = power = T.(1.)
	end

	Threads.@threads for (j,k) ∈ collect(Iterators.product(1:4,(FT,SC)))
		FCphasor = exp.(1im.*angle.(data[:,idx(k,j,FC)]))
		for i ∈ (D1,D2,D3,D4)
			d = view(data,:,idx(k,j,i)) 
			if !isnothing(faintparam) 
				(power,weight) = compute_mean_var_power(state[valid],d[valid])
			else
				weight = power = T.(1.)
			end
			p  = power.* FCphasor[valid]
			

			lkl = Chi2CostFunction(timestamp[valid],d[valid] ,weight,p,ω=M_2PI,offsets=fitoffsets)


			if init==:auto
				binit= 0.1#initialguess(d)
				ϕinit = ϕrange[argmin(map(ϕ -> lkl(binit,ϕ),ϕrange ))]
				xinit=[binit, ϕinit]
			end
			x = minimize!(lkl,xinit=xinit)
			lklval =  lkl(x)
			ϕπ = x[2]+ifelse(x[2]<0,+π,-π)
			
			if lklval > lkl(x[1],ϕπ)
				@info "bad minima at $k, $j, $i : $lklval > $(lkl(x[1],ϕπ))"
				x = minimize!(lkl,xinit=[x[1],ϕπ])
			end

			likelihood[idx(k,j,i)] = lkl(x)
			if recenter
				if fitoffsets
					@. output[:,idx(k,j,i)] = (d  - lkl.mod.c) * exp(-1im*( $(getphase(lkl.mod, timestamp)) - angle(lkl.mod.a)))
				else
					@. output[:,idx(k,j,i)] = (d  ) * exp(-1im*( $(getphase(lkl.mod, timestamp)) - angle(lkl.mod.a)))
				end
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

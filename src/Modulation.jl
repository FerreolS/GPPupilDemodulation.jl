
import OptimPackNextGen.Powell.Newuoa.newuoa

using LinearAlgebra
using Optim
using OptimPackNextGen
using StatsBase
using Zygote

@enum side FT=0 SC=16
@enum diode D1=1 D2=2 D3=3 D4=4 FC
const ȷ=im


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
	ϕ::T
	ω::T
end

function modulation(;T=Float32,c=0.,a=1.,b=1.,ϕ=0.,ω=2π)
	return modulation{T}(c,a,b,ϕ, ω)
end

function modulation{T}(mod::modulation) where{T<:AbstractFloat}
	return modulation{T}(convert.(Complex{T},(mod.c,mod.a))...,convert.(T,(mod.b,mod.ϕ, mod.ω))...)
end

function (self::modulation{T})(timestamp::Vector{T2}) where{T<:AbstractFloat,T2}
	timestamp = T.(timestamp)
	return self.c .+ self.a .* exp.(ȷ .* self.b .* sin.( self.ω .* timestamp .+ self.ϕ ))
end

function getphase(self::modulation{T},timestamp::Vector{T2}) where{T<:AbstractFloat,T2}
	timestamp = T.(timestamp)
	return self.b .* sin.( self.ω .* timestamp .+ self.ϕ ) .+ angle(self.a)
end

function updatemodulation(self::modulation{T}, timestamp::Vector{T}, data::Vector{Complex{T}},  b::T, ϕ::T) where{T<:AbstractFloat}
	self.b = b
	self.ϕ = ϕ
	if b==0.
		self.c = 0
		self.a = mean(data)
		return self.a .* ones(Complex{T},size(data))
	end
	model = exp.(ȷ .* (b .* sin.( self.ω .* timestamp .+ ϕ )))
	(self.c, self.a) = Zygote.@ignore linearregression( model, data)
	return 	self.c .+ self.a .* model
end

function updatemodulation!( model::Vector{Complex{T}}, self::modulation{T}, timestamp::Vector{T}, data::Vector{Complex{T}},  b::T, ϕ::T) where{T<:AbstractFloat}
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



function updatemodulation!( model::Vector{Complex{T}}, jac::Matrix{Complex{T}}, self::modulation{T}, timestamp::Vector{T}, data::Vector{Complex{T}},  b::T, ϕ::T) where{T<:AbstractFloat}
	self.b = b
	self.ϕ = ϕ
	@. model = exp(ȷ * (b * sin( self.ω * timestamp + ϕ )))

	@. jac[:,1] = ȷ * model * sin( self.ω * timestamp + ϕ )
	# D model / D ϕ
	@. jac[:,2] = ȷ * b* model * cos( self.ω * timestamp + ϕ )

	(self.c, self.a) = linearregression( model, data)
	@. model = 	self.c + self.a * model
	@. jac = self.a * model
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

function (self::Chi2CostFunction{T})(b::T,ϕ::T) where{T<:AbstractFloat}
	pupilmodulation = updatemodulation(self.mod, self.timestamp, self.data, b, ϕ)
	return sum(abs2,( pupilmodulation .- self.data))
end

function Chi2CostFunction(pupilmodulation::Vector{Complex{T}},timestamp::Vector,data::Vector{Complex{T}}; kwd...) where {T<:AbstractFloat}
	return Chi2CostFunction{T}(pupilmodulation,modulation(;T=T,kwd...),timestamp,data)
end

function (self::Chi2CostFunction{T})(pupilmodulation::Vector{Complex{T}},b::T,ϕ::T) where{T<:AbstractFloat}
	updatemodulation!(pupilmodulation, self.mod, self.timestamp, self.data, b, ϕ)
	return sum(abs2,( pupilmodulation .- self.data))
end

function (self::Chi2CostFunction{T})(pupilmodulation::Vector{Complex{T}},b::T,ϕ::T,gx::Vector{T}) where{T<:AbstractFloat}
	jac = similar(pupilmodulation, (axes(pupilmodulation)...,2))
	updatemodulation!(pupilmodulation, jac, self.mod, self.timestamp, self.data, b, ϕ)
   	residual = pupilmodulation .- self.data
	gx .= real.(jac' * residual)
	return 0.5 * sum(abs2,residual)
end

(self::Chi2CostFunction{T})() where{T<:AbstractFloat}  = self(self.mod.b,self.mod.ϕ)
(self::Chi2CostFunction{T})(x::Vector{T}) where{T<:AbstractFloat}  = self(x[1],x[2])
(self::Chi2CostFunction{T})(scratch::Vector{Complex{T}},x::Vector{T}) where{T<:AbstractFloat}  = self(scratch,x[1],x[2])

(self::Chi2CostFunction{T})(scratch::Vector{Complex{T}},x::Vector{T},gx::Vector{T}) where{T<:AbstractFloat}  = self(scratch,x[1],x[2],gx)


function minimize!(self::Chi2CostFunction{T};  kwd...) where {T<:AbstractFloat}
	minimize!(Val(:newoa),self::Chi2CostFunction{T};kwd...) 
end

function minimize!(::Val{:newoa},self::Chi2CostFunction{T};  xinit=[0.01,0]) where {T<:AbstractFloat}
	scratch =Vector{Complex{T}}(undef,self.N)
	xinit = T.(xinit)
	(status, x, χ2) =  newuoa(x ->self(scratch,x) , xinit,1,1e-3,maxeval=1500)
	return x	
end

function minimize!(::Val{:simplex},self::Chi2CostFunction{T};  xinit=[0.01,0]) where {T<:AbstractFloat}
	scratch =Vector{Complex{T}}(undef,self.N)
	xinit = T.(xinit)
	# using Nelder Mead simplex method from Optim.jl is 3x slower than newoa
	res = optimize(x ->self(scratch,x) , xinit, NelderMead())
	return Optim.minimizer(res)
	
end

function minimize!(::Val{:vmlmbZ},self::Chi2CostFunction{T};  xinit=[0.01,0]) where {T<:AbstractFloat}
	xinit = T.(xinit)
	x = vmlmb(x->self(x), xinit;autodiff=true)
	return x
end

function minimize!(::Val{:vmlmb},self::Chi2CostFunction{T}; xinit=[2,0]) where {T<:AbstractFloat}	
	scratch =Vector{Complex{T}}(undef,self.N)
	xinit = T.(xinit)
	gx = similar(xinit)
	cost(x,gx) = self(scratch,x,gx)
	x = vmlmb(cost, xinit)
	return x
end

function demodulateall( time::Vector{T},data::Matrix{Complex{T}}; algo=:newoa, xinit=[0.01,0],recenter=true)   where{T<:AbstractFloat}

	output = copy(data)
	param = Vector{modulation{T}}(undef,32) 
	likelihood =  Vector{T}(undef,32) 
	ϕrange= range(0,π,5)
	Threads.@threads for (j,k) ∈ collect(Iterators.product(1:4,(FT,SC)))
		FCphase = angle.(data[:,idx(k,j,FC)])
		for i ∈ (D1,D2,D3,D4)
			d = view(data,:,idx(k,j,i)) .*  exp.(-1im.*FCphase)
			lkl = Chi2CostFunction(time,d )
			if xinit==:auto
				binit= initialguess(d)
				ϕinit = ϕrange[argmin(map(ϕ -> lkl(binit,ϕ),ϕrange ))]
				xinit=[binit, ϕinit]
			end
			x = minimize!(Val(algo),lkl,xinit=xinit)
			likelihood[idx(k,j,i)] = lkl(x)
			if recenter
				@. output[:,idx(k,j,i)] = (d  - lkl.mod.c) * exp(-1im*( $(getphase(lkl.mod, time)) - FCphase- angle(lkl.mod.a)))
			else
				@. output[:,idx(k,j,i)] = data[:,idx(k,j,i)] * exp(-1im*( angle($(lkl.mod(time)))))
			end
			param[idx(k,j,i)] = lkl.mod
		end
	end
	return (output, param,likelihood)
end

function initialguess(data::Vector{Complex{T}}) where{T<:AbstractFloat}
	std(angle.(data .* exp.(-1ȷ .*angle(mean(data)))))
end

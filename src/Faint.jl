@enum MetState OFF=0 LOW=1 NORMAL=2 HIGH=3 TRANSIENT=-1
 
struct FaintStates{T<:AbstractFloat,A<:AbstractVector{T}}
	timer1::A
	timer2::A
	voltage1::T
	voltage2::T
	state1::MetState
	state2::MetState
end
 
function  FaintStates(state1::V,state2::V,voltage1,voltage2) where {T<:AbstractFloat, V<:AbstractVector{T}}
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


	premax = ceil( Int, preswitchdelay / timestep)
	postmax = ceil( Int, postwitchdelay / timestep)

	currentstate = NORMAL
	first1 = popfirst!(t1)
	first2 = popfirst!(t2)

	forget = 0
	@inbounds @simd for index ∈ 1:N
		time = timestamp[index]
		# HIGH
		if time >= first1  
			currentstate = faintstates.state1
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

function estimatelag(states::S ,data::D; range::R=-10:10) where {T<:AbstractFloat,D<:AbstractVector{Complex{T}},R<:AbstractVector{Int64},S<:AbstractVector{MetState}}
	m = [mean(abs,data[circshift(states,i) .== HIGH]) for i ∈ range];
	return range[argmax(m)]
end

function compute_mean_power(states::S ,data::D) where {T<:AbstractFloat,D<:AbstractVector{Complex{T}},S<:AbstractVector{MetState}}
	pow = zeros(T,length(data))
	 for st ∈ instances(MetState)
		idx = states.==st
		pow[idx] .= mean(abs,data[idx]);
	 end
	 return pow
end

function compute_mean_var_power(states::S ,data::D) where {T<:AbstractFloat,D<:AbstractVector{Complex{T}},S<:AbstractVector{MetState}}
	m = zeros(T,length(data))
	w = zeros(T,length(data))

	for st ∈ instances(MetState)
		idx = states.==st
		_m = mean(abs,data[idx])
		m[idx] .= _m
		w[idx] .= 1 ./ var(abs.(data[idx]);mean=_m)
	end
	return (m,w)
end
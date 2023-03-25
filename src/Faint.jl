
@enum MetState OFF=0 LOW=1 NORMAL=2 HIGH=3


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
	if  voltage1 > voltage2
	 	return FaintStates{T,A}(state1,state2,voltage1,voltage2,LOW,HIGH)
	end
	return FaintStates{T,A}(state1,state2,voltage1,voltage2,HIGH, LOW)

end

function buildstates(faintstates::FaintStates{T,A},timestamp::Vector{T}; lag::Integer=0) where {T<:AbstractFloat,A<:AbstractVector{T}}
	N = length(timestamp)
	laststate = NORMAL
	states = Vector{MetState}(undef,N)
	timestep = timestamp[2]-timestamp[1]
	t1 = collect(faintstates.timer1 .+ lag*timestep)
	t2 = collect(faintstates.timer2 .+ lag*timestep)
	first1 = popfirst!(t1)
	first2 = popfirst!(t2)

	@inbounds @simd for index ∈ 1:N
		time = timestamp[index]
		if time >= first1
			if isempty(t1) 
				first1 = last(timestamp)
				laststate = faintstates.state1
				if first2 == last(timestamp)
					laststate = NORMAL
				end
			else
				first1 = popfirst!(t1)
				laststate = faintstates.state1
			end
		end

		if time >= first2
			if isempty(t2) 
				first2 = last(timestamp)
				laststate = faintstates.state2
				if first1 == last(timestamp)
					laststate = NORMAL
				end
			else
				first2 = popfirst!(t2)
				laststate = faintstates.state2
			end
		end

		states[index] = laststate
	end	
	return states
end

function estimatelag(states::Vector{MetState} ,data::Vector{Complex{T}}; range::AbstractVector{Int64}=-10:10) where {T<:AbstractFloat}
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
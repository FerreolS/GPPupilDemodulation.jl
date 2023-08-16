function getmetrology(filename)
	getmetrology(FITS(filename));
end

function getmetrology(f::FITS)
	hdr = read_header(f[1])
	metrologyhdu = f["METROLOGY"];
	table = Dict(metrologyhdu);
	volt = Float64.(table["VOLT"])
	cmplxV = volt[1:2:end,:]' .+ im*volt[2:2:end,:]';

	times = (table["TIME"]).*1e-6 .+ (DAY_TO_SEC * hdr["MJD-OBS"] ) ;

	return (cmplxV,times,table)
end


function bright2states(bright::Vector{Int32})
	local maping::Vector{MetState} = [NORMAL,LOW,HIGH,OFF,OFF,OFF,OFF,OFF,OFF,OFF,TRANSIENT]
	N = length(bright)
	states = Vector{MetState}(undef,N)
	@inbounds @simd for i in eachindex(states, bright)
		states[i] = maping[bright[i]+1]
	end
	return states
end

struct avgVolts
    VX::Float64
    eVX::Float64
    VY::Float64
    eVY::Float64
end

function read_avg_v_values(filename::AbstractString)::Dict{String, avgVolts}
    avg_dict = Dict{String, avgVolts}()
    open(filename) do file
        for line in eachline(file)
            if occursin(r"^avg", line)
                values = split(line)
                key = values[2]
                float_values = avgVolts([1e-3 .* parse(Float64, values[i]) for i in 3:6]...)
                if haskey(avg_dict, key)
                    push!(avg_dict[key], float_values)
                else
                    avg_dict[key] = float_values
                end
            end
        end
    end
    avg_dict
end


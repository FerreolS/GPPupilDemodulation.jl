function getmetrology(filename)
	f =  FITS(filename);
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
		states[i] = maping[bright[i]]
	end
	return states
end
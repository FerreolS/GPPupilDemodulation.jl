module GPPupilDemodulation

using Logging, Pkg
include("FitsUtils.jl")

import .FitsUtils: FITScopy!

include("Modulation.jl")



using ArgParse, FITSIO

const SUFFIXES = [".fits", ".fits.gz","fits.Z"]
const MJD_1970_1_1 = 40587.0
const DAY_TO_SEC = 24*60*60 

include("Utils.jl")


const units = Dict("TIME"=>"usec","VOLT"  => "V", "POWER_LASER" => "mV", "LAMBDA_LASER" => "m", "FLAG"=>"-")
"""
	endswith(chain::Union{String,Vector{String}}, pattern::Vector{String})

Overloading of Base.endswith for vectors of string. Return `true` if `chain`
ends with one of the patterns given in `pattern`
"""
function Base.endswith(chain::Union{String,Vector{String}}, pattern::Vector{String})
	for str in pattern
		if endswith(chain,str)
			return true
		end
	end
	return false
end

"""
	endswith(chain::Vector{String}, pattern::String)

Overloading of Base.endswith for vectors of string. Return `true` if `chain`
ends with `pattern`
"""
function Base.endswith(chains::Vector{String},pattern::AbstractString)
	for chain in chains
		if endswith(chain,pattern)
			return true
		end
	end
	return false
end

"""
    buildfaintparameters(hdr::FITSHeader)

Extracts data from the FITS header `hdr` to build an instance of the `FaintStates` struct. 

Parameters:
- `hdr`: A FITS header object containing the necessary data.

Returns:
- An instance of the `FaintStates` struct initialized with the extracted data.
"""

function buildfaintparameters(hdr::FITSHeader)

	mjdobs = hdr["MJD-OBS"] # ESO PCR ACQ START (microsecond precision)
	rate1 = hdr["ESO INS ANLO3 RATE1"]
	rate2 = hdr["ESO INS ANLO3 RATE2"]

	repeat1 = hdr["ESO INS ANLO3 REPEAT1"]
	repeat2 = hdr["ESO INS ANLO3 REPEAT2"]

	start1 = hdr["ESO INS ANLO3 TIMER1"] +  MJD_1970_1_1 * DAY_TO_SEC
	start2 = hdr["ESO INS ANLO3 TIMER2"] +  MJD_1970_1_1 * DAY_TO_SEC

	voltage1 = hdr["ESO INS ANLO3 VOLTAGE1"]
	voltage2 = hdr["ESO INS ANLO3 VOLTAGE2"]
	timer1 = start1 .+ rate1 .* (0:(repeat1-1))
	timer2 = start2 .+ rate2 .* (0:(repeat2-1))
	return FaintStates(timer1,timer2,voltage1,voltage2)
end


function read_stefan_file()
	return read_stefan_file(joinpath(pkgdir(@__MODULE__),"data/Stefan_file.txt"))
end

function read_stefan_file(filename::AbstractString)
	offsets = Vector{ComplexF64}(undef,40)
    open(filename) do file
        for line in eachline(file)
            if occursin(r"^avg", line)
                values = split(line)
                name = values[2]
				side = @eval $(Symbol(name[1:2]))
				telescope = parse(Int,name[4])
				diode = @eval $(Symbol(name[5:6]))
				offsets[idx(side,telescope, diode)] = 1e-3 .* (parse(Float64, values[3]) + 1im * parse(Float64, values[5])) 
                
            end
        end
    end
    return offsets
end
function compute_offsets(cmplxV::AbstractMatrix{Complex{T}},state::Vector{MetState}) where{T<:AbstractFloat}
	offsets = zeros(ComplexF64,40)
	Threads.@threads for (i,j,k) ∈ collect(Iterators.product((D1,D2,D3,D4,FC),1:4,(FT,SC)))
		res = fit(Circle, reim(cmplxV[state .== HIGH,idx(k,j,i)])...)
		(x,y,r) = coef(res) 
		offsets[idx(k,j,i),1] = complex(x , y)
		#offsets[idx(k,j,i),2] = r
	end
	return offsets
end 


function compute_offsets(cmplxV::AbstractMatrix{Complex{T}},::Nothing) where{T<:AbstractFloat}
	offsets = zeros(ComplexF64,40)
	Threads.@threads for (i,j,k) ∈ collect(Iterators.product((D1,D2,D3,D4,FC),1:4,(FT,SC)))
		res = fit(Circle, reim(cmplxV[:,idx(k,j,i)])...)
		(x,y,r) = coef(res) 
		offsets[idx(k,j,i)] = complex(x , y)
	end
	return offsets
end 


function processmetrology(metrologyhdu::TableHDU, mjd::Float64; 
								window  = nothing,
								faintparam::Union{Nothing,FaintStates} = nothing, 
								keepraw = false,
								verb=false,
								onlyhigh=false,
								offsets::Union{Vector{ComplexF64},Bool}=true)


	hdr = read_header(metrologyhdu)
	table = Dict(metrologyhdu)
	times = Float64.(table["TIME"]).*1e-6 .+ (DAY_TO_SEC * mjd  )

	if isnothing(faintparam)
		state = nothing
	else
		state = buildstates(faintparam, times );
	end

	volt = Float64.(table["VOLT"])
	cmplxV = volt[1:2:end,:]' .+ im*volt[2:2:end,:]'
	
	fitoffsets = false
	if isa( offsets, Vector{ComplexF64})
		cmplxV .-= reshape(offsets,1,40)
	elseif offsets === true
		cmplxV .-= reshape(compute_offsets(cmplxV,state),1,40)
	else
		fitoffsets = true
	end


	if isnothing(window)
		(output, param,likelihood) = demodulateall(times, cmplxV; faintparam = state,onlyhigh=onlyhigh, fitoffsets=fitoffsets)
		
		if keepraw
			s = similar(volt,80+64,size(volt,2))
			s[1:80,:] .= volt
			s[81:2:end,:] .=  real(output[:,1:32])'
			s[82:2:end,:] .=  imag(output[:,1:32])'
			volt = s
		else
			volt[1:2:end,:] .=  real(output)'
			volt[2:2:end,:] .=  imag(output)'
		end
		
		for (k,j,i) ∈ Iterators.product((D1,D2,D3,D4),1:4,(FT,SC)) 
			b = param[idx(i,j,k)].b
			ϕ = param[idx(i,j,k)].ϕ
			if (b<0)
				b = -b
				ϕ = rem2pi(ϕ+π,RoundNearest) 
			end
			if fitoffsets
				setindex!(hdr,real(param[idx(i,j,k)].c),"DEMODULATION CENTER X0 $i T$j $k")
				setindex!(hdr,imag(param[idx(i,j,k)].c),"DEMODULATION CENTER Y0 $i T$j $k")
			end
			setindex!(hdr,abs(param[idx(i,j,k)].a),"DEMODULATION AMPLITUDE ABS $i T$j $k")
			setindex!(hdr,angle(param[idx(i,j,k)].a),"DEMODULATION AMPLITUDE ARG $i T$j $k")
			setindex!(hdr,b,"DEMODULATION SIN AMPLITUDE $i T$j $k")
			setindex!(hdr,ϕ,"DEMODULATION SIN PHASE $i T$j $k")
		end
		
	else
		nwindow = round(Int, window / (times[2] - times[1]))
		output = similar(cmplxV,length(times),40)

		if fitoffsets
			x0 = similar(times,32,length(times))
			y0 = similar(times,32,length(times))
		end
		absa = similar(times,32,length(times))
		arga = similar(times,32,length(times))
		b = similar(times,32,length(times))
		ϕ = similar(times,32,length(times))

		@views for I in Iterators.partition(axes(times, 1), nwindow)
			(_output, param,likelihood) = demodulateall( times[I], cmplxV[I,:]; faintparam = isnothing(state) ? state : state[I], onlyhigh=onlyhigh, fitoffsets=fitoffsets)
	
			output[I,:] .= _output
			
			for (k,j,i) ∈ Iterators.product((D1,D2,D3,D4),1:4,(FT,SC)) 
				tb = param[idx(i,j,k)].b
				tϕ = param[idx(i,j,k)].ϕ
				if (tb<0)
					tb = -tb
					tϕ = rem2pi(tϕ+π,RoundNearest) 
				end
				b[idx(i,j,k),I] .= tb
				ϕ[idx(i,j,k),I] .= tϕ 
				if fitoffsets
					x0[idx(i,j,k),I] .= real(param[idx(i,j,k)].c)
					y0[idx(i,j,k),I] .= imag(param[idx(i,j,k)].c)
				end
				absa[idx(i,j,k),I] .= abs(param[idx(i,j,k)].a)
				arga[idx(i,j,k),I] .= angle(param[idx(i,j,k)].a)
			end
		end
		#ioutput = reshape(ioutput,:,40)

		if keepraw
			s = similar(volt,80+64,size(volt,2))
			s[1:80,:] .= volt
			s[81:2:end,:] .=  real(output[:,1:32])'
			s[82:2:end,:] .=  imag(output[:,1:32])'
			volt = s
		else
			volt[1:2:end,:] .=  real(output)'
			volt[2:2:end,:] .=  imag(output)'
		end

		if fitoffsets
			table["X0"] = Float32.(x0)
			table["Y0"] = Float32.(y0)
		end
		table["ABSA"] = Float32.(absa)
		table["ARGA"] = Float32.(arga)
		table["B"] = Float32.(b)
		table["PHI"] = Float32.(ϕ)
		if !isnothing(state) 
			table["STATE"] = Int8.(state)
		end

	end
	setindex!(hdr,"GPPupilDemodulation.jl","PROCSOFT")
	table["VOLT"] = Float32.(volt)
	return (table, hdr) 
end

function main(args)

    settings = ArgParseSettings(prog = "GPPupilDemodulation",
						 #version = @project_version,
						 version = "0.1",
						 add_version = true)

	settings.description =  "Simple tool to demodulate  Gravity metrology table.\n\n"
							
    @add_arg_table! settings begin
        "--suffix", "-s"
			nargs = 1
			action = :store_arg
			arg_type = String
			default = [""]
            help = "Store the demodulated metrogoly in the INPUT.SUFFIX.fits file"
		"--onlyhigh", "-o"
			help = "Demodulate metrology using parameters estimated only on HIGH and NORMAL"
			action = :store_true
		"--nofaint", "-f"
			help = "Do no use the faint mode state to demodulate"
			action = :store_true
		"--recursive", "-r"
			help = "Recursively explore entire directories."
			action = :store_true
		"--verbose", "-v"
			help = "Verbose mode"
			action = :store_true
		"--keepraw", "-k"
			help = "keep raw"
			action = :store_true
		"--center", "-c"
			help = "center voltages: the centering methods are :\n
			- stefan : use stefan measured centers (default)\n \n
			- empirical : estimate the center by fitting a circle on the data\n \n
			- uncentered : no centering\n \n
			- fit : fit all parameters at once as it was done before"
			nargs = 1
			arg_type = String
			action = :store_arg
			default =  ["stefan"]
		# "--overwrite", "-w"
		# 	help = "overwrite the original file"
		# 	action = :store_true
        "--window", "-w"
			nargs = 1
			action = :store_arg
			arg_type = Float64
			default = [0.0]
            help = "Compute demodulation on non overlapping window of WINDOW second"
		"--dir", "-d"
			nargs = 1
			action = :store_arg
			arg_type = String
			default = [pwd()]
			help = "output folder"
		"INPUT"
			nargs = '*'
			arg_type = String
            help = "List of all TARGET to process. In conjunction with -r TARGET can contain directories."
			default = ["."]
    end

    parsed_args = parse_args(args, settings)

	args =parsed_args["INPUT"];

	files = Vector{String}()
	for arg in args
		if isdir(arg) && parsed_args["recursive"]
			files =  vcat(files,[root*"/"*filename for (root, dirs, TARGET) in walkdir(arg) for filename in TARGET  ])
		else
			files =  vcat(files,arg)
		end
	end

	folder = parsed_args["dir"][1]
	if folder[1] !='/'
		folder = pwd() * "/" * folder
		mkpath(folder)
	end

	if !parsed_args["verbose"]
		disable_logging(Logging.Info)
	else
		
		Pkg.status("FITSIO")
	end

	if parsed_args["center"][1] == "stefan"
		offsets = read_stefan_file()
	elseif parsed_args["center"][1] == "uncentered"
		offsets = zeros(ComplexF64,40)
	elseif parsed_args["center"][1] == "empirical"
		offsets = true
	elseif parsed_args["center"][1] == "fit"
		offsets = false
	end

	for filename in files
		if isfile(filename)
			if endswith(filename,SUFFIXES)
				pupmod=false
				metmod="ON"
				f= FITS(filename)
				try (pupmod,) = read_key(f[1],"ESO INS PMC1 MODULATE")
				catch
					@info	"no ESO INS PMC1 MODULATE keyword in $filename"
					continue
				end		
				if pupmod
					faintparam =  nothing
					@info "Processing  $filename"
					tstart = time()

					try (metmod,) = read_key(f[1],"ESO INS MET MODE")
						@info	"$filename use $metmod metrology mode"
						if metmod=="OFF"
							continue
						end
					catch
						@info "No ESO INS MET MODE keyword, mode set to $metmod"
					end

					if  metmod == "FAINT"
						if !parsed_args["nofaint"]
							faintparam = buildfaintparameters(read_header(f[1]));
						else 
							@info "FAINT mode deactivated"
						end
					end
					
					mjd = read_key(f[1],"MJD-OBS")[1]
					metrologyhdu = f["METROLOGY"]
					window = parsed_args["window"][1]
					if window ==0.
						window = nothing
					end
					(table, hdr) = processmetrology(metrologyhdu,mjd; 
										faintparam = faintparam, 
										verb=parsed_args["verbose"], 
										window=window,
										keepraw=parsed_args["keepraw"],
										onlyhigh=parsed_args["onlyhigh"],
										offsets=offsets)
					
					tend = time()
					@info "$filename processed in $(tend-tstart) s"
					
					fname = split(basename(filename),".fits")[1]
					outname = folder *"/" * fname * parsed_args["suffix"][1] *".fits"
					close(f)
					f= FITS(filename)
					g = FITS(outname, "w")
					
					FITScopy!(g,f,"METROLOGY"=>table, "METROLOGY"=>hdr)
					close(g)

					@info " $outname written"
				
				else 
					@info	"ESO INS PMC1 MODULATE set to false in  $filename"
					
				end
				close(f)
			end
		end
	end

end


end

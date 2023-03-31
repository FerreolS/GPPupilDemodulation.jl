module GPPupilDemodulation

using Logging, Pkg
include("FitsUtils.jl")

import .FitsUtils: FITScopy!

include("Modulation.jl")



using ArgParse, FITSIO

const SUFFIXES = [".fits", ".fits.gz","fits.Z"]
const MJD_1970_1_1 = 40587.0
const DAY_TO_SEC = 24*60*60 

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

	mjdobs = hdr["MJD-OBS"]
	rate1 = hdr["ESO INS ANLO3 RATE1"]
	rate2 = hdr["ESO INS ANLO3 RATE2"]

	repeat1 = hdr["ESO INS ANLO3 REPEAT1"]
	repeat2 = hdr["ESO INS ANLO3 REPEAT2"]

	start1 = hdr["ESO INS ANLO3 TIMER1"] - (mjdobs - MJD_1970_1_1)*DAY_TO_SEC
	start2 = hdr["ESO INS ANLO3 TIMER2"] - (mjdobs - MJD_1970_1_1)*DAY_TO_SEC

	voltage1 = hdr["ESO INS ANLO3 VOLTAGE1"]
	voltage2 = hdr["ESO INS ANLO3 VOLTAGE2"]
	timer1 = start1 .+ rate1 .* (0:(repeat1-1))
	timer2 = start2 .+ rate2 .* (0:(repeat2-1))
	return FaintStates(timer1,timer2,voltage1,voltage2)
end

function processmetrology(metrologyhdu::TableHDU;faintparam::Union{Nothing,FaintStates} = nothing, keepraw = false,verb=false,onlyhigh=onlyhigh)
	hdr = read_header(metrologyhdu)
	table = Dict(metrologyhdu)
	time = Float64.(table["TIME"]).*1e-6
	volt = Float64.(table["VOLT"])
	cmplxV = volt[1:2:end,:]' .+ im*volt[2:2:end,:]'
	(output, param,likelihood) = demodulateall(time, cmplxV; faintparam = faintparam,onlyhigh=onlyhigh)

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
	table["VOLT"] = Float32.(volt)

	setindex!(hdr,"GPPupilDemodulation.jl","PROCSOFT")
	for (k,j,i) ∈ Iterators.product((D1,D2,D3,D4),1:4,(FT,SC)) 
		b = param[idx(i,j,k)].b
		ϕ = param[idx(i,j,k)].ϕ
		if (b<0)
			b = -b
			ϕ = rem2pi(ϕ+π,RoundNearest) 
		end
		setindex!(hdr,real(param[idx(i,j,k)].c),"DEMODULATION CENTER X0 $i T$j $k")
		setindex!(hdr,imag(param[idx(i,j,k)].c),"DEMODULATION CENTER Y0 $i T$j $k")
		setindex!(hdr,abs(param[idx(i,j,k)].a),"DEMODULATION AMPLITUDE ABS $i T$j $k")
		setindex!(hdr,angle(param[idx(i,j,k)].a),"DEMODULATION AMPLITUDE ARG $i T$j $k")
		setindex!(hdr,b,"DEMODULATION SIN AMPLITUDE $i T$j $k")
		setindex!(hdr,ϕ,"DEMODULATION SIN PHASE $i T$j $k")
	end
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
			help = "Demodulate only on HIGH metrology  in faint mode"
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
		# "--overwrite", "-w"
		# 	help = "overwrite the original file"
		# 	action = :store_true
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
	for filename in files
		if isfile(filename)
			if endswith(filename,SUFFIXES)
				pupmod=false
				metmod="ON"
				f= FITS(filename)
				try (pupmod,) = read_key(f[1],"ESO INS PMC1 MODULATE")
				catch
					if parsed_args["verbose"]
						println("no ESO INS PMC1 MODULATE keyword in $filename")
					end
					continue
				end		
				if pupmod
					faintparam =  nothing
					if parsed_args["verbose"]
						println("Processing  $filename")
					end
					tstart = time()

					try (metmod,) = read_key(f[1],"ESO INS MET MODE")
						if parsed_args["verbose"]
							println("$filename use $metmod metrology mode")
						end
						if metmod=="OFF"
							continue
						end
					catch
						if parsed_args["verbose"]
							println("No ESO INS MET MODE keyword, mode set to $metmod")
						end
					end

					if metmod == "FAINT"
						faintparam = buildfaintparameters(read_header(f[1]));
					end

					metrologyhdu = f["METROLOGY"]
					(table, hdr) = processmetrology(metrologyhdu; faintparam = faintparam, verb=parsed_args["verbose"], keepraw=parsed_args["keepraw"],onlyhigh=parsed_args["onlyhigh"])
					
					tend = time()
					if parsed_args["verbose"]
						println("$filename processed in $(tend-tstart) s")
					end
					fname = split(basename(filename),".fits")[1]
					outname = folder *"/" * fname * parsed_args["suffix"][1] *".fits"
					close(f)
					f= FITS(filename)
					g = FITS(outname, "w")
					
					FITScopy!(g,f,"METROLOGY"=>table, "METROLOGY"=>hdr)
					close(g)

					if parsed_args["verbose"]
					println(" $outname written")
					end
				else 
					if parsed_args["verbose"]
						println("ESO INS PMC1 MODULATE set to false in  $filename")
					end
				end
				close(f)
			end
		end
	end

end


end

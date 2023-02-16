module GPPupilDemodulation


include("FitsUtils.jl")

import .FitsUtils: FITScopy!

include("Modulation.jl")


using ArgParse, FITSIO

const suffixes = [".fits", ".fits.gz","fits.Z"]

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


function processmetrology(metrologyhdu::TableHDU; verb=Bool)
	hdr = read_header(metrologyhdu)
	table = Dict(metrologyhdu)
	time = Float64.(table["TIME"]).*1e-6
	volt = Float64.(table["VOLT"])
	cmplxV = volt[1:2:end,:]' .+ im*volt[2:2:end,:]'

	(output, param) = demodulateall(time, cmplxV)
	volt[1:2:end,:] .=  real(output)'
	volt[2:2:end,:] .=  imag(output)'
	table["VOLT"] .= Float32.(volt)

	setindex!(hdr,"GPPupilDemodulation.jl","PROCSOFT")
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
		"--recursive", "-r"
			help = "Recursively explore entire directories."
			action = :store_true
		"--verbose", "-v"
			help = "Verbose mode"
			action = :store_true
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

			
	for filename in files
		if isfile(filename)
			if endswith(filename,suffixes)
				pupmod=false
				f= FITS(filename)
				try (pupmod,) = read_key(f[1],"ESO INS PMC1 MODULATE")
				catch
					if parsed_args["verbose"]
						println("no ESO INS PMC1 MODULATE keyword in $filename")
					end
					continue
				end		
				if pupmod
					if parsed_args["verbose"]
						println("Processing  $filename")
					end
					tstart = time()

					metrologyhdu = f["METROLOGY"]
					(table, hdr) = processmetrology(metrologyhdu; verb=parsed_args["verbose"])
					
					tend = time()
					if parsed_args["verbose"]
						println("$filename processed in $(tend-tstart) s")
					end
					fname = split(basename(filename),".fits")[1]
					outname = parsed_args["dir"][1] *"/" * fname * parsed_args["suffix"][1] *".fits"
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

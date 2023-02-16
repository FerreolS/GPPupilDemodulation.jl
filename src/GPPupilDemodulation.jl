module GPPupilDemodulation


include("FitsUtils.jl")

import .FitsUtils: FITScopy!

include("Modulation.jl")


using ArgParse, FITSIO

const suffixes = [".fits", ".fits.gz","fits.Z"]

function processmetrology(metrologyhdu::TableHDU; verb=Bool)
	hdr = read_header(metrologyhdu)
	table = Dict(metrologyhdu)
	time = Float64.(table["TIME"]).*1e-6
	volt = Float64.(table["VOLT"])
	cmplxV = volt[1:2:end,:]' .+ im*volt[2:2:end,:]'
	(output, param) = demodulateall(time, cmplxV)
	volt[1:2:end,:]' .=  real(output)
	volt[2:2:end,:]' .=  imag(output)
	table["VOLT"] .= Float32.(volt)
	setindex!(hdr,"GPPupilDemodulation.jl","PROCSOFT")
	return (table, hdr) 
end

function main(args)

    settings = ArgParseSettings(prog = "GPPupilDemodulation",
						 #version = @project_version,
						 version = "0.1",
						 add_version = true)

	settings.description =  "Simple tool to demodulate  Gravity metrology table.\n\n"*
							"Without any argument, it will display the name and the type of all HDU contained in the files TARGET."
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
				f= FITS(filename)
				try (pupmod,) = read_key(f[1],"ESO INS PMC1 MODULATE")
					if pupmod
						metrologyhdu = f["METROLOGY"]
						println("Processing  $filename")
						(table, hdr) = processmetrology(metrologyhdu; verb=parsed_args["verbose"])
						dir = dirname(filename)
						fname = split(basename(filename),".fits")[1]
						outname = tempdir() * fname *  parsed_args["suffix"]*".fits"
						g = FITS(outname, "w")
						FITScopy!(g,f,"METROLOGY"=>table, "METROLOGY"=>hdr)
					else 
						if parsed_args["verbose"]
							println("ESO INS PMC1 MODULATE set to false in  $filename")
						end
					end
				catch
					if parsed_args["verbose"]
						println("no ESO INS PMC1 MODULATE keyword in $filename")
					end
				end		
			end
		end
	end

end


end

module GPPupilDemodulation

using ArgParse
@enum side FT=0 SC=16
@enum diode D1=1 D2=2 D3=3 D4=4 FC
const ȷ=im

include("FitsUtils.jl")

import .FitsUtils: FITScopy!



include("Modulation.jl")



function main(args)

    settings = ArgParseSettings(prog = "FITSexplore",
						 #version = @project_version,
						 version = "0.2",
						 add_version = true)

	settings.description =  "Simple tool to explore the content of FITS files.\n\n"*
							"Without any argument, it will display the name and the type of all HDU contained in the files TARGET."
    @add_arg_table! settings begin
		"--header", "-d"
			help = "header"
			action = :store_true
			help = "Print the whole FITS header."
        "--stats", "-s"
			action = :store_true
            help = "Print the Statistics of the first HDU"
        "--keyword", "-k"
			nargs = 1
			action = :append_arg
			arg_type = String
            help = "Print the value of the FITS header KEYWORD. This argument can be set multiple times to display several FITS keyword"
        "--filter", "-f"
            help = "filter"
			arg_type = String
			nargs = 2
			metavar = ["KEYWORD", "VALUE"]
			help = "Print all files where the FITS header KEYWORD = VALUE."
		"--recursive", "-r"
			help = "Recursively explore entire directories."
			action = :store_true
		"TARGET"
			nargs = '*'
			arg_type = String
            help = "List of all TARGET to explore. In conjunction with -r TARGET can contain directories."
			default = ["."]
    end

    parsed_args = parse_args(args, settings)

	args =parsed_args["TARGET"];

	files = Vector{String}()
	for arg in args
		if isdir(arg) && parsed_args["recursive"]
			files =  vcat(files,[root*"/"*filename for (root, dirs, TARGET) in walkdir(arg) for filename in TARGET  ])
		else
			files =  vcat(files,arg)
		end
	end


	head::Bool =  parsed_args["header"];
	stats::Bool =  parsed_args["stats"];


end


end

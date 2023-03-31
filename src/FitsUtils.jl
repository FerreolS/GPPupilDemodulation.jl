module FitsUtils
import Base: Dict, names, write

using EasyFITS


name(hdu::FitsHDU) = hdu["EXTNAME"].value(String)

extver(hdu::FitsHDU) = hdu["EXTVER"].value(Integer)

function getunits(hdr::FitsHeader)
	col_units = Dict{String,String}()
	if !haskey(hdr, "TFIELDS") return nothing
	end
	ncol = hdr["TFIELDS"]
	for i ∈ 1:ncol
		if haskey(hdr, "TUNIT$i") 
			push!(col_units,hdr["TTYPE$i"] => hdr["TUNIT$i"])
		end
	end
	return col_units
end

# function rewind(file::FITSFile)
# 	fits_movabs_hdu(f::FITSFile, 1)
# end

function Base.Dict(hdu::FitsTableHDU) 
	D = Dict{String,Any}()
	for name ∈ names(hdu)
		push!(D,name => read(hdu,name))
	end
	return D
end

# Borrowed from https://github.com/JuliaAstro/CFITSIO.jl/issues/19
function fits_create_empty_hdu(f::FitsFile)

    status = Ref{Cint}(0)
    naxesr = C_NULL
    N = 0
    bitpix = bitpix_from_type(Int16)

    ccall(
         (:ffcrimll, libcfitsio),
         Cint,
         (Ptr{Cvoid}, Cint, Cint, Ptr{Int64}, Ref{Cint}),
         f.ptr,
         bitpix,
         N,
         naxesr,
         status,
     )

    fits_assert_ok(status[])
end

# Base.write(f::FitsFile, Nothing; kwds...) = write(f;kwds...)

# function Base.write(f::FitsFile;
# 	header::Union{Nothing, FitsHeader}=nothing,
# 	name::Union{Nothing, String}=nothing,
# 	ver::Union{Nothing, Integer}=nothing)
#     status = fits_create_empty_hdu(f.fitsfile)

#     if isa(header, FitsHeader)
#         FITSIO.fits_write_header(f.fitsfile, header, true)
#     end
#     if isa(name, String)
#         fits_update_key(f.fitsfile, "EXTNAME", name)
#     end
#     if isa(ver, Integer)
#         fits_update_key(f.fitsfile, "EXTVER", ver)
#     end
#     nothing
# end

FITScopy!(dst::FitsFile,src::FitsFile) = FITScopy!(dst,src,(),())
function FITScopy!(dst::FitsFile,src::FitsFile,content::Union{Pair{String,Union{T1,T2}},NTuple{N,Pair{String,Union{T1,T2}}}})  where {N,T1<:AbstractDict,T2<:AbstractArray} 
    FITScopy!(dst,src,content,())
end

function FITScopy!(dst::FitsFile,
		src::FitsFile,
		content::Union{Pair{String,Union{T1,T2}},NTuple{N,Pair{String,Union{T1,T2}}}},
		header::Union{Pair{String,FitsHeader},NTuple{M,Pair{String,FitsHeader}}}) where {M,N,T1<:AbstractDict,T2<:AbstractArray} 
 
		FITScopy!(dst,src,content,header,nothing)
	end
	
	
	function FITScopy!(dst::FitsFile,
			src::FitsFile,
			content::Union{Pair{String,Union{T1,T2}},NTuple{N,Pair{String,Union{T1,T2}}}},
			header::Union{Pair{String,FitsHeader},NTuple{M,Pair{String,FitsHeader}}},
			units::Union{Nothing,Pair{String,Dict{String,String}},NTuple{O,Pair{String,Dict{String,String}}}} ) where {M,N,O,T1<:AbstractDict,T2<:AbstractArray} 	
	 
	Dcontent = Dict(content)
	Dheader = Dict(header)
	if isnothing(units)
		Dunits = Dict()
	else
		Dunits = Dict(units)
	end

	hdr=nothing
	cunits = nothing

	for hdu ∈ src
		hduname =name(hdu)

		hdr = pop!(Dheader,hduname,read_header(hdu))
		
		if haskey(Dcontent,hduname)
			cntnt = pop!(Dcontent,hduname)
		else
			if isa(hdu,TableHDU)
				cntnt = Dict(hdu)
			elseif isa(hdu, ImageHDU) 
				if (size(hdu) == ())
					cntnt = nothing
				else
					cntnt = read(hdu)
				end
			end
		end
		if isa(hdu,TableHDU)
			if haskey(Dunits,hduname)
				cunits = pop!(Dcontent,hduname)
			else
				cunits = getunits(hdr)
			end
			write(dst,cntnt;header=hdr,name=hduname,ver=extver(hdu),units=cunits)
		else
			write(dst,cntnt;header=hdr,name=hduname,ver=extver(hdu))
		end
	end
	for (key,cntnt) ∈ Dcontent
		hdr = pop!(Dheader,key,nothing)
		write(dst,cntnt;header=hdr,name=key)
		delete!(Dcontent,key)
	end
	
	
	for (key,hdr) ∈ Dheader
		hdr = pop!(Dheader,key,nothing)
		write(dst;header=hdr,name=key)
		delete!(Dheader,key)
	end
	
	
	
end

end
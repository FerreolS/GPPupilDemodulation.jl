module FitsUtils
import Base: Dict, names, write

import CFITSIO: FITSFile, bitpix_from_type,libcfitsio,fits_assert_ok,fits_update_key
using FITSIO
#import FITSIO: fits_create_img

name(hdu::HDU) = FITSIO.fits_try_read_extname(hdu.fitsfile)

names(hdu::TableHDU) = FITSIO.colnames(hdu)

extver(hdu::HDU) = FITSIO.fits_try_read_extver(hdu.fitsfile)

function Dict(hdu::TableHDU) 
	D = Dict{String,Any}()
	for name ∈ names(hdu)
		push!(D,name => read(hdu,name))
	end
	return D
end

# Borrowed from https://github.com/JuliaAstro/CFITSIO.jl/issues/19
function fits_create_empty_hdu(f::FITSFile)

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

write(f::FITS, Nothing; kwds...) = write(f;kwds...)

function write(f::FITS;
	header::Union{Nothing, FITSHeader}=nothing,
	name::Union{Nothing, String}=nothing,
	ver::Union{Nothing, Integer}=nothing)
    status = fits_create_empty_hdu(f.fitsfile)

    if isa(header, FITSHeader)
        FITSIO.fits_write_header(f.fitsfile, header, true)
    end
    if isa(name, String)
        fits_update_key(f.fitsfile, "EXTNAME", name)
    end
    if isa(ver, Integer)
        fits_update_key(f.fitsfile, "EXTVER", ver)
    end
    nothing
end

FITScopy!(dst::FITS,src::FITS) = FITScopy!(dst,src,(),())
function FITScopy!(dst::FITS,src::FITS,content::Union{Pair{String,Union{T1,T2}},NTuple{N,Pair{String,Union{T1,T2}}}})  where {N,T1<:AbstractDict,T2<:AbstractArray} 
    FITScopy!(dst,src,content,())
end
function FITScopy!(dst::FITS,
		src::FITS,
		content::Union{Pair{String,Union{T1,T2}},NTuple{N,Pair{String,Union{T1,T2}}}},
		header::Union{Pair{String,FITSHeader},NTuple{M,Pair{String,FITSHeader}}} ) where {M,N,T1<:AbstractDict,T2<:AbstractArray} 
	 
	Dcontent = Dict(content)
	Dheader = Dict(header)

	for hdu ∈ src
		hdr = pop!(Dheader,name(hdu),read_header(hdu))
		if haskey(Dcontent,name(hdu))
			cntnt = pop!(Dcontent,name(hdu))
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
		
		write(dst,cntnt;header=hdr,name=name(hdu),ver=extver(hdu))
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
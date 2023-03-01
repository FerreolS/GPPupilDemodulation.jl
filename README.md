# GPPupilDemodulation.jl
Gravity+ Metrology metrology demodulation module


### To use it:
- add it to your julia using `add https://github.com/FerreolS/GPPupilDemodulation.jl.git`
- download the GPPupilDemodulation file in the bin folder
- make it executable using `chmod +x GPPupilDemodulation`  and add it to your path
- find the help by typing `GPPupilDemodulation -h`

### Example:
`GPPupilDemodulation -r -v  -d "output_folder" /data/ESO_archive_night/2023-01-05`
will demodulate all the file contained in folder `/data/ESO_archive_night/2023-01-05`  and place the demodulated files in `output_folder`

### Algo Benchmark

```julia


julia> @btime (out,param) = GPPupilDemodulation.demodulateall(algo=:newoa,times,cmplxV);
  1.706 s (4494 allocations: 2.94 GiB)

julia> @btime (out,param) = GPPupilDemodulation.demodulateall(algo=:simplex,times,cmplxV);
  7.107 s (17498 allocations: 8.74 GiB)

julia> @btime (out,param) = GPPupilDemodulation.demodulateall(algo=:vmlmbZ,times,cmplxV);
  7.101 s (170702 allocations: 24.10 GiB)

```

### Notebook
[Notebook](https://github.com/FerreolS/GPPupilDemodulation.jl/blob/notebooks/MetrologyModulation.ipynb) on the metrology data.
As Nbviewer cannot show plots properly, most of them can be seen 
[here](https://jovian.ml/ferreols/gravity-metrology-modulation/v/1)
### Notes demodulation
[![Notes about GRAVITY+ metrology demodulation](https://github.com/FerreolS/GPPupilDemodulation.jl/blob/gh-pages/GPPupilDemodulation.svg)](https://github.com/FerreolS/GPPupilDemodulation.jl/blob/gh-pages/GPPupilDemodulation.pdf)

# GPPupilDemodulation.jl
Gravity+ Metrology Pupil demodulation module


### To use it:
- add it to your julia using `add https://github.com/FerreolS/GPPupilDemodulation.jl.git`
- download the GPPupilDemodulation file in the bin folder
- make it executable using `chmod +x GPPupilDemodulation`  and add it to your path
- find the help by typing `GPPupilDemodulation -h`

### Example:
`GPPupilDemodulation -r -v  -d "output_folder" /data/ESO_archive_night/2023-01-05`
will demodulate all the file contained in folder `/data/ESO_archive_night/2023-01-05`  and place the demodulated files in `output_folder`

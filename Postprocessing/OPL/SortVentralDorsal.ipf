#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// script sorts scan field into dorsal and ventral based on color preference

function SortVentralDorsal()

// Set parameters

variable Ventral_Cutoff = -0.1 // spectral contrast
variable Dorsal_Cutoff = 0.1 // spectral contrast

// Define input

wave CenterSurround_SC, RoisPerField
variable nRois = DimSize(CenterSurround_SC,0)
variable nFields = DimSize(RoisPerField,0)
variable rr, ff

// Sort scan fields

make/o/n=(nFields) Fields_SC = nan

make/o/n=(nRois) ROIs_DorsalVentral = nan

for (ff=0; ff<nFields; ff+=1)
	// Find ROIs of this field
	if (ff==0)
		variable start = 0
		variable stop = start + RoisPerField[ff]
	else
		start = sum(RoisPerField,0,ff-1)
		stop = start + RoisPerField[ff]
	endif
	// Get mean spectral contrast of center responses
	make/o/n=(RoisPerField[ff]) CurrentWave = CenterSurround_SC[p+start][0]
	Wavestats/Q CurrentWave
	Fields_SC[ff] = V_avg
	// Assign ROIs to dorsal (=1) and ventral (=0) retina
	for (rr=start; rr<stop; rr+=1)
		if (V_avg > Dorsal_Cutoff)
			ROIs_DorsalVentral[rr] = 1
		elseif (V_avg < Ventral_Cutoff)
			ROIs_DorsalVentral[rr] = 0
		endif
	endfor
endfor

killwaves CurrentWave

end
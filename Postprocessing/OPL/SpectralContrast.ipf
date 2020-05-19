#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// script calculates spectral contrast based on full-field and center-surround flashes

function SpectralContrast()

// Set parameters

variable FullField_Threshold = 0.25 // quality threshold blue green white full-field flashes
variable CenterSurround_Threshold = 0.25 // quality threshold blue green center-surround flashes

// Define input

wave FullField_Traces, FullField_Quality, CenterSurround_Quality, CenterSurround_Traces
variable nF_BGW = DimSize(FullField_Traces,0)
variable nF_BG_CS = DimSize(CenterSurround_Traces,0)
variable nRois = DimSize(FullField_Traces,1)
variable rr

// Estimate spectral contrast based on full-field flashes

variable ReadOutLength = 600
variable Baseline_Start = 870
variable BaselineLength = 60
variable BlueStart = 0
variable GreenStart = 1000

make/o/n=(nRois) SC_FullField = nan
make/o/n=(nRois,2) FullField_Areas = nan

for (rr=0; rr<nRois; rr+=1)
	if (FullField_Quality[rr] > FullField_Threshold)
		// Subtract baseline
		make/o/n=(nF_BGW) CurrentTrace = FullField_Traces[p][rr]
		variable Baseline = mean(CurrentTrace,Baseline_Start,Baseline_Start+BaselineLength)
		CurrentTrace-=Baseline
		// Get areas
		variable AreaBlue = area(CurrentTrace,BlueStart,BlueStart+ReadOutLength)
		variable AreaGreen = area(CurrentTrace,GreenStart,GreenStart+ReadOutLength)		
		FullField_Areas[rr][0] = AreaBlue
		FullField_Areas[rr][1] = AreaGreen
		// Estimate spectral contrast		
		if (AreaBlue < 0 && AreaGreen < 0) // equation 6-a
			SC_FullField[rr] = (AreaGreen-AreaBlue)/(AreaGreen+AreaBlue)
		elseif (AreaBlue < 0 && AreaGreen > 0) // equation 6-c
			SC_FullField[rr] = -1 - ((AreaGreen)/abs(AreaBlue))
		elseif (AreaBlue > 0 && AreaGreen < 0) // equation 6-b
			SC_FullField[rr] = 1 + ((AreaBlue)/abs(AreaGreen))
		endif		
	endif
endfor

killwaves CurrentTrace

// Estimate spectral contrast based on center-surround flashes

BlueStart = 2000
GreenStart = 0
variable BlueStart_Surround = 3000
variable GreenStart_Surround = 1000
BaselineLength = 200
variable ThresholdSurroundArea = 0.1

make/o/n=(nRois,2) CenterSurround_SC = nan
make/o/n=(nRois,4) CenterSurround_Area = nan

for (rr=0; rr<nRois; rr+=1)
	if (CenterSurround_Quality[rr] > CenterSurround_Threshold)
		make/o/n=(nF_BG_CS) CurrentTrace = CenterSurround_Traces[p][rr]		
		// Get areas center
		Baseline = CurrentTrace[0]	
		CurrentTrace-=Baseline	// subtract baseline before green center response
		AreaGreen = area(CurrentTrace,GreenStart,GreenStart+ReadOutLength)
		Baseline = mean(CurrentTrace,BlueStart-BaselineLength, BlueStart)
		CurrentTrace-=Baseline	// subtract baseline before blue center response
		AreaBlue = area(CurrentTrace,BlueStart,BlueStart+ReadOutLength)		
		CenterSurround_Area[rr][0] = AreaGreen
		CenterSurround_Area[rr][1] = AreaBlue
		// Estimate spectral contrast center
		if (AreaBlue < 0 && AreaGreen < 0)
			CenterSurround_SC[rr][0] = (AreaGreen-AreaBlue)/(AreaGreen+AreaBlue) // equation 6-a
		elseif (AreaBlue < 0 && AreaGreen > 0)
			CenterSurround_SC[rr][0] = -1 - ((AreaGreen)/abs(AreaBlue)) // equation 6-c
		elseif (AreaBlue > 0 && AreaGreen < 0)
			CenterSurround_SC[rr][0] = 1 + ((AreaBlue)/abs(AreaGreen)) // equation 6-b
		endif
		// Get areas surround
		Baseline = mean(CurrentTrace,BlueStart_Surround-BaselineLength,BlueStart_Surround)
		CurrentTrace-=Baseline // subtract baseline before blue surround response
		variable AreaBlue_Surround = area(CurrentTrace,BlueStart_Surround,BlueStart_Surround+ReadOutLength)		
		Baseline = mean(CurrentTrace,GreenStart_Surround-BaselineLength,GreenStart_Surround)
		CurrentTrace-=Baseline // subtract baseline before green surround response
		variable AreaGreen_Surround = area(CurrentTrace,GreenStart_Surround,GreenStart_Surround+ReadOutLength)
		// Only consider increases in glutamate as surround response (=depolarization)		
		if (AreaGreen_Surround > 0) 
			CenterSurround_SC[rr][2] = AreaGreen_Surround 
		endif
		if (AreaBlue_Surround > 0)
			CenterSurround_SC[rr][3] = AreaBlue_Surround
		endif	
		// Get dominant center area
		variable CurrentCenterArea = 0
		if (CenterSurround_SC[rr][0] < CenterSurround_SC[rr][1])
			CurrentCenterArea = abs(CenterSurround_SC[rr][0])
		elseif (CenterSurround_SC[rr][0] > CenterSurround_SC[rr][1])
			CurrentCenterArea = abs(CenterSurround_SC[rr][1])
		endif
		// Estimate spectral contrast
		// Only consider cells with surround area > 0.1 of center area (ThresholdSurroundArea)
		if (AreaBlue_Surround/CurrentCenterArea > ThresholdSurroundArea || AreaGreen_Surround/CurrentCenterArea > ThresholdSurroundArea || (AreaBlue_Surround+AreaGreen_Surround)/CurrentCenterArea > ThresholdSurroundArea)
			if (AreaBlue_Surround > 0 && AreaGreen_Surround > 0)
				CenterSurround_SC[rr][1] = (AreaGreen_Surround-AreaBlue_Surround)/(AreaGreen_Surround+AreaBlue_Surround) // equation 6-a
			elseif (AreaBlue_Surround > 0 && AreaGreen_Surround < 0)
				CenterSurround_SC[rr][1] = -1 - (AreaBlue_Surround)/abs(AreaGreen_Surround) // equation 6-c
			elseif (AreaBlue_Surround < 0 && AreaGreen_Surround > 0) 
				CenterSurround_SC[rr][1] = 1 + (AreaGreen_Surround)/abs(AreaBlue_Surround) // equation 6-b
			endif
		endif		
	endif
endfor

killwaves CurrentTrace

end
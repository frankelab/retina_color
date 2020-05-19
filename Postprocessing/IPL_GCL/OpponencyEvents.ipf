#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// script calculates opponency of full-field events

function OpponencyEvents()

// Set parameters

variable FrStart = 400 // frame of kernel start correlation
variable FrEnd = 800  // frame of kernel end correlation
variable BaselineLength = 50 // number of frames to determine baseline
variable ThresholdKernels = 0.6 // quality threshold kernels
variable ThresholdEvent = 0.25 // quality threshold events

// Define input

wave CS_Events_On, AreaKernels, RoiInfo_corrected, Events_Qi_On, Kernels_Qi, CS_Events_Off, Events_Qi_Off
variable nF = DimSize(CS_Events_On,0)
variable nConditions = DimSize(CS_Events_On,1)
variable nRois = DimSize(CS_Events_On,2)
variable rr

// Estimate opponency of onset and offset full-field events

make/o/n=(nRois,2) Events_On_Areas = nan
make/o/n=(nRois,2) Events_Off_Areas = nan

make/o/n=(nRois) Events_On_Correlation = nan
make/o/n=(nRois) Events_Off_Correlation = nan

for (rr=0; rr<nRois; rr+=1)
	if (Kernels_Qi[rr][0] > ThresholdKernels || Kernels_Qi[rr][2] > ThresholdKernels) // if green or UV center above threshold
		// Onset events
		// Subtract baseline
		make/o/n=(FrEnd-FrStart) CurrentBlue = CS_Events_On[p+FrStart][2][rr]
		variable Baseline = mean(CurrentBlue,0,BaselineLength)
		CurrentBlue-=Baseline
		make/o/n=(FrEnd-FrStart) CurrentGreen = CS_Events_On[p+FrStart][5][rr]
		Baseline = mean(CurrentGreen,0,BaselineLength)
		CurrentGreen-=Baseline
		// Correlate
		StatsLinearCorrelationTest/Q CurrentBlue, CurrentGreen
		wave W_StatsLinearCorrelationTest
		// Get area
		CurrentBlue = abs(CurrentBlue)
		variable AreaBlue = area(CurrentBlue)
		CurrentGreen = abs(CurrentGreen)
		variable AreaGreen = area(CurrentGreen)
		Events_On_Areas[rr][0] = AreaBlue
		Events_On_Areas[rr][1] = AreaGreen
		// Save correlation
		if (Events_Qi_On[rr][2] > ThresholdEvent && Events_Qi_On[rr][5] > ThresholdEvent)
			Events_On_Correlation[rr] = W_StatsLinearCorrelationTest[1]
		endif
		
		// Offset events
		// Subtract baseline
		make/o/n=(FrEnd-FrStart) CurrentBlue = CS_Events_Off[p+FrStart][2][rr]
		Baseline = mean(CurrentBlue,0,BaselineLength)
		CurrentBlue-=Baseline
		make/o/n=(FrEnd-FrStart) CurrentGreen = CS_Events_Off[p+FrStart][5][rr]
		Baseline = mean(CurrentGreen,0,BaselineLength)
		CurrentGreen-=Baseline
		// Correlate
		StatsLinearCorrelationTest/Q CurrentBlue, CurrentGreen
		wave W_StatsLinearCorrelationTest
		// Get area
		CurrentBlue = abs(CurrentBlue)
		AreaBlue = area(CurrentBlue)
		CurrentGreen = abs(CurrentGreen)
		AreaGreen = area(CurrentGreen)
		Events_Off_Areas[rr][0] = AreaBlue
		Events_Off_Areas[rr][1] = AreaGreen
		// Save correlation
		if (Events_Qi_Off[rr][2] > ThresholdEvent && Events_Qi_Off[rr][5] > ThresholdEvent)
			Events_Off_Correlation[rr] = W_StatsLinearCorrelationTest[1]
		endif
	endif
endfor

killwaves CurrentBlue, CurrentGreen, W_StatsLinearCorrelationTest

// Combine onset and offset correlation (min(corr_onset, corr_offset))

make/o/n=(nRois) Events_Correlation = nan

for (rr=0; rr<nRois; rr+=1)
	if (numtype(Events_Off_Correlation[rr]) != 2 && numtype(Events_On_Correlation[rr]) != 2)
		Events_Correlation[rr] = min(Events_Off_Correlation[rr], Events_On_Correlation[rr])
	elseif (numtype(Events_Off_Correlation[rr]) != 2 && numtype(Events_On_Correlation[rr]) == 2)
		Events_Correlation[rr] = Events_Off_Correlation[rr]
	elseif (numtype(Events_Off_Correlation[rr]) == 2 && numtype(Events_On_Correlation[rr]) != 2)
		Events_Correlation[rr] = Events_On_Correlation[rr]
	endif
endfor

end
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// script calculates parameters related to quality of kernels and events

function KernelQuality()

// Set parameters

variable BaselineLength = 200 // number of points used to determine baseline
variable FrAreaStart = 400 // frame start kernel
variable FrAreaEnd = 900 // frame end kernel

// Define input

wave CS_Kernels, CS_Kernels_Drug, CS_Events_On, CS_Events_Off
variable nF = DimSize(CS_Kernels,0)
variable nKernels = DimSize(CS_Kernels,1)
variable nRois = DimSize(CS_Kernels,2)
variable rr, kk

// Estimate quality of stimulus kernels

duplicate/o CS_Kernels CS_Kernels_norm
make/o/n=(nRois, nKernels) AreaKernels = nan
make/o/n=(nRois, nKernels) Kernels_SD = nan
make/o/n=(nRois, nKernels) Kernels_Qi = nan

for (rr=0; rr<nRois; rr+=1)
	for (kk=0; kk<nKernels; kk+=1)
		make/o/n=(nF) CurrentKernel = CS_Kernels[p][kk][rr]
		// subtract baseline
		variable Baseline = mean(CurrentKernel,0,BaselineLength)
		CurrentKernel-=Baseline
		// save normalized kernels
		CS_Kernels_norm[][kk][rr] = CurrentKernel[p]
		// calculate area
		CurrentKernel = abs(CurrentKernel)		
		AreaKernels[rr][kk] = area(CurrentKernel,FrAreaStart,FrAreaEnd)	
		// calculate quality
		variable AreaCurrentBaseline = area(CurrentKernel,0,FrAreaStart)
		AreaCurrentBaseline+=area(CurrentKernel,FrAreaEnd,nF)
		AreaCurrentBaseline/=(FrAreaStart+nF-FrAreaEnd)		
		Kernels_Qi[rr][kk] = 1 - AreaCurrentBaseline/(area(CurrentKernel,FrAreaStart,FrAreaEnd)/(FrAreaEnd-FrAreaStart))
	endfor
endfor

killwaves CurrentKernel

// Estimate quality of event kernels

BaselineLength = 100
FrAreaStart = 400
FrAreaEnd = 800

nF = DimSize(CS_Events_On,0)
nKernels = DimSize(CS_Events_On,1)

duplicate/o CS_Events_On CS_Events_On_norm
duplicate/o CS_Events_Off CS_Events_Off_norm

make/o/n=(nRois, nKernels) AreaEvents_On = nan
make/o/n=(nRois, nKernels) Events_Qi_On = nan

make/o/n=(nRois, nKernels) AreaEvents_Off = nan
make/o/n=(nRois, nKernels) Events_Qi_Off = nan

for (rr=0; rr<nRois; rr+=1)
	for (kk=0; kk<nKernels; kk+=1)
		make/o/n=(nF) CurrentKernel = CS_Events_On[p][kk][rr]
		// subtract baseline
		Baseline = mean(CurrentKernel,0,BaselineLength)
		CurrentKernel-=Baseline
		// save normalized kernels
		CS_Events_On_norm[][kk][rr] = CurrentKernel[p]
		// calculate area
		CurrentKernel = abs(CurrentKernel)
		AreaEvents_On[rr][kk] = area(CurrentKernel,FrAreaStart,FrAreaEnd)
		AreaCurrentBaseline = area(CurrentKernel,0,FrAreaStart)
		AreaCurrentBaseline+=area(CurrentKernel,FrAreaEnd,nF)
		AreaCurrentBaseline/=(FrAreaStart+nF-FrAreaEnd)
		Events_Qi_On[rr][kk] = 1 - AreaCurrentBaseline/(area(CurrentKernel,FrAreaStart,FrAreaEnd)/(FrAreaEnd-FrAreaStart))
		
		make/o/n=(nF) CurrentKernel = CS_Events_Off[p][kk][rr]
		// subtract baseline
		Baseline = mean(CurrentKernel,0,BaselineLength)
		CurrentKernel-=Baseline
		// save normalized kernels
		CS_Events_Off_norm[][kk][rr] = CurrentKernel[p]
		// calculate area
		CurrentKernel = abs(CurrentKernel)
		AreaEvents_Off[rr][kk] = area(CurrentKernel,FrAreaStart,FrAreaEnd)
		AreaCurrentBaseline = area(CurrentKernel,0,FrAreaStart)
		AreaCurrentBaseline+=area(CurrentKernel,FrAreaEnd,nF)
		AreaCurrentBaseline/=(FrAreaStart+nF-FrAreaEnd)
		Events_Qi_Off[rr][kk] = 1 - AreaCurrentBaseline/(area(CurrentKernel,FrAreaStart,FrAreaEnd)/(FrAreaEnd-FrAreaStart))		
	endfor
	make/o/n=(nF, nKernels) CurrentKernel = CS_Events_On_norm[p][q][rr]
	Wavestats/Q CurrentKernel
	CurrentKernel-=V_min
	CurrentKernel/=V_max - V_min
	CS_Events_On_norm[][][rr] = CurrentKernel[p][q]
	
	make/o/n=(nF, nKernels) CurrentKernel = CS_Events_Off_norm[p][q][rr]
	Wavestats/Q CurrentKernel
	CurrentKernel-=V_min
	CurrentKernel/=V_max - V_min
	CS_Events_Off_norm[][][rr] = CurrentKernel[p][q]
endfor

killwaves CurrentKernel

end
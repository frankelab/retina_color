#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// script calculates spectral contrast based on stimulus kernels

function SpectralContrast()

// Set parameters

variable ThresholdKernels_Center = 0.6 // quality threshold center kernels
variable ThresholdKernels_Surround = 0.6 // quality threshold surround kernels
variable Threshold_Corr = 0.5 // correlation threshold of center and surround to determine which SC equation 

// Define input

wave AreaKernels, RoiInfo, AreaKernels_Drug, Kernels_Qi, Kernels_Correlations
variable nKernels = DimSize(AreaKernels,1)
variable nRois = DimSize(AreaKernels,0)
variable rr

// Estimate spectral contrast using stimulus kernels

make/o/n=(nRois,2) Kernels_BGi = nan

for (rr=0; rr<nRois; rr+=1)
	if ((Kernels_Qi[rr][0] > ThresholdKernels_Center) || (Kernels_Qi[rr][2] > ThresholdKernels_Center)) // if green or UV center kernel > threshold
		Kernels_BGi[rr][0] = (AreaKernels[rr][2] - AreaKernels[rr][0])/(AreaKernels[rr][2] + AreaKernels[rr][0]) // center spectral contrast
		// test that at least one surround kernel anticorrelated to center and > threshold
		if ((Kernels_Correlations[rr][1] < -Threshold_Corr && Kernels_Qi[rr][1] > ThresholdKernels_Surround) || (Kernels_Correlations[rr][3] < -Threshold_Corr && Kernels_Qi[rr][3] > ThresholdKernels_Surround))
			if (Kernels_Correlations[rr][1] < -Threshold_Corr && Kernels_Correlations[rr][3] > Threshold_Corr) // equation 6-c
				Kernels_BGi[rr][1] = -1 - (AreaKernels[rr][3]/AreaKernels[rr][1])
			elseif (Kernels_Correlations[rr][1] > Threshold_Corr && Kernels_Correlations[rr][3] < -Threshold_Corr) // equation 6-b
				Kernels_BGi[rr][1] = 1 + (AreaKernels[rr][1]/AreaKernels[rr][3])
			elseif (Kernels_Correlations[rr][1] < -Threshold_Corr && Kernels_Correlations[rr][3] < -Threshold_Corr) // equation 6-a
				Kernels_BGi[rr][1] = (AreaKernels[rr][3] - AreaKernels[rr][1])/(AreaKernels[rr][3] + AreaKernels[rr][1])		
			endif
		endif
	endif
endfor

end
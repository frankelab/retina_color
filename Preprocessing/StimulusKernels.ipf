#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// script calculates the temporal kernel (center and surround) for UV or green stimulation

function StimulusKernels()

// Set parameters

variable Color = 1 // 0: blue, 1: green

variable TriggerThreshold = 13000 // for trigger detection, arbitrary value
variable StimFreq = 5 // in Hz, 5 Hz for GCL and 10 Hz for IPL recordings
variable KernelLength = 1000 // in 2 and 1.6 ms points for GCL and IPL recordings
variable nF_Jump = 0 // number of frames jumped after detecting an event
variable StimulusOffset = 34 // time between trigger and stimulus onset, setup specific, in ms
variable LineDuration = 2 // in ms
variable nStimuli = 2 // center and surround
variable BaselineLength = 200 // number of frames for baseline of kernels

// Define input

string TracesName = "wDataCh0_Ave_QA_Nor" // traces of all ROIs
duplicate/o $TracesName Traces
string TriggerName = "wDataCh2" // data channel to detect triggers
duplicate/o $TriggerName Triggers
string RoiName = "ROIs" // name of ROI mask
duplicate/o $RoiName ROI

wave CS_Stimulus // load sequence of noise flicker

// Get data dimensions

variable nF = Dimsize(Traces,0)
variable nRois = Dimsize(Traces,1)
variable nRows = Dimsize(Triggers,0)
variable nLines = Dimsize(Triggers,1)
variable rr, ff, ll, tt, ss

// Detrend data

duplicate/o Traces TracesDetrended

for (rr=0;rr<nRois;rr+=1)
     make/o/n=(nF) CurrentTrace = Traces[p][rr]
     Smooth 2^15-1, CurrentTrace
     TracesDetrended[][rr] = Traces[p][rr] - CurrentTrace[p]
endfor

killwaves CurrentTrace, Traces

// Get ROI positions

setscale x, 0, nRows, ROI
setscale y, 0, nLines, ROI
GeometricCenter(ROI) // uses SARFIA built-in function 
wave GeoC

// Find triggers

variable nTriggers = 0
make/o/n=(5000) TriggerTimes = 0
for (ff=0;ff<nF;ff+=1)
	for (ll=0;ll<nLines;ll+=1)
		if (Triggers[0][ll][ff] > TriggerThreshold)
			TriggerTimes[nTriggers] = (ff*nLines) + ll + StimulusOffset/LineDuration // trigger times line precision, 2 or 1.6 ms
			nTriggers+=1
			ff+=1
			ll+=nLines/2
		endif
	endfor
endfor

Redimension/N=(nTriggers) TriggerTimes

// Make stimulus

make/o/n=(nF*nLines, nStimuli) Stimulus = 0

if (Color == 0)
	variable CurrentStart = 0
elseif (Color == 1)
	CurrentStart = nStimuli
endif

for (ss=0; ss<nStimuli; ss+=1)
	for (tt=0; tt<nTriggers-1; tt+=1)
		for (ff=1; ff<StimFreq+1; ff+=1)
			Multithread Stimulus[TriggerTimes[tt]+(ff-1)/StimFreq*(TriggerTimes[tt+1]-TriggerTimes[tt]),TriggerTimes[tt]+ff/StimFreq*(TriggerTimes[tt+1]-TriggerTimes[tt])][ss]=CS_Stimulus[tt*StimFreq+ff-1][ss+CurrentStart]
		endfor
	endfor
endfor

print "Done making stimuli"

// Calculate kernel

make/o/n=(KernelLength, nStimuli, nRois) Kernels = 0

for (rr=0;rr<nRois;rr+=1)
	// Scale trace
	make/o/n=(nF) CurrentTrace = TracesDetrended[p][rr]
	duplicate/o CurrentTrace CurrentTraceDif
	Differentiate CurrentTrace/D=CurrentTraceDif
	duplicate/o CurrentTraceDif CurrentTraceDifSort
	CurrentTraceDifSort = abs(CurrentTraceDifSort)
	Sort CurrentTraceDifSort, CurrentTraceDifSort
	variable SD = CurrentTraceDifSort[nF/2]/0.6745 // see Baden et al. Nature 2016
	// Find events and average stimulus
	make/o/n=(KernelLength, nStimuli) CurrentKernel = 0
	variable nEvents = 0
	for (ff=TriggerTimes[0]/nLines;ff<TriggerTimes[nTriggers-1]/nLines;ff+=1)
		if (CurrentTraceDif[ff] > SD)
			CurrentKernel[][]+=Stimulus[p-(0.75*KernelLength) + (ff*nLines)][q]*(CurrentTraceDif[ff])
			nEvents+=1
			ff+=nF_Jump
		endif
	endfor
	// Normalize by number of events and subtract baseline
	CurrentKernel/=nEvents
	variable Baseline = mean(CurrentKernel,0,BaselineLength)
	CurrentKernel-=Baseline
	Kernels[][][rr] = CurrentKernel[p][q]
	
	print "Roi", rr, ":", nEvents
endfor

killwaves CurrentTrace, CurrentTraceDif, CurrentTraceDifSort, CurrentKernel

end
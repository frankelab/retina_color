#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// script calculates the events (center, surround and full-field) for UV or green stimulation

function EventKernels()

// Set parameters

variable KernelLength = 2000 // in 2 and 1.6 ms points for GCL and IPL recordings
variable nConditions = 3 // center, surround and full-field
variable StimLength = 100 // in ms, 100 and 200 ms for IPL and GCL recordings
variable LineDuration = 2 // in ms 
variable Threshold = (StimLength/LineDuration) + 1 // theshold for detecting On times
variable StimLength_pp = StimLength/LineDuration

// Make step stimulus to convolve

make/o/n=(StimLength_pp+50) StimOn = 0
StimOn[25,25+StimLength_pp] = 1

// Define input

wave Stimulus, TracesDetrended, TriggerTimes
variable nF = DimSize(Stimulus,0)
variable nRois = DimSize(TracesDetrended,1)
variable nTriggers = DimSize(TriggerTimes,0)
variable rr, ff

// Convolve stimulus with step stimulus

make/o/n=(nF) CentreConvolved = Stimulus[p][0]
Convolve StimOn, CentreConvolved
WaveStats/q CentreConvolved
make/o/n=(nF) SurroundConvolved = Stimulus[p][1]
Convolve StimOn, SurroundConvolved

// Find stimulus conditions and estimate events

make/o/n=(KernelLength, nConditions, nRois) Events_On = 0 // Onset events
make/o/n=(KernelLength, nConditions, nRois) Events_Off = 0 // Offset events

for (rr=0; rr<nRois; rr+=1)
	make/o/n=(DimSize(TracesDetrended,0)) CurrentTrace = TracesDetrended[p][rr]
	Interpolate2/n=(nF)/t=1/y=CurrentTrace_L CurrentTrace
	variable nEvents = 0
	// On center
	for (ff=TriggerTimes[0]+KernelLength; ff<TriggerTimes[nTriggers-1]-KernelLength; ff+=1)
		if (CentreConvolved[ff] == Threshold)
			Events_On[][0][rr]+=CurrentTrace_L[p+ff-(0.5*KernelLength)]
			nEvents+=1
			ff+=StimLength_pp
		endif
	endfor
	Events_On[][0][rr]/=nEvents
	// Off center
	for (ff=TriggerTimes[0]+KernelLength; ff<TriggerTimes[nTriggers-1]-KernelLength; ff+=1)
		if (CentreConvolved[ff] < 1)
			Events_Off[][0][rr]+=CurrentTrace_L[p+ff-(0.5*KernelLength)]
			nEvents+=1
			ff+=StimLength_pp
		endif
	endfor
	Events_Off[][0][rr]/=nEvents
	// On surround
	nEvents = 0
	for (ff=TriggerTimes[0]+KernelLength; ff<TriggerTimes[nTriggers-1]-KernelLength; ff+=1)
		if (SurroundConvolved[ff] == Threshold)
			Events_On[][1][rr]+=CurrentTrace_L[p+ff-(0.5*KernelLength)]
			nEvents+=1
			ff+=StimLength_pp
		endif
	endfor
	Events_On[][1][rr]/=nEvents
	// Off surround
	nEvents = 0
	for (ff=TriggerTimes[0]+KernelLength; ff<TriggerTimes[nTriggers-1]-KernelLength; ff+=1)
		if (SurroundConvolved[ff] < 1)
			Events_Off[][1][rr]+=CurrentTrace_L[p+ff-(0.5*KernelLength)]
			nEvents+=1
			ff+=StimLength_pp
		endif
	endfor
	Events_Off[][1][rr]/=nEvents
	// On full-field
	nEvents = 0
	for (ff=TriggerTimes[0]+KernelLength; ff<TriggerTimes[nTriggers-1]-KernelLength; ff+=1)
		if (CentreConvolved[ff] == Threshold && SurroundConvolved[ff] == Threshold)
			Events_On[][2][rr]+=CurrentTrace_L[p+ff-(0.5*KernelLength)]
			nEvents+=1
			ff+=StimLength_pp
		endif
	endfor
	Events_On[][2][rr]/=nEvents
	// Off full-field
	nEvents = 0
	for (ff=TriggerTimes[0]+KernelLength; ff<TriggerTimes[nTriggers-1]-KernelLength; ff+=1)
		if (CentreConvolved[ff] < 1 && SurroundConvolved[ff] < 1)
			Events_Off[][2][rr]+=CurrentTrace_L[p+ff-(0.5*KernelLength)]
			nEvents+=1
			ff+=StimLength_pp
		endif
	endfor
	Events_Off[][2][rr]/=nEvents
	
	print "Done with ROI", rr
endfor

end

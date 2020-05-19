#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// script calculates the mean response for each ROI (e.g. chirps, moving bars, sine, flashes)

function AverageLoops()

// Set parameters

variable SkipFrames = 100 // number of frames skipped after finding a trigger, bars: 10, chirp: 100

variable TriggerThreshold = 13000 // for trigger detection, arbitrary value
variable StimulusOffset = 34 // time between trigger and stimulus onset, setup specific, in ms
variable LineDuration = 2 // in ms
variable TriggersPerLoop = 1 // number of triggers per stimulus presentation


// Define input

string traces_name = "wDataCh0_ave_QA_Nor" // traces of all ROIs
string Roi_name = "ROIs" // name of ROI mask
string triggerchannel_name = "wDataCh2" // data channel to detect triggers
duplicate /o $traces_name Traces
duplicate /o $Roi_name Roi
duplicate /o $triggerchannel_name triggerchannel

// Get data dimensions

variable nF = dimsize($triggerchannel_name,2)
variable nL = dimsize($triggerchannel_name,1)
variable nX = dimsize($triggerchannel_name,0)
variable nRois = dimsize($traces_name,1)
variable rr, ff, xx, yy, ll

// Detrend data

duplicate /o Traces Traces_Detrended

variable r
for (r=0;r<nRois;r+=1)
     make /o/n=(nF) CurrentTrace = Traces[p][r]
     Smooth 2^15-1, CurrentTrace
     Traces_Detrended[][r]=Traces[p][r]-CurrentTrace[p]
endfor

killwaves CurrentTrace,Traces

// Get ROI position

setscale x, 0, nX, ROI
setscale y, 0, nL, ROI
GeometricCenter(ROI) // uses SARFIA built-in function 
wave GeoC

// Find triggers

make /o/n=1000 Triggertimes_Accurate = nan // saves accurate trigger position, line-precision
make /o/n=1000 Triggertimes_Frame = nan // saves triggers in frame-precision
variable nTriggers

for (ff=0; ff<nF; ff+=1)
	for (ll=0; ll<nL; ll+=1)
		for (xx=0; xx<nX; xx+=1)
			if (triggerchannel[xx][ll][ff]>TriggerThreshold)
				Triggertimes_Frame[nTriggers]=ff
				Triggertimes_Accurate[nTriggers]=(ff*nL)+ll
				ff+=SkipFrames
				ll=nL-1
				xx=0
				nTriggers+=1
			endif
		endfor
	endfor
endfor

variable nLoops = nTriggers/TriggersPerLoop
print "Stimulus loops: ",nLoops

Redimension/N=(nTriggers) Triggertimes_Frame, Triggertimes_Accurate

variable Snippet_Duration_pp = 16000 // number of points per stimulus presentation for interpolation, arbitrary (here set to 2 ms points, as stimulus is 32 s)
variable nFrames_Per_Cycle = Triggertimes_frame[3] - Triggertimes_frame[2] // number of frames per stimulus presentation

// Snippeting

make/o/n=(Snippet_duration_pp,nLoops,nRois) Responses = 0
make/o/n=(Snippet_duration_pp,nRois) ResponsesAverage = 0

make/o/n=(nFrames_Per_Cycle*nL) CurrentWave = 0
make/o/n=(nF*nL,nRois) Traces_Upsampled = 0

variable AddFrames = 6 // number of frames to add around current snippet 
make/o/n=(nFrames_Per_Cycle+AddFrames) CurrentTrace = 0

for (rr=0;rr<nRois;rr+=1)
	make/o/n=(nF) TraceToNormalize = Traces_Detrended[p][rr]
	for (ll=0;ll<(nLoops);ll+=1)
		variable LineOffset = Triggertimes_Accurate[ll] - (floor(triggertimes_accurate[ll]/nL)*nL) // number of line the trigger was found
		variable Response_Delay = GeoC[rr][1] - LineOffset // offset between trigger and ROI position, ROI-specific
		CurrentTrace[] = TraceToNormalize[p+(Triggertimes_Frame[ll]-(AddFrames/2))]
		Interpolate2/T=1/N=((nFrames_Per_Cycle+AddFrames)*nL)/Y=CurrentTrace_L CurrentTrace
		CurrentWave[]=CurrentTrace_L[p+nL*(AddFrames/2) + Response_Delay + StimulusOffset/LineDuration] // adjust for ROI position and stimulus offset
		Interpolate2/T=1/N=(Snippet_Duration_pp)/Y=CurrentWave_l CurrentWave
		Responses[][ll][rr]=CurrentWave_l[p]
		ResponsesAverage[][rr]+=CurrentWave_l[p]/nLoops
	endfor
endfor

killwaves CurrentWave, TraceToNormalize, CurrentTrace_L, CurrentTrace, CurrentWave_l

Setscale x,0,32,"s" Responses
Setscale x,0,32,"s" ResponsesAverage

// Estimate quality value (see Baden et al. Nature 2016)

make/o/n=(nRois) Quality = 0

for (rr=0;rr<nRois;rr+=1)
	variable VarianceOfMean = 0
	variable MeanOfVariance = 0
	for (ll=0;ll<(nLoops);ll+=1)
		make/o/n=(nF) CurrentTrace = Responses[p][ll][rr]
		WaveStats/Q CurrentTrace
		MeanOfVariance+=(V_sdev^2)/nLoops
	endfor
	make/o/n=(nF) CurrentTrace = ResponsesAverage[p][rr]
	WaveStats/Q CurrentTrace
	VarianceOfMean = V_sdev^2
	Quality[rr] = VarianceOfMean/MeanOfVariance
endfor

killwaves CurrentTrace

end
/*
Header file for IIR filter design functions
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#ifndef IIR_H
#define IIR_H

#include "constants.h"

// constants to specify filter response
#define IIR_LOWPASS  1
#define IIR_HIGHPASS 2
#define IIR_BANDPASS 4
#define IIR_BANDSTOP 8

// constants to specify filter design method
#define IIR_BUTTERWORTH 1
#define IIR_CHEBYCHEV1  2

// represents the cutoff frequency for a lowpass or highpass filter
typedef struct
{
	double dw; // transition bandwidth (rad)
	double wc; // cutoff frequency (rad)
} IIR_OneFreqSpec;

// represents the cutoff frequencies for a bandpass or bandstop filter
typedef struct
{
	double dw; // transition bandwidth (rad)
	double wl; // lower cutoff frequency (rad)
	double wu; // upper cutoff frequency (rad)
} IIR_TwoFreqSpec;

// represents the cutoff frequency / frequencies for any filter
typedef union
{
	IIR_OneFreqSpec lp; // lowpass spec
	IIR_OneFreqSpec hp; // highpass spec
	IIR_TwoFreqSpec bp; // bandpass spec
	IIR_TwoFreqSpec bs; // bandstop spec
} IIR_FreqSpec;

// represents a filter specification
typedef struct
{
	int response;      // any of IIR_LOWPASS, IIR_HIGHPASS, IIR_BANDPASS, or IIR_BANDSTOP
	int design;        // either IIR_BUTTERWORTH or IIR_CHEBYCHEV1
	double Rp;         // passband ripple (dB power)
	double As;         // stopband attenuation (dB power)
	IIR_FreqSpec freq; // cutoff frequency / cutoff frequencies
} IIR_FilterSpec;

// represents a second order section of a filter by its rational expression coefficients
typedef struct
{
	double b0;
	double b1;
	double b2;
	double a1;
	double a2;
} IIR_SectionSOS;

// represents a filter by second order sections
typedef struct
{
	int count;
	IIR_SectionSOS *sections;
	double gain;
} IIR_FilterSOS;

IIR_FilterSOS *IIR_DesignFilter(IIR_FilterSpec *spec);

#endif

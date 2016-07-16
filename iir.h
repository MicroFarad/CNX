/*
Header file for IIR filter design functions
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#ifndef IIR_H
#define IIR_H

#include <complex.h>

// represents the roots of a polynomial distinguishing between unpaired real and paired complex roots
typedef struct
{
	int n_unpaired;
	int n_paired;
	double *unpaired;
	double complex *paired;
} IIR_PairedRoots;

// a representation of a filter's transfer function by its poles and zeros and a constant
typedef struct
{
	double k;
	IIR_PairedRoots poles;
	IIR_PairedRoots zeros;
} IIR_PairedPZK;

// represents a polynomial with real coefficients of maximum order 2
typedef struct
{
	int order;
	double coeffs[3];
} IIR_RealQuadratic;

// represents a polynomial with complex coefficients of maximum order 2
typedef struct
{
	int order;
	double complex coeffs[3];
} IIR_ComplexQuadratic;

// a representation of a second order filter's transfer function by its rational expression
typedef struct
{
	IIR_RealQuadratic numerator;
	IIR_RealQuadratic denominator;
} IIR_QuadraticRationalExpression;

/*
// a second order filter section represented by coefficients
typedef struct
{
	double b0;
	double b1;
	double b2;
	double a1;
	double a2;
} IIR_SecondOrderSection;

// a first order filter section represented by coefficients
typedef struct
{
	double b0;
	double b1;
	double a1;
} IIR_FirstOrderSection;

// a representation of a filter by its first and second order sections
typedef struct
{
	int n_order1;
	int n_order2;
	IIR_FirstOrderSection *order1;
	IIR_SecondOrderSection *order2;
} IIR_SectionedFilter;

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

typedef struct
{
	double Wp;
	double Ws;
	double Rp;
	double As;
} IIR_AnalogDesignSpecification;
*/
typedef struct
{
	double Rp;
	double As;
	double Wp;
	double Ws;
} IIR_DesignSpec;

typedef struct
{
	int n_paired;
	int n_unpaired;
	double Wc;
} IIR_ButterworthSpec;

typedef struct
{
	int n_paired;
	int n_unpaired;
	double e;
	double Wc;
} IIR_ChebychevSpec;

void IIR_SubstituteRationalIntoPZK(IIR_PairedPZK *original, IIR_QuadraticRationalExpression *substitution, IIR_PairedPZK *result);

double complex IIR_ResponsePairedRoots(double complex x, IIR_PairedRoots *roots);

double complex IIR_ResponsePairedPZK(double complex x, IIR_PairedPZK *pzk);

void IIR_AnalogSpecificationButterworth(IIR_DesignSpec *dspec, IIR_ButterworthSpec *fspec);

void IIR_AnalogPrototypeButterworth(IIR_ButterworthSpec *fspec, IIR_PairedPZK *result);

void IIR_AnalogSpecificationChebychevI(IIR_DesignSpec *dspec, IIR_ChebychevSpec *fspec);

void IIR_AnalogPrototypeChebychevI(IIR_ChebychevSpec *fspec, IIR_PairedPZK *result);

void IIR_AnalogSpecificationChebychevII(IIR_DesignSpec *dspec, IIR_ChebychevSpec *fspec);

void IIR_AnalogPrototypeChebychevII(IIR_ChebychevSpec *fspec, IIR_PairedPZK *result);

#endif

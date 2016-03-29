/*
Header file for IIR filter design functions
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#ifndef IIR_H
#define IIR_H

#include <complex.h>
#include "constants.h"

// constants to specify filter response
#define IIR_LOWPASS  1
#define IIR_HIGHPASS 2
#define IIR_BANDPASS 3
#define IIR_BANDSTOP 4

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

// represents a second order section of a filter
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
	CNXIIR_SectionSOS *sections;
	double gain;
} IIR_FilterSOS;

// substitutes a rational expression as a second order section (substitution) into the independent variable of another second order section (original)
// produces, as a result, a second order section equivalent to the result of the substitution
// note that since this function is used for first order rational expressions, it is assumed that a2 = b2 = 0 for the substitution
void SubstituteFirstOrderExpression(IIR_FilterSOS *original, IIR_FilterSOS *substitution, IIR_FilterSOS *result)
{
	double dd0 = 1.0;
	double dd1 = 2.0 * substitution->a1;
	double dd2 = substitution->a1 * substitution->a1;
	double nd0 = substitution->b0;
	double nd1 = substitution->b0 * substitution->a1 + substitution->b1;
	double nd2 = substitution->b1 * substitution->a1;
	double nn0 = substitution->b0 * substitution->b0;
	double nn1 = 2.0 * substitution->b0 * substitution->b1;
	double nn2 = substitution->b1 * substitution->b1;
	double a0 = dd0 + original->a1 * nd0 + original->a2 * nn0;
	result->a1 = (dd1 + original->a1 * nd1 + original->a2 * nd1) / a0;
	result->a2 = (dd2 + original->a1 * nd2 + original->a2 * nd2) / a0;
	result->b0 = (original->b0 * dd0 + original->b1 * nd0 + original->b2 * nn0) / a0;
	result->b1 = (original->b0 * dd1 + original->b1 * nd1 + original->b2 * nn1) / a0;
	result->b2 = (original->b0 * dd2 + original->b1 * nd2 + original->b2 * nn2) / a0;
}

// substitutes a rational expression as a second order section (substitution) into the independent variable of another second order section (original)
// produces, as a result, two second order sections which, multiplied together, are equivalent to the result of the substitution
void SubstituteSecondOrderExpression(IIR_FilterSOS *original, IIR_FilterSOS *substitution, IIR_FilterSOS *result1, IIR_FilterSOS *result2)
{
	double r1, r2, s1, s2, a0_1, a0_2;

	r1 = original->a1;
	r2 = original->a2;
	s1 = (r1 + sqrt(r1 * r1 - 4.0 * r2)) / 2.0;
	s2 = r2 / s1;
	a0_1 = 1.0 + s1 * substituiton->b0;
	result1->a1 = (substitution->a1 + s1 * substitution->b1) / a0_1;
	result1->a2 = (substitution->a2 + s2 * substitution->b2) / a0_1;
	a0_2 = 1.0 + s2 * substitution->b0;
	result2->a1 = (substitution->a1 + s2 * substitution->b1) / a0_2;
	result2->a2 = (substitution->a2 + s2 * substitution->b2) / a0_2;

	r1 = original->b1 / original->b0;
	r2 = original->b2 / original->b0;
	s1 = (r1 + sqrt(r1 * r1 - 4.0 * r2)) / 2.0;
	s2 = r2 / s1;
	result1->b0 = original->b0 * (1.0 + s1 * substitution->b0) / a0_1;
	result1->b1 = original->b0 * (substitution->a1 + s1 * substitution->b1) / a0_1;
	result1->b2 = original->b0 * (substitution->a2 + s1 * substitution->b2) / a0_1;
	result2->b0 = (1.0 + s2 * substitution->b0) / a0_2;
	result2->b1 = (substitution->a1 + s2 * substitution->b1) / a0_2;
	result2->b2 = (substitution->a2 + s2 * substitution->b2) / a0_2;
}

// transforms an analog filter to a digital filter by the bilinear transform
void Bilinear(IIR_FilterSOS *analog, IIR_FilterSOS *digital)
{
	if(digital->count != analog->count) return;
	digital->gain = analog->gain;
	for(int k = 0; k < analog->count; k++)
	{
		double b0 = analog->sections[k].b0;
		double b1 = analog->sections[k].b1;
		double b2 = analog->sections[k].b2;
		double a1 = analog->sections[k].a1;
		double a2 = analog->sections[k].a2;
		double a0 = 1.0 + 2.0*a1 + 4.0*a2;
		digital->sections[k].b0 = (b0 + 2.0*b1 + 4.0*b2) / a0;
		digital->sections[k].b1 = (2.0*b0 - 8.0*b2) / a0;
		digital->sections[k].b2 = (b0 - 2.0*b1 + 4.0*b2) / a0;
		digital->sections[k].numorder = 2;
		digital->sections[k].a1 = (2.0 - 8.0*a2) / a0;
		digital->sections[k].a2 = (1.0 - 2.0*a1 + 4.0*a2) / a0;
		digital->sections[k].denorder = 2;
	}
}

// prewarp a frequency
double Prewarp(double w)
{
	return 2.0*tan(w/2.0);
}

// takes the specification for a digital filter with passband ripple Rp, stopband ripple As, passband up to wp, stopband beyond ws
// produces the order N and cutoff Wc of an analog butterworth prototype which can be bilinear transformed to the desired digital filter
void ButterworthSpecification(double Rp, double As, double wp, double ws, int *N, double *Wc)
{
	double Wp = Prewarp(wp);
	double Ws = Prewarp(ws);
	*N = ceil(log10((pow(10, Rp/10.0)-1.0)/(pow(10, As/10.0)-1.0))/(2.0*log10(Wp/Ws)));
	*Wc = ((Wp/pow(pow(10, Rp/10.0)-1.0, 1.0/(2.0*(*N))))+(Ws/pow(pow(10, As/10.0)-1.0, 1.0/(2.0*(*N)))))/2.0;
}

// takes the order N and cutoff Wc of an analog butterworth prototype and uses the bilinear transform to yeild a digital butterworth filter
void ButterworthPrototype(int N, double Wc, CNXIIR_FilterSOS *digital)
{
	int count = (N+1)/2;
	if(digital->count != count) return;
	CNXIIR_FilterSOS analog;
	analog.count = count;
	CNXIIR_SectionSOS analog_sections[count];
	analog.sections = analog_sections;
	for(int k = 0; k < analog.count; k++)
	{
		double t = TAU*(2.0*k+N+1.0)/(4.0*N);
		double r = Wc*cos(t);
		analog.sections[k].b0 = 1.0;
		analog.sections[k].b1 = 0.0;
		analog.sections[k].b2 = 0.0;
		if((k == nquads-1) && (N&1))
		{
			analog.sections[k].a1 = 1.0 / r;
			analog.sections[k].a2 = 0.0;
		}
		else
		{
			double i = Wc*sin(t);
			double a0 = r*r+i*i;
			analog.sections[k].a1 = -2.0*r/a0;
			analog.sections[k].a2 = 1.0/a0;
		}
	}
	analog.gain = 1.0;
	Bilinear(&analog, digital);
}

void DigitalLowpassTransform(wcp, wc, IIR_SectionSOS *transform)
{
	double a = sin((wcp-wc)/2.0)/sin((wcp+wc)/2.0);
	transform->b0 = -a;
	transform->b1 = 1.0;
	transform->b2 = 0.0;
	transform->a1 = -a;
	transform->a2 = 0.0;
}

void DigitalHighpassTransform(wcp, wc, IIR_SectionSOS *transform)
{
	double a = -cos((wcp+wc)/2.0)/cos((wcp-wc)/2.0);
	transform->b0 = -a;
	transform->b1 = -1.0;
	transform->b2 = 0.0;
	transform->a1 = a;
	transform->a2 = 0.0;
}

void DigitalBandpassTransform(wcp, wl, wu, IIR_SectionSOS *transform)
{
	double K = tan(wcp/2.0)/tan((wu-wl)/2.0);
	double B = cos((wu+wl)/2.0)/cos((wu-wl)/2.0);
	double a1 = -2.0*B*K/(K+1.0);
	double a2 = (K-1.0)/(K+1.0);
	transform->b0 = -a2;
	transform->b1 = a1;
	transform->b2 = -1.0;
	transform->a1 = -a1;
	transform->a2 = a2;
}

void DigitalBandstopTransform(wcp, wl, wu, IIR_SectionSOS *transform)
{
	double K = tan(wcp/2.0)*tan((wu-wl)/2.0);
	double B = cos((wu+wl)/2.0)/cos((wu-wl)/2.0);
	double a1 = -2.0*B/(K+1.0);
	double a2 = (K-1.0)/(K+1.0);
	transform->b0 = a2;
	transform->b1 = -a1;
	transform->b2 = 1.0;
	transform->a1 = -a1;
	transform->a2 = a2;
}

double complex SecondOrderSectionResponse(IIR_SectionSOS *section, double complex x)
{
	return (section->b0 + section->b1*x + section->b2*x*x)/(1.0 + section->a1*x + section->a2*x*x);
}

IIR_FilterSOS *DesignFilter(IIR_FilterSpec *spec)
{
	double wcp = 1.0;
	double wsp;
	IIR_SectionSOS transform;
	switch(spec->response)
	{
		case IIR_LOWPASS:
			DigitalLowpassTransform(wcp, spec->freq.lp.wc, &transform);
			wsp = carg(SecondOrderSectionResponse(&transform, cexp(-I*(spec->freq.lp.wc+spec->freq.lp.dw))));
			break;
		case IIR_HIGHPASS:
			DigitalHighpassTransform(wcp, spec->freq.hp.wc, &transform);
			wsp = carg(SecondOrderSectionResponse(&transform, cexp(-I*(spec->freq.hp.wc-spec->freq.hp.dw))));
			break;
		case IIR_BANDPASS:
			DigitalBandpassTransform(wcp, spec->freq.bp.wl, spec->freq.bp.wu, &transform);
			double wsp1 = carg(SecondOrderSectionResponse(&transform, cexp(-I*(spec->freq.bp.wu+spec->freq.bp.dw))));
			double wsp2 = carg(SecondOrderSectionResponse(&transform, cexp(-I*(spec->freq.bp.wl-spec->freq.bp.dw))));
			if(wsp1 < wsp2) wsp = wsp1;
			else wsp = wsp2;
			break;
		case IIR_BANDSTOP:
			DigitalBandpassTransform(wcp, spec->freq.bs.wl, spec->freq.bs.wu, &transform);
			double wsp1 = carg(SecondOrderSectionResponse(&transform, cexp(-I*(spec->freq.bs.wu-spec->freq.bs.dw))));
			double wsp2 = carg(SecondOrderSectionResponse(&transform, cexp(-I*(spec->freq.bs.wl+spec->freq.bs.dw))));
			if(wsp1 < wsp2) wsp = wsp1;
			else wsp = wsp2;
			break;
		default:
			return NULL;
	}
	int N;
	double Wc;
	switch(spec->design)
	{
		case IIR_BUTTERWORTH:
			ButterworthSpecification(spec->Rp, spec->As, wcp, ws, &N, &Wc);
			break;
		case IIR_CHEBYCHEV1:
			break;
		default:
			return NULL;
	}
	int count
	IIR_SectionSOS prototype[(N+1)/2];
	IIR_SectionSOS *transformed;
	switch(
	switch(spec->design)
	{
		case IIR_BUTTERWORTH:
			ButterworthPrototype(N, Wc, prototype);
			break;
		case IIR_CHEBYCHEV1:
			break;
	}
	IIR_SectionSOS transformed[
}

#endif

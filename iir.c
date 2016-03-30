/*
Source file for IIR filter design functions
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#include "iir.h"

#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

#define PLOTRES 1000

double complex SecondOrderSectionResponse(IIR_SectionSOS *section, double complex x)
{
	return (section->b0 + section->b1*x + section->b2*x*x)/(1.0 + section->a1*x + section->a2*x*x);
}

double complex FilterResponse(IIR_FilterSOS *filter, double w)
{
	double complex x = cexp(-I*w);
	double complex y = filter->gain;
	for(int i = 0; i < filter->count; i++) y *= SecondOrderSectionResponse(filter->sections+i, x);
	return y;
}

void BodeCSV(IIR_FilterSOS *filter, int points, int analog, double fmax, char *filename)
{
	FILE *fp = fopen(filename, "w");
	for(int n = 0; n < points; n++)
	{
		double f = fmax*n/points;
		double complex x;
		if(analog)
		{
			x = I*f;
		}
		else
		{
			x = cexp(-I*f);
		}
		double complex y = filter->gain;
		for(int i = 0; i < filter->count; i++) y *= SecondOrderSectionResponse(filter->sections+i, x);
		fprintf(fp, "%lf, %lf, %lf, %lf\n", f, cabs(y), 20*log10(cabs(y)), carg(y));
	}
	fclose(fp);
}

// substitutes a rational expression as a second order section (substitution) into the independent variable of another second order section (original)
// produces, as a result, a second order section equivalent to the result of the substitution
// note that since this function is used for first order rational expressions, it is assumed that a2 = b2 = 0 for the substitution
void SubstituteFirstOrderExpression(IIR_SectionSOS *original, IIR_SectionSOS *substitution, IIR_SectionSOS *result)
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
	result->a1 = (dd1 + original->a1 * nd1 + original->a2 * nn1) / a0;
	result->a2 = (dd2 + original->a1 * nd2 + original->a2 * nn2) / a0;
	result->b0 = (original->b0 * dd0 + original->b1 * nd0 + original->b2 * nn0) / a0;
	result->b1 = (original->b0 * dd1 + original->b1 * nd1 + original->b2 * nn1) / a0;
	result->b2 = (original->b0 * dd2 + original->b1 * nd2 + original->b2 * nn2) / a0;
}

void PolyRoots(double complex *polynomial, double complex *roots)
{
	double complex radical = csqrt(polynomial[1]*polynomial[1] - 4*polynomial[0]*polynomial[2]);
	roots[0] = (-polynomial[1]+radical)/(2*polynomial[0]);
	roots[1] = (-polynomial[1]-radical)/(2*polynomial[0]);
}

void PolyScale(double complex s, double complex *original, double complex *result)
{
	result[0] = s * original[0];
	result[1] = s * original[1];
	result[2] = s * original[2];
}

void PolyAdd(double complex *left, double complex *right, double complex *result)
{
	result[0] = left[0] + right[0];
	result[1] = left[1] + right[1];
	result[2] = left[2] + right[2];
}

void PolyFromRoots(double complex s, double complex r1, double complex r2, double complex *result)
{
	result[0] = s;
	result[1] = -s * (r1 + r2);
	result[2] = s * r1 * r2;
}

// substitutes a rational expression as a second order section (substitution) into the independent variable of another second order section (original)
// produces, as a result, two second order sections which, multiplied together, are equivalent to the result of the substitution
void SubstituteSecondOrderExpression(IIR_SectionSOS *original, IIR_SectionSOS *substitution, IIR_SectionSOS *result1, IIR_SectionSOS *result2)
{
	/*
	double r1, r2;
	double complex s1, s2;
	r1 = original->a1;
	r2 = original->a2;
	s1 = (r1 + csqrt(r1 * r1 - 4.0 * r2)) / 2.0;
	s2 = conj(s1);
	*/
	double complex ON[3], OD[3], SN[3], SD[3], P1[3], P2[3], RN1[3], RD1[3], RN2[3], RD2[3];
	double complex O[2];
	double complex R[4];
	double min;
	int minn, mini;
	int otrn, otri;
	double a0;

	ON[0] = original->b2;
	ON[1] = original->b1;
	ON[2] = original->b0;
	OD[0] = original->a2;
	OD[1] = original->a1;
	OD[2] = 1.0;
	SN[0] = substitution->b2;
	SN[1] = substitution->b1;
	SN[2] = substitution->b0;
	SD[0] = substitution->a2;
	SD[1] = substitution->a1;
	SD[2] = 1.0;
	
	PolyRoots(ON, O);
	PolyScale(-O[0], SN, P1);
	PolyAdd(P1, SD, P1);
	PolyScale(ON[0], P1, P1);
	PolyScale(-O[1], SN, P2);
	PolyAdd(P2, SD, P2);
	PolyRoots(P1, R+0);
	PolyRoots(P2, R+2);
	min = INFINITY;
	for(int n = 0; n < 4; n++)
	{
		double rn = creal(R[n]);
		for(int i = 3; i > n; i--)
		{
			double ri = creal(R[i]);
			double diff = fabs(rn-ri);
			if(diff < min)
			{
				min = diff;
				minn = n;
				mini = i;
			}
		}
	}
	for(otrn = 0; (otrn == minn)||(otrn == mini); otrn++);
	for(otri = 3; (otri == minn)||(otri == mini); otri--);
	PolyFromRoots(P1[0]*P2[0], R[minn], R[mini], RN1);
	PolyFromRoots(1.0, R[otrn], R[otri], RN2);
	
	PolyRoots(OD, O);
	PolyScale(-O[0], SN, P1);
	PolyAdd(P1, SD, P1);
	PolyScale(OD[0], P1, P1);
	PolyScale(-O[1], SN, P2);
	PolyAdd(P2, SD, P2);
	PolyRoots(P1, R+0);
	PolyRoots(P2, R+2);
	min = INFINITY;
	for(int n = 0; n < 4; n++)
	{
		double rn = creal(R[n]);
		for(int i = 3; i > n; i--)
		{
			double ri = creal(R[i]);
			double diff = fabs(rn-ri);
			if(diff < min)
			{
				min = diff;
				minn = n;
				mini = i;
			}
		}
	}
	for(otrn = 0; (otrn == minn)||(otrn == mini); otrn++);
	for(otri = 3; (otri == minn)||(otri == mini); otri--);
	PolyFromRoots(P1[0]*P2[0], R[minn], R[mini], RD1);
	PolyFromRoots(1.0, R[otrn], R[otri], RD2);

	a0 = creal(RD1[2]);
	result1->a1 = creal(RD1[1]) / a0;
	result1->a2 = creal(RD1[0]) / a0;
	result1->b0 = creal(RN1[2]) / a0;
	result1->b1 = creal(RN1[1]) / a0;
	result1->b2 = creal(RN1[0]) / a0;

	a0 = creal(RD2[2]);
	result2->a1 = creal(RD2[1]) / a0;
	result2->a2 = creal(RD2[0]) / a0;
	result2->b0 = creal(RN2[2]) / a0;
	result2->b1 = creal(RN2[1]) / a0;
	result2->b2 = creal(RN2[0]) / a0;
	/*
	printf("ROOTS:\n%f %f\n%f %f\n%f %f\n%f %f\n\n",
		creal(R[0]), cimag(R[0]), creal(R[1]), cimag(R[1]),
		creal(R[2]), cimag(R[2]), creal(R[3]), cimag(R[3]));
	printf("INDICES:\nminn: %d\nmini: %d\notrn: %d\notri: %d\n\n", minn, mini, otrn, otri);
	printf("POLY:\n%f %f\n%f %f\n%f %f\n\n%f %f\n%f %f\n%f %f\n\n",
		creal(P1[0]), cimag(P1[0]), creal(P1[1]), cimag(P1[1]), creal(P1[2]), cimag(P1[2]),
		creal(P2[0]), cimag(P2[0]), creal(P2[1]), cimag(P2[1]), creal(P2[2]), cimag(P2[2]));
	*/
	
	
	/*
	double r1, r2, s1, s2, a0_1, a0_2;

	r1 = original->a1;
	r2 = original->a2;
	s1 = (r1 + sqrt(r1 * r1 - 4.0 * r2)) / 2.0;
	s2 = r2 / s1;
	printf("%f %f %f %f\n", r1, r2, s1, s2);
	a0_1 = 1.0 + s1 * substitution->b0;
	result1->a1 = (substitution->a1 + s1 * substitution->b1) / a0_1;
	result1->a2 = (substitution->a2 + s1 * substitution->b2) / a0_1;
	a0_2 = 1.0 + s2 * substitution->b0;
	result2->a1 = (substitution->a1 + s2 * substitution->b1) / a0_2;
	result2->a2 = (substitution->a2 + s2 * substitution->b2) / a0_2;
	printf("%f %f\n", a0_1, a0_2);
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
	printf("%f %f %f %f %f\n", result1->b0, result1->b1, result1->b2, result1->a1, result1->a2);
	printf("%f %f %f %f %f\n", result2->b0, result2->b1, result2->b2, result2->a1, result2->a2);
	*/
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
		digital->sections[k].a1 = (2.0 - 8.0*a2) / a0;
		digital->sections[k].a2 = (1.0 - 2.0*a1 + 4.0*a2) / a0;
		printf("%f %f %f %f %f\n", digital->sections[k].b0, digital->sections[k].b1, digital->sections[k].b2, digital->sections[k].a1, digital->sections[k].a2);
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
	printf("wp: %f\nws: %f\n", wp, ws);
	double Wp = Prewarp(wp);
	double Ws = Prewarp(ws);
	printf("Wp: %f\nWs: %f\n", Wp, Ws);
	*N = ceil(log10((pow(10, Rp/10.0)-1.0)/(pow(10, As/10.0)-1.0))/(2.0*log10(Wp/Ws)));
	*Wc = ((Wp/pow(pow(10, Rp/10.0)-1.0, 1.0/(2.0*(*N))))+(Ws/pow(pow(10, As/10.0)-1.0, 1.0/(2.0*(*N)))))/2.0;
	printf("N: %d\nWc: %f\n", *N, *Wc);
}

// takes the order N and cutoff Wc of an analog butterworth prototype and uses the bilinear transform to yeild a digital butterworth filter
void ButterworthPrototype(int N, double Wc, IIR_FilterSOS *digital)
{
	int count = (N+1)/2;
	if(digital->count != count) return;
	IIR_FilterSOS analog;
	analog.count = count;
	IIR_SectionSOS analog_sections[count];
	analog.sections = analog_sections;
	for(int k = 0; k < analog.count; k++)
	{
		double t = TAU*(2.0*k+N+1.0)/(4.0*N);
		double r = Wc*cos(t);
		analog.sections[k].b0 = 1.0;
		analog.sections[k].b1 = 0.0;
		analog.sections[k].b2 = 0.0;
		if((k == analog.count-1) && (N&1))
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
	BodeCSV(&analog, PLOTRES, 1, 4*Wc, "analog_proto.csv");
	Bilinear(&analog, digital);
	BodeCSV(digital, PLOTRES, 0, TAU/2, "digital_proto.csv");
}

void DigitalLowpassTransform(double wcp, double wc, IIR_SectionSOS *transform)
{
	double a = sin((wcp-wc)/2.0)/sin((wcp+wc)/2.0);
	transform->b0 = -a;
	transform->b1 = 1.0;
	transform->b2 = 0.0;
	transform->a1 = -a;
	transform->a2 = 0.0;
}

void DigitalHighpassTransform(double wcp, double wc, IIR_SectionSOS *transform)
{
	double a = -cos((wcp+wc)/2.0)/cos((wcp-wc)/2.0);
	transform->b0 = -a;
	transform->b1 = -1.0;
	transform->b2 = 0.0;
	transform->a1 = a;
	transform->a2 = 0.0;
}

void DigitalBandpassTransform(double wcp, double wl, double wu, IIR_SectionSOS *transform)
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

void DigitalBandstopTransform(double wcp, double wl, double wu, IIR_SectionSOS *transform)
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

IIR_FilterSOS *IIR_DesignFilter(IIR_FilterSpec *spec)
{
	double wcp = 1.0;
	double wsp, wsp1, wsp2;
	IIR_SectionSOS transform;
	switch(spec->response)
	{
		case IIR_LOWPASS:
			DigitalLowpassTransform(wcp, spec->freq.lp.wc, &transform);
			wsp = carg(SecondOrderSectionResponse(&transform, cexp(I*(spec->freq.lp.wc+spec->freq.lp.dw))));
			break;
		case IIR_HIGHPASS:
			DigitalHighpassTransform(wcp, spec->freq.hp.wc, &transform);
			wsp = carg(SecondOrderSectionResponse(&transform, cexp(-I*(spec->freq.hp.wc-spec->freq.hp.dw))));
			break;
		case IIR_BANDPASS:
			DigitalBandpassTransform(wcp, spec->freq.bp.wl, spec->freq.bp.wu, &transform);
			wsp1 = carg(SecondOrderSectionResponse(&transform, cexp(I*(spec->freq.bp.wl+spec->freq.bp.dw))));
			wsp2 = carg(SecondOrderSectionResponse(&transform, cexp(I*(spec->freq.bp.wu-spec->freq.bp.dw))));
			printf("%f %f\n", wsp1, wsp2);
			if(wsp1 > wsp2) wsp = wsp1;
			else wsp = wsp2;
			break;
		case IIR_BANDSTOP:
			DigitalBandpassTransform(wcp, spec->freq.bs.wl, spec->freq.bs.wu, &transform);
			wsp1 = carg(SecondOrderSectionResponse(&transform, cexp(-I*(spec->freq.bs.wu-spec->freq.bs.dw))));
			wsp2 = carg(SecondOrderSectionResponse(&transform, cexp(-I*(spec->freq.bs.wl+spec->freq.bs.dw))));
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
			ButterworthSpecification(spec->Rp, spec->As, wcp, wsp, &N, &Wc);
			break;
		case IIR_CHEBYCHEV1:
			break;
		default:
			return NULL;
	}
	int count = (N+1)/2;
	IIR_SectionSOS prototype_sections[count];
	IIR_FilterSOS prototype;
	prototype.sections = prototype_sections;
	prototype.count = count;
	IIR_FilterSOS *filter = malloc(sizeof(IIR_FilterSOS));
	if(spec->response & (IIR_LOWPASS | IIR_HIGHPASS))
	{
		filter->count = count;
		filter->sections = malloc(count*sizeof(IIR_SectionSOS));
	}
	else
	{
		filter->count = 2*count;
		filter->sections = malloc(2*count*sizeof(IIR_SectionSOS));
	}
	filter->gain = 1.0;
	switch(spec->design)
	{
		case IIR_BUTTERWORTH:
			ButterworthPrototype(N, Wc, &prototype);
			break;
		case IIR_CHEBYCHEV1:
			break;
	}
	for(int n = 0; n < count; n++)
	{
		if(spec->response & (IIR_LOWPASS | IIR_HIGHPASS))
		{
			SubstituteFirstOrderExpression(prototype.sections+n, &transform, filter->sections+n);
		}
		else
		{
			SubstituteSecondOrderExpression(prototype.sections+n, &transform, filter->sections+2*n, filter->sections+2*n+1);
		}
	}
	BodeCSV(filter, PLOTRES, 0, TAU/2, "filter.csv");
	switch(spec->response)
	{
		case IIR_LOWPASS:
			printf("%f %f\n",
				-20*log10(cabs(FilterResponse(filter, spec->freq.lp.wc))),
				-20*log10(cabs(FilterResponse(filter, spec->freq.lp.wc+spec->freq.lp.dw))));
			break;
		case IIR_HIGHPASS:
			printf("%f %f\n",
				-20*log10(cabs(FilterResponse(filter, spec->freq.lp.wc-spec->freq.lp.dw))),
				-20*log10(cabs(FilterResponse(filter, spec->freq.lp.wc))));
			break;
		case IIR_BANDPASS:
			printf("CHECK: %f %f %f %f\n",
				-20*log10(cabs(FilterResponse(filter, spec->freq.bp.wl))),
				-20*log10(cabs(FilterResponse(filter, spec->freq.bp.wl+spec->freq.bp.dw))),
				-20*log10(cabs(FilterResponse(filter, spec->freq.bp.wu-spec->freq.bp.dw))),
				-20*log10(cabs(FilterResponse(filter, spec->freq.bp.wu))));
			break;
		case IIR_BANDSTOP:
			
			break;
		default:
			return NULL;
	}
	return filter;
}

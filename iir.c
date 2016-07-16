/*
Source file for IIR filter design functions
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#include "iir.h"
#include "constants.h"
#include "commacro.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PLOTRES 1000
/*
IIR_ComplexQuadratic *IIR_RealComplexQuadraticSum(IIR_RealQuadratic *left, IIR_ComplexQuadratic *right, IIR_ComplexQuadratic *result)
{
	result->order = left->order > right->order ? left->order : right->order;
	for(int n = 0; n <= result->order; n++)
		result->coeffs[n] = (n <= left->order ? left->coeffs[n] : 0) + (n <= right->order ? right->coeffs[n] : 0);
	return result;
}

IIR_RealQuadratic *IIR_RealQuadraticSum(IIR_RealQuadratic *left, IIR_RealQuadratic *right, IIR_RealQuadratic *result)
{
	result->order = left->order > right->order ? left->order : right->order;
	for(int n = 0; n <= result->order; n++)
		result->coeffs[n] = (n <= left->order ? left->coeffs[n] : 0) + (n <= right->order ? right->coeffs[n] : 0);
	return result;
}

IIR_ComplexQuadratic *IIR_RealQuadraticComplexScale(IIR_RealQuadratic *left, double complex right, IIR_ComplexQuadratic *result)
{
	result->order = left->order;
	for(int n = 0; n <= result->order; n++) result->coeffs[n] = left->coeffs[n] * right;
	return result;
}

IIR_RealQuadratic *IIR_RealQuadraticRealScale(IIR_RealQuadratic *left, double right, IIR_RealQuadratic *result)
{
	result->order = left->order;
	for(int n = 0; n <= result->order; n++) result->coeffs[n] = left->coeffs[n] * right;
	return result;
}
*/

#define square(x) ((x) * (x))

#define QuadraticDiscriminant(coeffs) ((coeffs)[1] * (coeffs)[1] - 4 * (coeffs)[2] * (coeffs)[0])

#define QuadraticRoot1(coeffs, radical) ((-(coeffs)[1] + radical) / (4 * (coeffs)[2]))

#define QuadraticRoot2(coeffs, radical) ((-(coeffs)[1] - radical) / (4 * (coeffs)[2]))

#define MonomialRoot(coeffs) (-(coeffs)[0] / (coeffs)[1])

// computes x * conj(x)
double ConjugateProduct(double complex x)
{
	return square(creal(x)) + square(cimag(x));
}

// adds the roots to a result paired roots polynomial due to substituting a rational expression into a single paired root term
// returns scaling coefficient
double PairedRootSubstitution(IIR_QuadraticRationalExpression *substitution, double complex root, IIR_PairedRoots *result)
{
	IIR_ComplexQuadratic term;
	term.order = MAX(substitution->numerator.order, substitution->denominator.order);
	for(int n = 0; n <= term.order; n++) term.coeffs[n] =
			(n <= substitution->numerator.order ? substitution->numerator.coeffs[n] : 0) -
			(n <= substitution->denominator.order ? root * substitution->denominator.coeffs[n] : 0);
	double complex radical;
	switch(term.order)
	{
		case 1:
			result->paired[result->n_paired++] = MonomialRoot(term.coeffs);
			return ConjugateProduct(term.coeffs[1]);
		case 2:
			// works under the assumption that the discriminant is not both real and positive
			radical = csqrt(QuadraticDiscriminant(term.coeffs));
			result->paired[result->n_paired++] = QuadraticRoot1(term.coeffs, radical);
			result->paired[result->n_paired++] = QuadraticRoot2(term.coeffs, radical);
			return ConjugateProduct(term.coeffs[2]);
		default:
			// something is broken
			return 0;
	}
}

// adds the roots to a result paired roots polynomial due to substituting a rational expression into a single paired root term
// returns scaling coefficient
double UnpairedRootSubstitution(IIR_QuadraticRationalExpression *substitution, double root, IIR_PairedRoots *result)
{
	IIR_RealQuadratic term;
	term.order = MAX(substitution->numerator.order, substitution->denominator.order);
	for(int n = 0; n <= term.order; n++) term.coeffs[n] =
			(n <= substitution->numerator.order ? substitution->numerator.coeffs[n] : 0) -
			(n <= substitution->denominator.order ? root * substitution->denominator.coeffs[n] : 0);
	double discriminant;
	switch(term.order)
	{
		case 1:
			result->unpaired[result->n_paired++] = MonomialRoot(term.coeffs);
			return term.coeffs[1];
		case 2:
			discriminant = QuadraticDiscriminant(term.coeffs);
			if(discriminant > 0)
			{
				double radical = sqrt(discriminant);
				result->unpaired[result->n_unpaired++] = QuadraticRoot1(term.coeffs, radical);
				result->unpaired[result->n_unpaired++] = QuadraticRoot2(term.coeffs, radical);
			}
			else result->paired[result->n_paired++] = QuadraticRoot1(term.coeffs, csqrt(discriminant));
			return term.coeffs[2];
		default:
			// something is broken
			return 0;
	}
}

// adds the roots of a quadratic multiple times to a paired roots polynomial
// returns scaling coefficient
double AddQuadraticRootsMultiple(int power, IIR_RealQuadratic *quadratic, IIR_PairedRoots *result)
{
	double discriminant;
	switch(quadratic->order)
	{
		case 1:
			result->unpaired[result->n_unpaired] = MonomialRoot(quadratic->coeffs);
			return pow(quadratic->coeffs[1], power);
		case 2:
			discriminant = QuadraticDiscriminant(quadratic->coeffs);
			if(discriminant > 0)
			{
				double radical = sqrt(discriminant);
				for(int n = 0; n < power; n++)
				{
					result->unpaired[result->n_unpaired++] = QuadraticRoot1(quadratic->coeffs, radical);
					result->unpaired[result->n_unpaired++] = QuadraticRoot2(quadratic->coeffs, radical);
				}
			}
			else result->paired[result->n_paired++] = QuadraticRoot1(quadratic->coeffs, csqrt(discriminant));
			return pow(quadratic->coeffs[2], power);
		default:
			// something is broken
			return 0;
	}
}

// substitutes a rational expression into a paired PZK
void IIR_SubstituteRationalIntoPZK(IIR_PairedPZK *original, IIR_QuadraticRationalExpression *substitution, IIR_PairedPZK *result)
{
	result->poles.n_unpaired = 0;
	result->zeros.n_unpaired = 0;
	result->poles.n_paired = 0;
	result->zeros.n_paired = 0;
	result->k = original->k;
	for(int n = 0; n < original->poles.n_unpaired; n++) result->k /= UnpairedRootSubstitution(substitution, original->poles.unpaired[n], &result->poles);
	for(int n = 0; n < original->zeros.n_unpaired; n++) result->k *= UnpairedRootSubstitution(substitution, original->zeros.unpaired[n], &result->zeros);
	for(int n = 0; n < original->poles.n_paired; n++) result->k /= PairedRootSubstitution(substitution, original->poles.paired[n], &result->poles);
	for(int n = 0; n < original->zeros.n_paired; n++) result->k *= PairedRootSubstitution(substitution, original->zeros.paired[n], &result->zeros);
	int imbalance = (original->zeros.n_unpaired + 2 * original->zeros.n_paired) - (original->poles.n_unpaired + 2 * original->poles.n_unpaired);
	if(imbalance > 0) result->k/= AddQuadraticRootsMultiple(imbalance, &substitution->denominator, &original->poles);
	else if(imbalance < 0) result->k*= AddQuadraticRootsMultiple(-imbalance, &substitution->denominator, &original->zeros);
}

double complex IIR_ResponsePairedRoots(double complex x, IIR_PairedRoots *roots)
{
	double complex response = 1.0;
	for(int n = 0; n < roots->n_paired; n++) response *= (x - roots->paired[n]) * (x - conj(roots->paired[n]));
	for(int n = 0; n < roots->n_unpaired; n++) response *= (x - roots->unpaired[n]);
	return response;
}

double complex IIR_ResponsePairedPZK(double complex x, IIR_PairedPZK *pzk)
{
	return pzk->k * IIR_ResponsePairedRoots(x, &pzk->zeros) / IIR_ResponsePairedRoots(x, &pzk->poles);
}

void IIR_AnalogSpecificationButterworth(IIR_DesignSpec *dspec, IIR_ButterworthSpec *fspec)
{
	int N = ceil(log10((pow(10.0, dspec->Rp / 10.0) - 1.0) / (pow(10.0, dspec->As / 10.0)-1.0)) / (2.0 * log10(dspec->Wp / dspec->Ws)));
	fspec->n_paired = N / 2;
	fspec->n_unpaired = N & 1;
	fspec->Wc = ((dspec->Wp/pow(pow(10.0, dspec->Rp / 10.0) - 1.0, 1.0 / (2.0 * N))) + (dspec->Ws / pow(pow(10.0, dspec->As / 10.0) - 1.0, 1.0 / (2.0 * N)))) / 2.0;
}

int ButterworthN(IIR_ButterworthSpec *fspec)
{
	return 2 * fspec->n_paired + fspec->n_unpaired;
}

void IIR_AnalogPrototypeButterworth(IIR_ButterworthSpec *fspec, IIR_PairedPZK *result)
{
	result->zeros.n_paired = 0;
	result->zeros.n_unpaired = 0;
	result->poles.n_paired = fspec->n_paired;
	result->poles.n_unpaired = fspec->n_unpaired;
	int N = ButterworthN(fspec);
	result->k = pow(-fspec->Wc, N);
	for(int k = 0; k < result->poles.n_paired; k++) result->poles.paired[k] = fspec->Wc * cexp(I * TAU * (2.0 * k + N + 1.0) / (4.0 * N));
	if(result->poles.n_unpaired) result->poles.unpaired[0] = -fspec->Wc;
}

int ChebychevN(IIR_ChebychevSpec *fspec)
{
	return 2 * fspec->n_paired + fspec->n_unpaired;
}

void IIR_AnalogSpecificationChebychevI(IIR_DesignSpec *dspec, IIR_ChebychevSpec *fspec)
{
	fspec->e = sqrt(pow(10.0, dspec->Rp / 10.0) - 1.0);
	int N = ceil(cacosh(sqrt(pow(10.0, dspec->As / 10.0) - 1.0) / fspec->e) / cacosh(dspec->Ws / dspec->Wp));
	fspec->n_paired = N / 2;
	fspec->n_unpaired = N & 1;
	fspec->Wc = dspec->Wp;
}

void IIR_AnalogPrototypeChebychevI(IIR_ChebychevSpec *fspec, IIR_PairedPZK *result)
{
	int N = ChebychevN(fspec);
	double phi = asinh(1.0 / fspec->e) / N;
	double a = sinh(phi);
	double b = cosh(phi);
	result->poles.n_paired = fspec->n_paired;
	result->poles.n_unpaired = fspec->n_unpaired;
	result->zeros.n_paired = 0;
	result->zeros.n_unpaired = 0;
	result->k = fspec->n_unpaired ? 1.0 : 1.0 / sqrt(pow(fspec->e, 2.0) + 1.0);
	for(int k = 0; k < result->poles.n_paired; k++)
	{
		double theta = (2.0 * k + N + 1.0) * TAU / (4.0 * N);
		double complex pole = fspec->Wc * (a * cos(theta) + I * b * sin(theta));
		result->poles.paired[k] = pole;
		result->k *= ConjugateProduct(pole);
	}
	if(result->poles.n_unpaired)
	{
		double pole = -fspec->Wc * a;
		result->poles.unpaired[0] = pole;
		result->k *= pole;
	}
}

void IIR_AnalogSpecificationChebychevII(IIR_DesignSpec *dspec, IIR_ChebychevSpec *fspec)
{
	fspec->e = 1.0 / sqrt(pow(10.0, dspec->As / 10.0) - 1.0);
	int N = ceil(cacosh(1.0 / (fspec->e * sqrt(pow(10.0, dspec->Rp / 10.0) - 1.0))) / cacosh(dspec->Ws / dspec->Wp));
	fspec->n_paired = N / 2;
	fspec->n_unpaired = N & 1;
	fspec->Wc = dspec->Ws;
}

void IIR_AnalogPrototypeChebychevII(IIR_ChebychevSpec *fspec, IIR_PairedPZK *result)
{
	int N = ChebychevN(fspec);
	double phi = asinh(1.0 / fspec->e) / N;
	double a = sinh(phi);
	double b = cosh(phi);
	result->poles.n_paired = fspec->n_paired;
	result->poles.n_unpaired = fspec->n_unpaired;
	result->zeros.n_paired = fspec->n_paired;
	result->zeros.n_unpaired = 0;
	result->k = 1.0;
	for(int k = 0; k < result->poles.n_paired; k++)
	{
		double theta = (2.0 * k + N + 1.0) * TAU / (4.0 * N);
		double complex pole = fspec->Wc / (a * cos(theta) + I * b * sin(theta));
		result->poles.paired[k] = pole;
		result->k *= ConjugateProduct(pole);
		double W = fspec->Wc / sin(theta);
		result->zeros.paired[k] = I * W;
		result->k /= -W * W;
	}
	if(result->poles.n_unpaired)
	{
		printf("unpaired");
		double pole = -fspec->Wc / a;
		result->poles.unpaired[0] = pole;
		result->k *= pole;
	}
	//result->k /= 10;
}

/*
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

void BodeCSV(IIR_SectionSOS *substitute, IIR_FilterSOS *filter, int points, int analog, double fmax, char *filename)
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
		if(substitute) x = SecondOrderSectionResponse(substitute, x);
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
/*
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
		double complex rn = R[n];
		for(int i = 3; i > n; i--)
		{
			double complex ri = R[i];
			double diff = cabs(rn-conj(ri));
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
	result2->b2 = creal(RN2[0]) / a0;*/
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
	*//*
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
	transform->a2 = a2;*/
	/*
	double a = sin((2.0*wcp-wu+wl)/4.0)/sin((2.0*wcp+wu-wl)/4.0);
	double B = cos((wl+wu)/2.0);
	transform->b0 = -1.0/a;
	transform->b1 = B*(1+a)/a;
	transform->b2 = -1.0;
	transform->a1 = -B*(1+a)/a;
	transform->a2 = 1.0/a;
	*/
	/*
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
			wsp1 = carg(SecondOrderSectionResponse(&transform, cexp(I*(spec->freq.bp.wu+spec->freq.bp.dw))));
			wsp2 = carg(SecondOrderSectionResponse(&transform, cexp(-I*(spec->freq.bp.wl-spec->freq.bp.dw))));
			printf("%f %f\n", wsp1, wsp2);
			if(wsp1 > wsp2) wsp = wsp1
			else wsp = wsp2;
			//wsp=0.3;
			break;
		case IIR_BANDSTOP:
			DigitalBandstopTransform(wcp, spec->freq.bs.wl, spec->freq.bs.wu, &transform);
			wsp1 = carg(SecondOrderSectionResponse(&transform, cexp(-I*(spec->freq.bs.wu-spec->freq.bs.dw))));
			wsp2 = carg(SecondOrderSectionResponse(&transform, cexp(I*(spec->freq.bs.wl+spec->freq.bs.dw))));
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
	BodeCSV(&transform, &prototype, PLOTRES, 0, TAU/2, "filter_substituted.csv");
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
	BodeCSV(NULL, filter, PLOTRES, 0, TAU/2, "filter.csv");
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
*/

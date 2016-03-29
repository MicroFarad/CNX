/*
Source file for CNX adaptive filtering
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#include "cnxadapt.h"
#include "prng.h"

void CNXADAPT_ConstructAdaptive(
	CNXADAPT_Adaptive *adaptive, int bufferlen, int coeffslen, kiss_fft_scalar alpha, kiss_fft_scalar beta,
	CNX_Connector *inbound_x, CNX_Connector *inbound_d,
	CNX_Connector *outbound_y, CNX_Connector *outbound_e, CNX_Connector *outbound_h)
{
	int nfft = 2*coeffslen + bufferlen - 2;
	if(nfft & 1) nfft++;
	adaptive->bufferlen = bufferlen;
	adaptive->coeffslen = coeffslen;
	adaptive->nfft = nfft;
	adaptive->fwdcfg = kiss_fftr_alloc(nfft, 0, NULL, NULL);
	adaptive->backcfg = kiss_fftr_alloc(nfft, 1, NULL, NULL);
	adaptive->in_x = NULL;
	adaptive->in_d = NULL;
	adaptive->inbound_x = inbound_x;
	adaptive->inbound_d = inbound_d;
	adaptive->outbound_y = outbound_y;
	adaptive->outbound_e = outbound_e;
	adaptive->outbound_h = outbound_h;
	adaptive->alpha = alpha;
	adaptive->beta = beta;
	adaptive->coefficients = malloc(nfft*sizeof(kiss_fft_scalar));
	adaptive->history = malloc(nfft*sizeof(kiss_fft_scalar));
	unsigned long long smstate = 1337;
	unsigned long long xorstate[2];
	xorstate[0] = PRNG_splitmix64(&smstate);
	xorstate[1] = PRNG_splitmix64(&smstate);
	kiss_fft_scalar s = 1/((CNX_Scalar)0x8000000000000000ull);
	adaptive->tmp1 = malloc((nfft/2+1)*sizeof(kiss_fft_cpx));
	adaptive->tmp2 = malloc((nfft/2+1)*sizeof(kiss_fft_cpx));
	adaptive->tmp3 = malloc((nfft/2+1)*sizeof(kiss_fft_cpx));
	adaptive->tmp4 = malloc(nfft*sizeof(kiss_fft_scalar));
	adaptive->tmp5 = malloc(nfft*sizeof(kiss_fft_scalar));
	printf(">>>>>> %d %p\n",nfft,adaptive->tmp5);
	for(int n = 0; n < nfft; n++)
	{
		kiss_fft_scalar u, v, w;
		do
		{
			u = (CNX_Scalar)PRNG_xorshift128plus(xorstate)*s-1;
			v = (CNX_Scalar)PRNG_xorshift128plus(xorstate)*s-1;
			w = u * u + v * v;
		}
		while(w >= 1.0);
		w = sqrt(-2.0 * log(w) / w);
		adaptive->coefficients[n] = u*w;
		adaptive->history[n] = 0;
		adaptive->tmp4[n] = 0;
		adaptive->tmp5[n] = 0;
	}
	for(int n = 0; n <= nfft / 2; n++)
	{
		adaptive->tmp1[n].r = 0;
		adaptive->tmp1[n].i = 0;
		adaptive->tmp2[n].r = 0;
		adaptive->tmp2[n].i = 0;
		adaptive->tmp3[n].r = 0;
		adaptive->tmp3[n].i = 0;
	}
}

int CNXADAPT_ProcessAdaptive(void *adaptive, int count)
{
	int n;
	CNXADAPT_Adaptive *adapt = (CNXADAPT_Adaptive*)adaptive;
	CNX_Buffer *in_x = adapt->in_x;
	CNX_Buffer *in_d = adapt->in_d;
	adapt->in_x = NULL;
	adapt->in_d = NULL;
	CNX_Buffer *out;
	int bufferlen = adapt->bufferlen;
	int coeffslen = adapt->coeffslen;
	int nfft = adapt->nfft;
	kiss_fft_scalar alpha = adapt->alpha;
	kiss_fft_scalar beta = adapt->beta;
	kiss_fft_cpx *tmp1 = adapt->tmp1;
	kiss_fft_cpx *tmp2 = adapt->tmp2;
	kiss_fft_cpx *tmp3 = adapt->tmp3;
	kiss_fft_scalar *x = adapt->history + coeffslen - 1;
	kiss_fft_scalar *e = adapt->tmp5;
	kiss_fft_scalar *y = adapt->tmp4 + coeffslen - 1;
	kiss_fft_scalar *a = adapt->tmp4 + bufferlen - 1;
	kiss_fft_scalar *h = adapt->coefficients;
	kiss_fft_scalar scl = 1.0/(nfft*nfft);
	for(n = 0; n < count; n++)
	{
		if(!in_x) in_x = CNX_TakeQueue(adapt->inbound_x);
		if(!in_d) in_d = CNX_TakeQueue(adapt->inbound_d);
		if(!in_x || !in_d)
		{
			adapt->in_x = in_x;
			adapt->in_d = in_d;
			break;
		}
		memcpy(adapt->history, adapt->history+bufferlen, (coeffslen-1)*sizeof(kiss_fft_scalar));
		memcpy(adapt->history+coeffslen-1, in_x->data, bufferlen*sizeof(kiss_fft_scalar));
		CNX_PoolBuffer(adapt->inbound_x, in_x);
		kiss_fftr(adapt->fwdcfg, adapt->history, adapt->tmp1);
		kiss_fftr(adapt->fwdcfg, adapt->coefficients, adapt->tmp2);
		for(int i = 0; i <= nfft/2; i++)
		{
			tmp3[i].r = tmp1[i].r * tmp2[i].r - tmp1[i].i * tmp2[i].i;
			tmp3[i].i = tmp1[i].r * tmp2[i].i + tmp1[i].i * tmp2[i].r;
		}
		kiss_fftri(adapt->backcfg, tmp3, adapt->tmp4);
		for(int i = 0; i < bufferlen; i++) y[i] *= scl;
		if(adapt->outbound_y)
		{
			out = CNX_TakePool(adapt->outbound_y);
			memcpy(out->data, y, bufferlen);
			CNX_QueueBuffer(adapt->outbound_y, out);
		}
		kiss_fft_scalar *d = in_d->data;
		for(int i = 0; i < bufferlen; i++) e[i] = d[i] - y[i];
		CNX_PoolBuffer(adapt->inbound_d, in_d);
		if(adapt->outbound_e)
		{
			out = CNX_TakePool(adapt->outbound_e);
			memcpy(out->data, e, bufferlen);
			CNX_QueueBuffer(adapt->outbound_e, out);
		}
		kiss_fft_scalar p = 0;
		//for(int i = 1-coeffslen; i < 0; i++) p += x[i] * x[i];
		for(int i = 0; i < bufferlen; i++)
		{
			//p += x[i] * x[i];
			e[i] *= alpha / bufferlen;// / p;
			//p -= x[i-coeffslen+1] * x[i-coeffslen+1];
		}
		kiss_fftr(adapt->fwdcfg, e, adapt->tmp2);
		for(int i = 0; i <= nfft/2; i++)
		{
			tmp3[i].r = tmp1[i].r * tmp2[i].r - tmp1[i].i * tmp2[i].i;
			tmp3[i].i = tmp1[i].r * tmp2[i].i + tmp1[i].i * tmp2[i].r;
		}
		kiss_fftri(adapt->backcfg, tmp3, adapt->tmp4);
		for(int i = 0; i < coeffslen; i++) h[i] = beta*h[i] + a[i]*scl;
		if(adapt->outbound_h)
		{
			out = CNX_TakePool(adapt->outbound_h);
			memcpy(out->data, h, coeffslen);
			CNX_QueueBuffer(adapt->outbound_h, out);
		}
		in_x = NULL;
		in_d = NULL;
	}
	return n;
}

void CNXADAPT_DestroyAdaptive(CNXADAPT_Adaptive *adaptive)
{
	free(adaptive->fwdcfg);
	free(adaptive->backcfg);
	if(adaptive->in_x) CNX_PoolBuffer(adaptive->inbound_x, adaptive->in_x);
	if(adaptive->in_d) CNX_PoolBuffer(adaptive->inbound_d, adaptive->in_d);
	free(adaptive->coefficients);
	free(adaptive->history);
	free(adaptive->tmp1);
	free(adaptive->tmp2);
	free(adaptive->tmp3);
	free(adaptive->tmp4);
	free(adaptive->tmp5);
}

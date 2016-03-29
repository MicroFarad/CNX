/*
Header file for CNX adaptive filtering
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#ifndef CNXADAPT_H
#define CNXADAPT_H

#include "cnx.h"
#include "kiss_fftr.h"

typedef struct
{
	kiss_fftr_cfg fwdcfg;
	kiss_fftr_cfg backcfg;
	CNX_Connector *inbound_x;
	CNX_Connector *inbound_d;
	CNX_Connector *outbound_y;
	CNX_Connector *outbound_e;
	CNX_Connector *outbound_h;
	CNX_Buffer *in_x;
	CNX_Buffer *in_d;
	kiss_fft_scalar *coefficients;
	kiss_fft_scalar *history;
	kiss_fft_cpx *tmp1;
	kiss_fft_cpx *tmp2;
	kiss_fft_cpx *tmp3;
	kiss_fft_scalar *tmp4;
	kiss_fft_scalar *tmp5;
	int bufferlen;
	int coeffslen;
	int nfft;
	kiss_fft_scalar alpha;
	kiss_fft_scalar beta;
} CNXADAPT_Adaptive;


void CNXADAPT_ConstructAdaptive(
	CNXADAPT_Adaptive *adaptive, int bufferlen, int coeffslen, kiss_fft_scalar alpha, kiss_fft_scalar beta,
	CNX_Connector *inbound_x, CNX_Connector *inbound_d,
	CNX_Connector *outbound_y, CNX_Connector *outbound_e, CNX_Connector *outbound_h);

int CNXADAPT_ProcessAdaptive(void *adaptive, int count);

void CNXADAPT_DestroyAdaptive(CNXADAPT_Adaptive *adaptive);

#endif

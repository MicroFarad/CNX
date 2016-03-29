/*
Header file for CNX wrapper of Mark Borgerding's KISS FFT
Copyright (C) Kyle Gagner 2016
All rights reserved
*/

#ifndef CNXKISS_H
#define CNXKISS_H

#include "cnx.h"
#include "kiss_fft.h"

typedef struct
{
	int length;
	kiss_fft_cfg cfg;
} CNXKISS_Prototype;

typedef struct
{
	CNXKISS_Prototype *prototype;
	CNX_Connector *inbound;
	CNX_Connector *outbound;
} CNXKISS_Instance;

void CNXKISS_ConstructPrototype(CNXKISS_Prototype *prototype, int length, int inverse);

void CNXKISS_ConstructInstance(CNXKISS_Prototype *prototype, CNXKISS_Instance *instance, CNX_Connector *inbound, CNX_Connector *outbound);

int CNXKISS_Process(void *instance, int count);

void CNXKISS_DestroyPrototype(CNXKISS_Prototype *prototype);

#endif

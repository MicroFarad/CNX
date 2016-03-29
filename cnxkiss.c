/*
Source file for CNX wrapper of Mark Borgerding's KISS FFT
Copyright (C) Kyle Gagner 2016
All rights reserved
*/

#include "cnxkiss.h"

void CNXKISS_ConstructPrototype(CNXKISS_Prototype *prototype, int length, int inverse)
{
	prototype->length = length;
	prototype->cfg = kiss_fft_alloc(length, inverse, NULL, NULL);
}

void CNXKISS_ConstructInstance(CNXKISS_Prototype *prototype, CNXKISS_Instance *instance, CNX_Connector *inbound, CNX_Connector *outbound)
{
	instance->prototype = prototype;
	instance->inbound = inbound;
	instance->outbound = outbound;
}

int CNXKISS_Process(void *instance, int count)
{
	int n;
	CNX_Connector *inbound = ((CNXKISS_Instance*)instance)->inbound;
	CNX_Connector *outbound = ((CNXKISS_Instance*)instance)->outbound;
	kiss_fft_cfg cfg = ((CNXKISS_Instance*)instance)->prototype->cfg;
	CNX_Buffer *in;
	for(n = 0; (n < count) && (in = CNX_TakeQueue(inbound)); n++)
	{
		CNX_Buffer *out = CNX_TakePool(outbound);
		kiss_fft(cfg, in->data, out->data);
		CNX_QueueBuffer(outbound, out);
		CNX_PoolBuffer(inbound, in);
	}
	return n;
}

void CNXKISS_DestroyPrototype(CNXKISS_Prototype *prototype)
{
	free(prototype->cfg);
}

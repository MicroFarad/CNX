/*
Source file for basic CNX algorithms
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#include <string.h>
#include "cnxbasic.h"

int CNXBASIC_ProcessAmplifier(void *amplifier, int count)
{
	int n;
	CNX_Connector *inbound = ((CNXBASIC_Amplifier*)amplifier)->inbound;
	CNX_Connector *outbound = ((CNXBASIC_Amplifier*)amplifier)->outbound;
	CNX_Scalar gain = ((CNXBASIC_Amplifier*)amplifier)->gain;
	CNX_Buffer *in;
	for(n = 0; (n < count) && (in = CNX_TakeQueue(inbound)); n++)
	{
		CNX_Buffer *out = CNX_TakePool(outbound);
		CNX_Scalar *inp = in->data;
		CNX_Scalar *outp = out->data;
		for(int i = inbound->bufferlength; i; i--)
		{
			*outp = *inp * gain;
			inp++;
			outp++;
		}
		CNX_QueueBuffer(outbound, out);
		CNX_PoolBuffer(inbound, in);
	}
	return n;
}

int CNXBASIC_ProcessConsumer(void *inbound, int count)
{
	int n;
	CNX_Buffer *in;
	for(n = 0; (n < count) && (in = CNX_TakeQueue(inbound)); n++) CNX_PoolBuffer(inbound, in);
	return n;
}

int CNXBASIC_ProcessRoundRobin(void *splitter, int count)
{
	int n;
	CNX_Connector *inbound = ((CNXBASIC_Splitter*)splitter)->inbound;
	CNX_Connector **outbounds = ((CNXBASIC_Splitter*)splitter)->outbounds;
	CNX_Buffer *buf;
	for(n = 0; (n < count) && (buf = CNX_TakeQueue(inbound)); n++)
	{
		CNX_Connector **outbound;
		for(outbound = outbounds; *outbound; outbound++)
		{
			CNX_QueueBuffer(*outbound, buf);
			buf = CNX_TakePool(*outbound);
		}
		CNX_PoolBuffer(inbound, buf);
	}
	return n;
}

int CNXBASIC_ProcessCopySplit(void *splitter, int count)
{
	int n;
	CNX_Connector *inbound = ((CNXBASIC_Splitter*)splitter)->inbound;
	CNX_Connector **outbounds = ((CNXBASIC_Splitter*)splitter)->outbounds;
	CNX_Buffer *in;
	for(n = 0; (n < count) && (in = CNX_TakeQueue(inbound)); n++)
	{
		CNX_Connector **outbound;
		for(outbound = outbounds; *outbound; outbound++)
		{
			CNX_Buffer *out = CNX_TakePool(*outbound);
			memcpy(out->data, in->data, inbound->bufferwidth * inbound->bufferlength);
			CNX_QueueBuffer(*outbound, out);
		}
		CNX_PoolBuffer(inbound, in);
	}
	return n;
}

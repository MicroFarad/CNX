/*
Source file for CNX wrapper of Port Audio
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#include "cnxpa.h"
#include <string.h>

int pacallback(const void* inbuffer, void* outbuffer, unsigned long frames, const PaStreamCallbackTimeInfo* timeinfo, PaStreamCallbackFlags status, void* data)
{
	CNX_Buffer *buffer;
	buffer = CNX_TakeQueue(((CNXPA_Stream*)data)->inbound);
	if(buffer) memcpy(outbuffer, buffer->data, frames*sizeof(float));
	else memset(outbuffer, 0, frames*sizeof(float));
	CNX_PoolBuffer(((CNXPA_Stream*)data)->inbound, buffer);
	buffer = CNX_TakePool(((CNXPA_Stream*)data)->outbound);
	memcpy(buffer->data, inbuffer, frames*sizeof(float));
	CNX_QueueBuffer(((CNXPA_Stream*)data)->outbound, buffer);
	return paNoError;
}

void CNXPA_Initialize(void)
{
	Pa_Initialize();
}

void CNXPA_Terminate(void)
{
	Pa_Terminate();
}

void CNXPA_ConstructStream(CNXPA_Stream *stream, CNX_Connector *inbound, CNX_Connector *outbound, int fs, int bufferlength)
{
	stream->inbound = inbound;
	stream->outbound = outbound;
	Pa_OpenDefaultStream(&stream->stream, 1, 1, paFloat32, fs, bufferlength, pacallback, stream);
}

void CNXPA_StartStream(CNXPA_Stream *stream)
{
	Pa_StartStream(stream->stream);
}

void CNXPA_StopStream(CNXPA_Stream *stream)
{
	Pa_StopStream(stream->stream);
}

void CNXPA_DestroyStream(CNXPA_Stream *stream)
{
	Pa_CloseStream(stream->stream);
}

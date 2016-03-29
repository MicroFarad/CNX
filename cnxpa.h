/*
Header file for CNX wrapper of Port Audio
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#ifndef CNXPA_H
#define CNXPA_H

#include "cnx.h"
#include "portaudio.h"

#define CNXPA_8K    8000
#define CNXPA_22K05 22050
#define CNXPA_44K1  44100
#define CNXPA_96K   96000

typedef struct
{
	PaStream *stream;
	CNX_Connector *inbound;
	CNX_Connector *outbound;
} CNXPA_Stream;

void CNXPA_Initialize(void);

void CNXPA_Terminate(void);

void CNXPA_ConstructStream(CNXPA_Stream *stream, CNX_Connector *inbound, CNX_Connector *outbound, int fs, int bufferlength);

void CNXPA_StartStream(CNXPA_Stream *stream);

void CNXPA_StopStream(CNXPA_Stream *stream);

void CNXPA_DestroyStream(CNXPA_Stream *stream);

#endif

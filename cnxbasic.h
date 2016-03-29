/*
Header file for basic CNX algorithms
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#ifndef CNXBASIC_H
#define CNXBASIC_H

#include "cnx.h"

typedef struct
{
	CNX_Scalar gain;
	CNX_Connector *inbound;
	CNX_Connector *outbound;
} CNXBASIC_Amplifier;

typedef struct
{
	CNX_Connector *inbound;
	CNX_Connector **outbounds;
} CNXBASIC_Splitter;

int CNXBASIC_ProcessAmplifier(void *amplifier, int count);

int CNXBASIC_ProcessConsumer(void *inbound, int count);

int CNXBASIC_ProcessRoundRobin(void *splitter, int count);

int CNXBASIC_ProcessCopySplit(void *splitter, int count);

#endif

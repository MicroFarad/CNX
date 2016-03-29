/*
Source file for CNX of xorshift128+ and splitmix64 PRNGs
This is a derivative work adapted in 2016 from Sebastiano Vigna's code at http://xorshift.di.unimi.it/ by Kyle Gagner
Copyright (C) 2016 Kyle Gagner
All rights reserved

Sebastiano Vigna's original implementation is in the public domain and may be obtained at the above url
*/

#ifndef CNXRAND_H
#define CNXRAND_H

#include "cnx.h"

typedef struct
{
	unsigned long long state[2];
	CNX_Connector *outbound;
} CNXRAND_Generator;

void CNXRAND_ConstructGenerator(CNXRAND_Generator *generator, CNX_Connector *outbound, unsigned long long seed);

int CNXRAND_Generate(void *generator, int count);

#endif

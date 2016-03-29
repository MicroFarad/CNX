/*
Source file for CNX gaussian white noise generator
Copyright (C) 2016 Kyle Gagner
All rights reserved
*/

#include <math.h>
#include "cnxrand.h"
#include "prng.h"

void CNXRAND_ConstructGenerator(CNXRAND_Generator *generator, CNX_Connector *outbound, unsigned long long seed)
{
	generator->outbound = outbound;
	generator->state[0] = PRNG_splitmix64(&seed);
	generator->state[1] = PRNG_splitmix64(&seed);
}

int CNXRAND_Generate(void *generator, int count)
{
	CNX_Connector *outbound = ((CNXRAND_Generator*)generator)->outbound;
	unsigned long long *state = ((CNXRAND_Generator*)generator)->state;
	CNX_Scalar s = 1/((CNX_Scalar)0x8000000000000000ull);
	for(int n = 0; n < count; n++)
	{
		CNX_Buffer *out = CNX_TakePool(outbound);
		CNX_Scalar *p = out->data;
		for(int i = outbound->bufferlength; i; i--)
		{
			CNX_Scalar u, v, w;
			do
			{
				u = (CNX_Scalar)PRNG_xorshift128plus(state)*s-1;
				v = (CNX_Scalar)PRNG_xorshift128plus(state)*s-1;
				w = u * u + v * v;
			}
			while(w >= 1.0);
			w = sqrt(-2.0 * log(w) / w);
			*p = u*w;
			p++;
			i--;
			if(i) *p = v*w;
			p++;
		}
		CNX_QueueBuffer(outbound, out);
	}
	return count;
}

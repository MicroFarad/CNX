/*
Header file for xorshift128+ and plitmix64 PRNGs
This is a derivative work adapted in 2016 from Sebastiano Vigna's code at http://xorshift.di.unimi.it/ by Kyle Gagner
Copyright (C) 2016 Kyle Gagner
All rights reserved

Sebastiano Vigna's original implementation is in the public domain and may be obtained at the above url
*/

#ifndef PRNG_H
#define PRNG_H

unsigned long long PRNG_splitmix64(unsigned long long *state);

unsigned long long PRNG_xorshift128plus(unsigned long long *state);

#endif

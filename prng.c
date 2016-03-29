/*
Source file for xorshift128+ and plitmix64 PRNGs
This is a derivative work adapted in 2016 from Sebastiano Vigna's code at http://xorshift.di.unimi.it/ by Kyle Gagner
Copyright (C) 2016 Kyle Gagner
All rights reserved

Sebastiano Vigna's original implementation is in the public domain and may be obtained at the above url
*/

unsigned long long PRNG_splitmix64(unsigned long long *state)
{
	unsigned long long z = ((*state) += 0x9E3779B97F4A7C15ull);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
	return z ^ (z >> 31);
}

unsigned long long PRNG_xorshift128plus(unsigned long long *state)
{
	unsigned long long s1 = state[0];
	unsigned long long s0 = state[1];
	state[0] = s0;
	s1 ^= s1 << 23;
	state[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
	return state[1] + s0;
}

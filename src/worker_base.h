#pragma once

#include "defs.h"
#include "params.h"
#include <string>

#include "xmmintrin.h"
#include <nmmintrin.h>
#include <immintrin.h>
//#include <intrin.h>

using namespace std;

class CWorkerBase
{
protected:
	CParams& params;

	int codes[256];

	int hts_mask;
	uint32_t htl_size;
	uint32_t htl_mask;

	const double htl_max_fill_factor = 0.1;

	vector<CFactor> v_parsing;

	vector<pair<int64_t, int64_t>> v_kmers_rl, v_kmers_rs;
	vector<pair<int64_t, int64_t>> v_kmers_dl, v_kmers_ds;

	int hash_mm(uint64_t x, int mask)
	{
		x ^= x >> 33;
		x *= 0xff51afd7ed558ccdLL;
		x ^= x >> 33;
		x *= 0xc4ceb9fe1a85ec53LL;
		x ^= x >> 33;

		return (int)(x & mask);
	}

	int lzcnt(uint64_t x)
	{
		x |= x >> 1;
		x |= x >> 2;
		x |= x >> 4;
		x |= x >> 8;
		x |= x >> 16;
		x |= x >> 32;

		return (int)_mm_popcnt_u64(~x);
	}

	int lzcnt32(uint32_t x)
	{
		x |= x >> 1;
		x |= x >> 2;
		x |= x >> 4;
		x |= x >> 8;
		x |= x >> 16;

		return (int)_mm_popcnt_u32(~x);
	}



public:
	CWorkerBase(CParams& params) : params(params)
	{
		fill_n(codes, 256, 4);
		codes['A'] = 0;
		codes['C'] = 1;
		codes['G'] = 2;
		codes['T'] = 3;

		hts_mask = (int)(1u << (2 * (params.min_distant_match_len - params.min_match_len))) - 1;
	};

};
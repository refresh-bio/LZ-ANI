#pragma once

#include "defs.h"

#include <cstdint>
#include <vector>
#include <utility>

#include "xmmintrin.h"
#include <nmmintrin.h>
#include <immintrin.h>

class BaseWorker {
public:

	BaseWorker();

	bool load_file(const string &file_name, seq_t &seq, uint32_t &n_parts, int separator);

	void duplicate_rev_comp(seq_t &seq);

	void export_parsing();

protected:

	int codes[256];
	int hts_mask;

	std::vector<CFactor> v_parsing;

	void prepare_kmers(std::vector<std::pair<int64_t, int>> &v_kmers, const seq_t &seq, int len, bool store_all = false);

	// ****************************************************************************
	// !!! To moze byc szybsze jesli CPU ma instrukcje _lzcnt. Ona niestety nie zawsze jest obecna.
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

	// ****************************************************************************
	// !!! To moze byc szybsze jesli CPU ma instrukcje _lzcnt. Ona niestety nie zawsze jest obecna.
	int lzcnt32(uint32_t x)
	{
		x |= x >> 1;
		x |= x >> 2;
		x |= x >> 4;
		x |= x >> 8;
		x |= x >> 16;

		return (int)_mm_popcnt_u32(~x);
	}

	// ****************************************************************************
	int hash_mm(uint64_t x, int mask)
	{
		x ^= x >> 33;
		x *= 0xff51afd7ed558ccdLL;
		x ^= x >> 33;
		x *= 0xc4ceb9fe1a85ec53LL;
		x ^= x >> 33;

		return (int)(x & mask);
	}

};
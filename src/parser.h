#pragma once
#include <vector>
#include <array>

#include "xmmintrin.h"
#include <nmmintrin.h>
#include <immintrin.h>

#include "defs.h"
#include "seq_reservoir.h"

using namespace std;

class CParser
{
	CParams2 params;

	seq_t seq_ref;
	seq_t seq_data;
	uint32_t n_ref_seqs;
	uint32_t n_data_seqs;

	seq_view sv_data;

	vector<CFactor> v_parsing;

	vector<pair<int64_t, int64_t>> v_kmers_ref_long, v_kmers_ref_short;
	vector<pair<int64_t, int64_t>> v_kmers_data_long, v_kmers_data_short;

	const double ht_long_max_fill_factor = 0.1;

	int ht_short_mask;
	uint32_t ht_long_size;
	uint32_t ht_long_mask;

	int est_len_correction;

	vector<int> ht_long;
	vector<pair<int, int>> ht_short;
	vector<pair<int, int>> ht_short_desc;


	void append(seq_t& seq, uint32_t len, uint8_t x)
	{
		seq.resize(seq.size() + len, x);
//		for (uint32_t i = 0; i < len; ++i)
//			seq.emplace_back(x);
	}

	void append(seq_t& seq, const seq_view& sv, const uint8_t allowed_N, const uint8_t forbidden_N)
	{
		for (uint32_t i = 0; i < sv.size(); ++i)
		{
			auto c = sv[i];
			if (c == forbidden_N)
				seq.emplace_back(allowed_N);
			else
				seq.emplace_back(c);
		}
	}

	void append_rc(seq_t& seq, const seq_view& sv, const uint8_t allowed_N, const uint8_t forbidden_N)
	{
		for (uint32_t i = 0; i < sv.size(); ++i)
		{
			auto c = sv[sv.size() - 1 - i];
			if (c == forbidden_N)
				seq.emplace_back(allowed_N);
			else if (c < code_N)
				seq.emplace_back(code_T - c);
			else
				seq.emplace_back(c);
		}
	}

	uint64_t hash_mm(uint64_t x)
	{
		x ^= x >> 33;
		x *= 0xff51afd7ed558ccdLL;
		x ^= x >> 33;
		x *= 0xc4ceb9fe1a85ec53LL;
		x ^= x >> 33;

		return x;
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

	void prefetch(int pos)
	{
		{
#ifdef _WIN32
			_mm_prefetch((const char*)(seq_ref.data() + pos), _MM_HINT_T0);
#else
			__builtin_prefetch(seq_ref.data() + pos);
#endif
		}
	}

	void prefetch_hts1(int pos)
	{
#ifdef _WIN32
		_mm_prefetch((const char*)(ht_short_desc.data() + pos), _MM_HINT_T0);
#else
		__builtin_prefetch(ht_short_desc.data() + pos);
#endif
	}

	void prefetch_hts2(int pos)
	{
#ifdef _WIN32
		_mm_prefetch((const char*)(ht_short.data() + ht_short_desc[pos].first), _MM_HINT_T0);
#else
		__builtin_prefetch((const char*)(ht_short.data() + ht_short_desc[pos].first), 0);
#endif
	}

	void prefetch_htl(int pos)
	{
#ifdef _WIN32
		_mm_prefetch((const char*)(ht_long.data() + pos), _MM_HINT_T0);
#else
		__builtin_prefetch(ht_long.data() + pos);
#endif
	}

	void prepare_kmers(vector<pair<int64_t, int64_t>>& v_kmers, const seq_t& seq, int len, bool store_all = false);

	void prepare_ht_short();
	void prepare_ht_long();

	int equal_len(int ref_pos, int data_pos, int starting_pos = 0);
	int est_equal_len(int64_t x, int64_t y);

	void compare_ranges(int data_start_pos, int ref_start_pos, int len, bool backward);
	int try_extend_forward(int data_start_pos, int ref_start_pos);
	int try_extend_forward2(int data_start_pos, int ref_start_pos);
	int try_extend_backward(int data_start_pos, int ref_start_pos, int max_len);
	int try_extend_backward2(int data_start_pos, int ref_start_pos, int max_len);

public:
	CParser(const CParams2 &params) :
		params(params)
	{
		est_len_correction = params.min_match_len - (16 - (params.min_distant_match_len - params.min_match_len));
	}

	bool prepare_reference(const seq_view ref_view, uint32_t n_seqs);
	bool prepare_data(const seq_view data_view, uint32_t n_seqs);

	void parse();
	CResults calc_stats();

};
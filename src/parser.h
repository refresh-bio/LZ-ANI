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
	const int HT_EMPTY = -1;
	const int HT_FAIL = -1;
	const double ht_long_max_fill_factor = 0.1;

	CParams params;

	seq_t seq_ref;
	seq_t seq_data;
	uint32_t n_ref_seqs;
	uint32_t n_data_seqs;

	seq_view sv_data;

	vector<factor_t> v_parsing;

	vector<pair<int64_t, int64_t>> v_kmers_ref_long, v_kmers_ref_short;
	vector<pair<int64_t, int64_t>> v_kmers_data_long, v_kmers_data_short;

	int ht_short_mask;
	uint32_t ht_long_size;
	uint32_t ht_long_mask;

	vector<int> ht_long;
	vector<pair<int, int>> ht_short;
	vector<pair<int, int>> ht_short_desc;

	int est_len_correction;

	void append(seq_t& seq, const uint32_t len, const uint8_t x)
	{
		seq.resize(seq.size() + len, x);
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

	uint64_t hash_mm(uint64_t x) const
	{
		x ^= x >> 33;
		x *= 0xff51afd7ed558ccdLL;
		x ^= x >> 33;
		x *= 0xc4ceb9fe1a85ec53LL;
		x ^= x >> 33;

		return x;
	}

	int lzcnt(uint64_t x) const
	{
		x |= x >> 1;
		x |= x >> 2;
		x |= x >> 4;
		x |= x >> 8;
		x |= x >> 16;
		x |= x >> 32;

		return (int)_mm_popcnt_u64(~x);
	}

	int lzcnt32(uint32_t x) const
	{
		x |= x >> 1;
		x |= x >> 2;
		x |= x >> 4;
		x |= x >> 8;
		x |= x >> 16;

		return (int)_mm_popcnt_u32(~x);
	}

	void _my_prefetch(const char* ptr) const
	{
#ifdef _WIN32
		_mm_prefetch(ptr, _MM_HINT_T0);
#else
		__builtin_prefetch(ptr);
#endif
	}

	void prefetch(const int pos) const
	{
		_my_prefetch((const char*)(seq_ref.data() + pos));
	}

	void prefetch_hts1(const int pos) const
	{
		_my_prefetch((const char*)(ht_short_desc.data() + pos));
	}

	void prefetch_hts2(const int pos) const
	{
		_my_prefetch((const char*)(ht_short.data() + ht_short_desc[pos].first));
	}

	void prefetch_htl(const int pos) const
	{
		_my_prefetch((const char*)(ht_long.data() + pos));
	}

	void prepare_kmers(vector<pair<int64_t, int64_t>>& v_kmers, const seq_t& seq, int len, bool store_all = false);

	void prepare_ht_short();
	void prepare_ht_long();

	int equal_len(const int ref_pos, const int data_pos, const int starting_pos = 0) const;
	int est_equal_len(const int64_t x, const int64_t y) const;

	void compare_ranges(const int data_start_pos, const int ref_start_pos, const int len, const bool backward);
//	int try_extend_forward(int data_start_pos, int ref_start_pos);								// OLD
	int try_extend_forward2(const int data_start_pos, const int ref_start_pos);
//	int try_extend_backward(int data_start_pos, int ref_start_pos, int max_len);				// OLD
	int try_extend_backward2(const int data_start_pos, const int ref_start_pos, const int max_len);

public:
	CParser(const CParams &params) :
		params(params)
	{
		est_len_correction = params.min_match_len - (16 - (params.min_distant_match_len - params.min_match_len));
	}

	bool prepare_reference(const seq_view ref_view, uint32_t n_seqs);
	bool prepare_data(const seq_view data_view, uint32_t n_seqs);

	void parse();
	results_t calc_stats();
};

// EOF
#pragma once

#include "defs.h"
#include <string>

using namespace std;

class CWorker
{
	seq_t s_reference;
	seq_t s_data;

	uint32_t n_reference;
	uint32_t n_data;

	uint32_t htl_size;
	uint32_t htl_mask;
	const double htl_max_fill_factor = 0.1;
	vector<int> htl;
	vector<vector<int>> hts;
	vector<CFactor> v_parsing;

	vector<pair<int64_t, int>> v_kmers_l, v_kmers_s;

	void prepare_kmers(vector<pair<int64_t, int>> &v_kmers, const seq_t &seq, int len, bool store_all = false);
	int hash_mm(uint64_t x, int mask);

	void prefetch(int pos);
	void prefetch_hts(int pos);
	void prefetch_htl(int pos);
	int equal_len(int ref_pos, int data_pos, int starting_pos = 0);
	void duplicate_rev_comp(seq_t &seq);

	void compare_ranges(int data_start_pos, int ref_start_pos, int len, bool backward);
	int try_extend_forward(int data_start_pos, int ref_start_pos);
	int try_extend_backward(int data_start_pos, int ref_start_pos, int max_len);

public:
	bool load_data(string fn_ref, string fn_data);
	void swap_data();
	void parse();
	void parse_new();
	void parse_new2();
	void parsing_postprocess();
	void export_parsing();
	void prepare_ht_short();
	void prepare_ht_long();
	void prepare_pf();

	void calc_ani(CResults &res, int mode);
	void calc_ani_thr(CResults &res, int mode);
	bool load_file(const string &file_name, seq_t &seq, uint32_t &n_parts);

	void clear();
};

// EOF

#pragma once

#include "defs.h"
#include <string>

using namespace std;

class CWorker
{
	int codes[256];

	seq_t s_reference;
	seq_t s_data;

	uint32_t n_reference;
	uint32_t n_data;

	uint32_t htl_size;
	uint32_t htl_mask;
	const double htl_max_fill_factor = 0.1;
	int hts_mask;

	vector<int> htl;
	vector<vector<int>> hts;
	vector<vector<pair<int, int>>> hts2;
	vector<pair<int, int>> hts3;
	vector<pair<int, int>> hts3_desc;

	vector<CFactor> v_parsing;

	vector<pair<int64_t, int>> v_kmers_rl, v_kmers_rs;
	vector<pair<int64_t, int>> v_kmers_dl, v_kmers_ds;

	void prepare_kmers(vector<pair<int64_t, int>> &v_kmers, const seq_t &seq, int len, bool store_all = false);
	int hash_mm(uint64_t x, int mask);

	int lzcnt(uint64_t x);
	int lzcnt32(uint32_t x);
	void prefetch(int pos);
	void prefetch_hts1(int pos);
	void prefetch_hts2(int pos);
	void prefetch_htl(int pos);
	int equal_len(int ref_pos, int data_pos, int starting_pos = 0);
	int est_equal_len(int64_t x, int64_t y);

	void duplicate_rev_comp(seq_t &seq);

	void compare_ranges(int data_start_pos, int ref_start_pos, int len, bool backward);
	int try_extend_forward(int data_start_pos, int ref_start_pos);
	int try_extend_forward2(int data_start_pos, int ref_start_pos);
	int try_extend_backward(int data_start_pos, int ref_start_pos, int max_len);
	int try_extend_backward2(int data_start_pos, int ref_start_pos, int max_len);

public:
	CWorker();

	bool load_data(string fn_ref, string fn_data);
	void swap_data();
	void parse();
	void export_parsing();
	void prepare_ht_short();
	void prepare_ht_long();
	void prepare_kmers();

	void calc_ani(CResults &res, int mode);
	bool load_file(const string &file_name, seq_t &seq, uint32_t &n_parts, int separator);

	void clear();
};

// EOF

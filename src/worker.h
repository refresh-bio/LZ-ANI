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

public:
	void prefetch(int pos);
	bool load_data(string fn_ref, string fn_data);
	void parse();
	void parsing_postprocess();
	void export_parsing();
	void prepare_ht_short();
	void prepare_ht_long();
	int my_hash_s(seq_t::iterator p, int len);
	int my_hash_l(seq_t::iterator p, int len);
	int equal_len(int ref_pos, int data_pos, int starting_pos = 0);

	void calc_ani(CResults &res);
	void duplicate_rev_comp(seq_t &seq);
	bool load_file(const string &file_name, seq_t &seq, uint32_t &n_parts);

	void clear();
};

// EOF

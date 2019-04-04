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

	uint32_t ht_size;
	uint32_t ht_mask;
	const double ht_max_fill_factor = 0.5;
	vector<int> ht;
	vector<CFactor> v_parsing;

public:
	bool load_data(string fn_ref, string fn_data);
	void parse();
	void parsing_postprocess();
	void export_parsing();
	void prepare_ht();
	int my_hash(seq_t::iterator p, int len);
	int equal_len(int ref_pos, int data_pos);

	void calc_ani(CResults &res);
	void duplicate_rev_comp(seq_t &seq);
	bool load_file(const string &file_name, seq_t &seq, uint32_t &n_parts);

	void clear();
};

// EOF

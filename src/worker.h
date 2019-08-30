#pragma once

#include "defs.h"
#include "base_worker.h"
#include <string>

using namespace std;

class CWorker : public BaseWorker
{
	seq_t raw_reference; // added to unify s_worker and worker interfaces

	vector<int> htl;
	vector<vector<int>> hts;
	vector<vector<pair<int, int>>> hts2;
	vector<pair<int, int>> hts3;
	vector<pair<int, int>> hts3_desc;

	vector<pair<int64_t, int>> v_kmers_rl, v_kmers_rs;
	vector<pair<int64_t, int>> v_kmers_dl, v_kmers_ds;

	void prefetch(int pos);
	void prefetch_hts1(int pos);
	void prefetch_hts2(int pos);
	void prefetch_htl(int pos);
	int equal_len(int ref_pos, int data_pos, int starting_pos = 0);
	int est_equal_len(int64_t x, int64_t y);

	void compare_ranges(int data_start_pos, int ref_start_pos, int len, bool backward);
	int try_extend_forward(int data_start_pos, int ref_start_pos);
	int try_extend_forward2(int data_start_pos, int ref_start_pos);
	int try_extend_backward(int data_start_pos, int ref_start_pos, int max_len);
	int try_extend_backward2(int data_start_pos, int ref_start_pos, int max_len);

public:
	
	CWorker() : BaseWorker() {
		s_reference = &raw_reference; 
	}

	bool load_data(string fn_ref, string fn_data);
	void swap_data();
	void parse();

	void prepare_ht_short();
	void prepare_ht_long();
	void prepare_kmers();

	void clear();
};

// EOF

#pragma once

#include "worker_base.h"

class CWorker : public CWorkerBase
{
	seq_t s_reference;
	seq_t s_data;

	uint32_t n_reference;
	uint32_t n_data;

	vector<int> htl;
	vector<vector<int>> hts;
	vector<vector<pair<int, int>>> hts2;
	vector<pair<int, int>> hts3;
	vector<pair<int, int>> hts3_desc;


	void prepare_kmers(vector<pair<int64_t, int64_t>> &v_kmers, const seq_t &seq, int len, bool store_all = false);

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
	CWorker(CParams& params) : CWorkerBase(params)
	{
	};

	bool load_data(string fn_ref, string fn_data);
	void swap_data();
	void parse();
	void export_parsing();
	void prepare_ht_short();
	void prepare_ht_long();
	void prepare_kmers();

	void calc_ani(CFatResults&res, int mode);
	bool load_file(const string &file_name, seq_t &seq, uint32_t &n_parts, int separator);

	void clear();
};

// EOF

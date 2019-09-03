#pragma once

#include "defs.h"
#include "base_worker.h"
#include <string>

using namespace std;

class CSharedWorker : public BaseWorker
{
	
	int est_len_correction;

	vector<int> *htl;
	vector<vector<int>> *hts;
	vector<vector<pair<int, int>>> *hts2;
	vector<pair<int, int>> *hts3;
	vector<pair<int, int>> *hts3_desc;

	vector<pair<int64_t, int>> v_kmers_rl, v_kmers_rs;
	vector<pair<int64_t, int>> v_kmers_dl, v_kmers_ds;

	void prefetch_hts1(int pos);
	void prefetch_hts2(int pos);
	void prefetch_htl(int pos);
	int equal_len(int ref_pos, int data_pos, int starting_pos = 0);
	int est_equal_len(int64_t x, int64_t y);

public:
	CSharedWorker();
	~CSharedWorker();

	bool load_reference(string fn_ref, pair<seq_t, int>* buffered_data);
	bool load_data(string fn_data, pair<seq_t, int>* buffered_data);

	bool share_from(CSharedWorker* base);

	void parse();
	void prepare_ht_short();
	void prepare_ht_long();
	void prepare_kmers_data();
	void prepare_kmers_ref_short();
	void prepare_kmers_ref_long();

	void clear_ref();
	void clear_data();
};

// EOF

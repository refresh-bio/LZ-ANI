#pragma once

#include <chrono>

#include "params.h"
#include "defs.h"
#include "seq_reservoir.h"
#include "filter.h"

using namespace std::chrono;

class CLZMatcher
{
	CParams params;

	vector<pair<high_resolution_clock::time_point, string>> times;

	CSeqReservoir seq_reservoir;
	CFilter filter;

	vector<vec_id_results_t> results;

public:
	CLZMatcher(CParams& params) :
		params(params)
	{}

	bool load_sequences();
	bool load_filter();
	bool compare_sequences();
	void reorder_sequences();

	void show_timinigs_info();

	void run_all2all();
	
	bool store_results();
};
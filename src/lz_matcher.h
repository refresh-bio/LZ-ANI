#pragma once

#include <chrono>

#include "params.h"
#include "defs.h"
#include "seq_reservoir.h"
#include "filter.h"

using namespace std::chrono;

class CLZMatcher
{
	CParams2 params;

	vector<pair<high_resolution_clock::time_point, string>> times;

	CSeqReservoir seq_reservoir;
	CFilter filter;

	vector<VecIdResults> results;

public:
	CLZMatcher(CParams2& params) :
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
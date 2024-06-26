// *******************************************************************************************
// This file is a part of LZ-ANI software distributed under GNU GPL 3 license.
// The homepage of the LZ-ANI project is https://github.com/refresh-bio/LZ-ANI
//
// Copyright(C) 2024-2024, S.Deorowicz, A.Gudys
//
// Version: 1.0.0
// Date   : 2024-06-26
// *******************************************************************************************

#pragma once

#include <chrono>

#include "params.h"
#include "defs.h"
#include "seq_reservoir.h"
#include "filter.h"

using namespace std;
using namespace std::chrono;

class CLZMatcher
{
	CParams params;

	vector<pair<high_resolution_clock::time_point, string>> times;

	CSeqReservoir seq_reservoir;
	CFilter filter;

	vector<vec_id_results_t> results;

	bool load_sequences();
	bool load_filter();
	bool compare_sequences();
	void reorder_sequences();

	void show_timinigs_info();

	void do_matching();

	bool store_results();

public:
	CLZMatcher(CParams& params) :
		params(params)
	{}

	bool run_all2all();
};

// EOF
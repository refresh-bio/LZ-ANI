// *******************************************************************************************
// This file is a part of LZ-ANI software distributed under GNU GPL 3 license.
// The homepage of the LZ-ANI project is https://github.com/refresh-bio/LZ-ANI
//
// Copyright(C) 2024-2024, S.Deorowicz, A.Gudys
//
// Version: 1.1.0
// Date   : 2024-09-05
// *******************************************************************************************

#pragma once

#include <chrono>
#include <mutex>

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

	FILE* f_alignment = nullptr;
	char alignment_buffer[1024];
	mutex mtx_alignment;

	bool load_sequences();
	bool load_filter();
	bool compare_sequences();
	void reorder_sequences();

	void show_timinigs_info();

	void do_matching();

	bool store_results();

	void store_alignment(uint64_t idx1, uint64_t idx2, const vector<region_t>& parsing);

public:
	CLZMatcher(CParams& params) :
		params(params),
		seq_reservoir(params.internal_packing)
	{
		if (!params.output_alignment_file_name.empty())
		{
			f_alignment = fopen(params.output_alignment_file_name.c_str(), "wb");

			if (!f_alignment)
			{
				cerr << "Cannot open output file for alignment storage: " << params.output_alignment_file_name << endl;
				exit(1);
			}

			setvbuf(f_alignment, nullptr, _IOFBF, 16 << 20);
			fputs("query\treference\tpident\talnlen\tqstart\tqend\trstart\trend\tnt_match\tnt_mismatch\n", f_alignment);
		}
	}

	~CLZMatcher()
	{
		if (f_alignment)
			fclose(f_alignment);
	}

	bool run_all2all();
};

// EOF
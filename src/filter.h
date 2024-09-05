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

#include <vector>
#include <string>

#include "defs.h"
#include "params.h"
#include "../libs/refresh/compression/lib/file_wrapper.h"

using namespace std;

class CFilter
{
	vector<vector<id_t>> filter;
	vector<string> sequence_names;

	long int local_strtol(const char* str, char** endptr) {
		long int val = 0;
		char* p = (char*)str;
		bool is_negative = false;

		if (*p == '-')
		{
			is_negative = true;
			++p;
		}

		while (*p >= '0' && *p <= '9')
		{
			val = val * 10 + (*p++ - '0');
		}

		if (endptr)
			*endptr = p;

		return is_negative ? -val : val;
	}

public:
	CFilter() = default;

	bool load_filter(const string &fn, double thr, uint32_t no_threads, uint32_t verbosity_level);

	vector<string>& get_sequence_names()
	{
		return sequence_names;
	}

	void clear_sequence_names()
	{
		sequence_names.clear();
		sequence_names.shrink_to_fit();
	}

	bool is_empty()
	{
		return filter.empty();
	}

	void reorder_items(const vector<uint32_t>& reordering_map, uint32_t no_threads, uint32_t verbosity_level);

	vector<id_t> &get_row(size_t i) 
	{ 
		return filter[i];
	}

	void clear_row(size_t i)
	{
		filter[i].clear();
		filter[i].shrink_to_fit();
	}
};

// EOF

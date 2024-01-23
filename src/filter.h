#pragma once

#include <vector>
#include <string>

#include "defs.h"
#include "params.h"
#include "../libs/refresh/file_wrapper.h"

using namespace std;

class CFilter
{
	vector<vector<id_t>> filter;
	vector<string> sequence_names;

public:
	CFilter() = default;

	bool load_filter(const string &fn, double thr, uint32_t no_threads);

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

	void reorder_items(const vector<uint32_t>& reordering_map, uint32_t no_threads);

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
#pragma once

#include <vector>
#include <string>
#include <map>
//#include <unordered_map>

#include <thread>
#include <future>

#include "params.h"
#include "defs.h"


using namespace std;


class CLZMatcher
{
	CParams& params;

	vector<string> input_file_names;

	map<pair<int, int>, CResults> results;
	vector<int> seq_len;

	vector<pair<seq_t, int>> v_buffer_seqs;

	vector<thread> v_threads;
	vector<future<void>> v_fut;

	bool prefetch_input_files();

public:
	CLZMatcher(CParams& params) : params(params)
	{};

	bool run_all2all(vector<string>& _input_file_names);
};
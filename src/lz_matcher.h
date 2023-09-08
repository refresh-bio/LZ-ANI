#pragma once

#include <cinttypes>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <thread>
#include <future>
#include <chrono>

#include "params.h"
#include "defs.h"
#include "s_worker.h"
#include "data_storage.h"

using namespace std;
using namespace std::chrono;

using filter_dict_t = unordered_set<pair_id_t>;
//using filter_dict_t = set<pair_id_t>;

//using results_dict_t = unordered_map<pair_id_t, CResults>;
using results_dict_t = map<pair_id_t, CResults>;

class CLZMatcher
{
	CParams& params;

	CDataStorage data_storage;

	vector<pair<high_resolution_clock::time_point, string>> times;

//	vector<string> input_file_names;
	vector<pair<string, int>> filter_genome_names;
	string filter_name;
	uint32_t filter_thr;
	filter_dict_t filter_set;
	vector<pair_id_t> filter_vec;

	vector<file_desc_t> input_file_desc;

	results_dict_t results;
	vector<int> seq_len;

	vector<thread> v_threads;
	vector<future<void>> v_fut;

	uint64_t encode_pair_id(uint64_t x, uint64_t y)
	{
		return (x << 32) + y;
	}

	uint64_t encode_pair_id_mm(uint64_t x, uint64_t y)
	{
		if (x < y)
			swap(x, y);

		return (x << 32) + y;
	}

	uint64_t swap_pair_id(uint64_t x)
	{
		return (x << 32) + (x >> 32);
	}

	bool is_pair_id_mm(uint64_t x)
	{
		return (x >> 32) <= (x & 0xffffffffull);
	}

	pair<uint32_t, uint32_t> decode_pair_id(uint64_t x)
	{
		return pair<uint32_t, uint32_t>((uint32_t) (x >> 32), (uint32_t) (x & 0xffffffffull));
	}

	bool prefetch_input_files();
	bool load_filter();
	bool reorder_input_files();

	bool prepare_worker_base(CSharedWorker* wb, uint32_t id);

	bool store_results(const string &output_file_name);

	void show_timinigs_info();

public:
	CLZMatcher(CParams& params) : 
		params(params), 
		filter_thr(0),
		data_storage(params.buffer_input_data)
	{};

	bool set_filter(const string& _filter_name, const uint32_t _filter_thr);
	bool init_data_storage(const vector<string>& input_file_names);
	bool init_data_storage(const string & input_file_name);

//	bool run_all2all(vector<string>& _input_file_names, const string& output_file_name);
	bool run_all2all(const string& output_file_name);
};
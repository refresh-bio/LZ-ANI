#pragma once

#include <sstream>
#include <string>
#include <vector>
#include <thread>

using namespace std;

enum class output_type_t {single_file, split_files};

struct CParams
{
public:
	uint32_t verbosity_level = 1;

	uint32_t no_threads = 0;

	int min_match_len = 8;
	int min_close_match_len = 8;
	int min_distant_match_len = 20;
	int close_dist = 256;
	int max_lit_run_in_match = 32;
	double min_coverage = 0.10;
	int min_region_len = 128;
	int approx_window = 16;
	int approx_mismatches = 5;
	int approx_run_len = 3;
	bool multisample_fasta = false;
	double filter_thr = 0.0;

	output_type_t output_type = output_type_t::single_file;

	vector<string> input_file_names;
	string output_file_name;
	string filter_file_name;

	bool store_total_ani = false;
	bool store_global_ani = false;
	bool store_local_ani = false;
	bool store_coverage = false;
	bool store_regions = false;
	bool store_full_seq_ids = false;
	bool store_condensed = false;

	string str()
	{
		stringstream ss;

		ss << "[params]" << endl;
		ss << "min_match_len         : " << min_match_len << endl;
		ss << "min_close_match_len   : " << min_close_match_len << endl;
		ss << "min_distant_match_len : " << min_distant_match_len << endl;
		ss << "close_dist            : " << close_dist << endl;
		ss << "max_lit_run_in_match  : " << max_lit_run_in_match << endl;
		ss << "min_coverage          : " << min_coverage << endl;
		ss << "min_region_len        : " << min_region_len << endl;
		ss << "approx_window         : " << approx_window << endl;
		ss << "approx_mismatches     : " << approx_mismatches << endl;
		ss << "approx_run_len        : " << approx_run_len << endl;
		ss << "multisample_fasta     : " << boolalpha << multisample_fasta << noboolalpha << endl;
		ss << "filter_thr            : " << filter_thr << endl;
		ss << "store_total_ani       : " << boolalpha << store_total_ani << endl;
		ss << "store_global_ani      : " << boolalpha << store_global_ani << endl;
		ss << "store_local_ani       : " << boolalpha << store_local_ani << endl;
		ss << "store_coverage        : " << boolalpha << store_coverage << endl;
		ss << "store_regions         : " << boolalpha << store_regions << endl;
		ss << "store_full_seq_ids    : " << boolalpha << store_full_seq_ids << endl;
		ss << "store_condensed       : " << boolalpha << store_condensed << endl;
		ss << noboolalpha;

		ss << "no_threads            : " << no_threads << endl;

		ss << "output_file_name      : " << output_file_name << endl;
		ss << "filter_file_name      : " << filter_file_name << endl;
		ss << "input_file_names      : ";
		for (size_t i = 0; i + 1 < input_file_names.size(); ++i)
			ss << input_file_names[i] << ", ";
		ss << input_file_names.back() << endl;

		return ss.str();
	}

	void adjust_threads()
	{
		if (no_threads == 0)
		{
			no_threads = thread::hardware_concurrency();
			if (!no_threads)
				no_threads = 1;
		}
	}
};

// EOF
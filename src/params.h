#pragma once

#include <sstream>
#include <string>
#include <vector>
#include <thread>

using namespace std;

enum class output_type_t {single_file, split_files};

struct CParams2
{
public:

	int verbosity_level = 1;

	int no_threads = 0;

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
	bool store_ani = false;
	bool store_cov = false;
	bool store_shorter_ani = false;
	bool store_shorter_cov = false;
	bool store_full_seq_ids = false;


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
		ss << "multisample_fasta     : " << multisample_fasta << endl;
		ss << "no_threads            : " << no_threads << endl;

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

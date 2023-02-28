#pragma once

#include <sstream>
#include <string>

using namespace std;

enum class output_mode_t {total_ani, genomic_ani_cov, sym_in_matches_literals};

struct CParams
{
private:
	string output_mode_str()
	{
		if (output_mode == output_mode_t::total_ani)
			return "total_ani";
		if (output_mode == output_mode_t::genomic_ani_cov)
			return "genomic_ani_cov";
		if (output_mode == output_mode_t::sym_in_matches_literals)
			return "sym_in_matches_literals";

		return "unknown";
	}

public:

	int verbosity_level = 1;

	int no_threads = 0;

	output_mode_t output_mode = output_mode_t::sym_in_matches_literals;

	int min_match_len = 8;
	int min_close_match_len = 8;
	int min_distant_match_len = 20;
	int close_dist = 256;
	//int LONG_LITERAL_RUN_LEN = DEF_LONG_LITERAL_RUN_LEN;
	int max_lit_tun_in_match = 32;
	double min_coverage = 0.10;
	int min_region_len = 128;
	int approx_window = 16;
	int approx_mismatches = 5;
	int approax_run_len = 3;

	bool buffer_input_data = true;	// In the current implementation must be true!
	bool output_dense_matrix = false;

	string str()
	{
		stringstream ss;

		ss << "[params]" << endl;
		ss << "min_match_len         : " << min_match_len << endl;
		ss << "min_close_match_len   : " << min_close_match_len << endl;
		ss << "min_distant_match_len : " << min_distant_match_len << endl;
		ss << "close_dist            : " << close_dist << endl;
		ss << "max_lit_tun_in_match  : " << max_lit_tun_in_match << endl;
		ss << "min_coverage          : " << min_coverage << endl;
		ss << "min_region_len        : " << min_region_len << endl;
		ss << "approx_window         : " << approx_window << endl;
		ss << "approx_mismatches     : " << approx_mismatches << endl;
		ss << "approx_run_len        : " << approax_run_len << endl;
		ss << "output_mode           : " << output_mode_str() << endl;
		ss << "no_threads            : " << no_threads << endl;

		return ss.str();
	}
};

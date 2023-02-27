#pragma once

struct CParams
{
	int verbosity_level = 1;

	int no_threads = 0;

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


};

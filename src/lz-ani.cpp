// ani-entropy.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

// all2all --in-file-names fl/ds50_0.txt -out aaa_new -t 32 -bs -cd 128 -reg 48 -mml 6 -mdl 11 -aw 16 -am 6 -ar 2  -filter kmer-db.db/ds50_0.a2a 3 --verbose 2

#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>

#include "defs.h"
#include "lz_matcher.h"

using namespace refresh;
using namespace std;

enum class working_mode_t {none, all2all};

working_mode_t working_mode;

CParams params;

bool parse_params(int argc, char** argv);
vector<string> load_input_names(const string& fn);
void usage();

// ****************************************************************************
void usage()
{
	cerr << "lz-ani <mode> [options]\n";
	cerr << "Modes:\n";
	cerr << "   all2all                - all to all\n";
	cerr << "Options:\n";
	cerr << "   --in-fasta <file_name> - FASTA file (for multisample-fasta mode)\n";
	cerr << "   --in-txt <file_name>   - text file with FASTA file names\n";

	cerr << "   -out <file_name>       - output file name\n";
	cerr << "   -t <val>               - no of threads (default: " << params.no_threads << ")\n";
	cerr << "   -mml <val>             - min. match length (default: " << params.min_match_len << ")\n";
	cerr << "   -mdl <val>             - min. distant length (default: " << params.min_distant_match_len << ")\n";
	cerr << "   -cd <val>              - max. dist. between close matches (default: " << params.close_dist << ")\n";
	cerr << "   -mlrin <val>           - max. literal run len. in match (default: " << params.max_lit_run_in_match << ")\n";
	cerr << "   -cov <val>             - min. coverage threshold (default: " << params.min_coverage << ")\n";
	cerr << "   -reg <val>             - min. considered region length (default: " << params.min_region_len << ")\n";
	cerr <<	"   -aw <val>              - approx. window length (default: " << params.approx_window << ")\n";
	cerr << "   -am <val>              - max. no. of mismatches in approx. window (default: " << params.approx_mismatches << ")\n";
	cerr << "   -ar <val>              - min. length of run ending approx. extension (default: " << params.approx_run_len << ")\n";
	cerr << "   -filter <file_name> <min_val> - filtering file (kmer-db output) and threshold\n";
	cerr << "   --verbose <int>        - verbosity level (default: " << params.verbosity_level << ")\n";
	cerr << "   --multisample-fasta    - multi sample FASTA input (default: " << params.multisample_fasta << ")\n";

	cerr << "   --output-type <type>   - one of: 'single-file', 'split-files' (default: " << "single-file" << ")\n";
	cerr << "   --store-total-ani      - store total ANI\n";
	cerr << "   --store-ani            - store ANI in both directions\n";
	cerr << "   --store-cov            - store coverage in both directions\n";
	cerr << "   --store-regions        - store no. of regions in both directions\n";
	cerr << "   --store-shorter-ani    - store ANI for shorter sequence\n";
	cerr << "   --store-shorter-cov    - store coverage for shorter sequence\n";
	cerr << "   --store-full-seq-ids   - store full sequence ids in main file\n";
}

// ****************************************************************************
vector<string> load_input_names(const string& fn)
{
	ifstream ifs(fn);
	vector<string> vec;

	if (!ifs.is_open())
	{
		cerr << "Cannot open file: " << fn << endl;
		return vec;
	}

	vec.assign(istream_iterator<string>(ifs), istream_iterator<string>());

	return vec;		
}

// ****************************************************************************
bool parse_params(int argc, char** argv)
{
	if (argc < 3)
	{
		usage();
		return false;
	}

	working_mode = working_mode_t::none;

	if (argv[1] == "all2all"s)
		working_mode = working_mode_t::all2all;
	else
	{
		cerr << "Unknown mode: " << argv[1] << endl;
		usage();
		return false;
	}

	for (int i = 2; i < argc;)
	{
		string par = string(argv[i]);

		if (par == "--in-txt"s && i + 1 < argc)
		{
			params.input_file_names = load_input_names(argv[i+1]);
			if (params.input_file_names.empty())
				return false;

			i += 2;
		}
		else if (par == "--in-fasta"s && i + 1 < argc)
		{
			params.input_file_names.clear();
			params.input_file_names.emplace_back(argv[i + 1]);
			i += 2;
		}
		else if (par == "-out"s)
		{
			if (i + 1 >= argc)
			{
				cerr << "Unknown out name\n";
				return false;
			}

			params.output_file_name = argv[i + 1];
			i += 2;
		}
		else if (par == "-t")
		{
			params.no_threads = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-mml")
		{
			params.min_match_len = atoi(argv[i + 1]);
			params.min_close_match_len = params.min_match_len;
			i += 2;
		}
		else if (par == "-mdl")
		{
			params.min_distant_match_len = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-cd")
		{
			params.close_dist = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-mlrim")
		{
			params.max_lit_run_in_match = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-cov")
		{
			params.min_coverage = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-reg")
		{
			params.min_region_len = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-aw")
		{
			params.approx_window = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-am")
		{
			params.approx_mismatches = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-ar")
		{
			params.approx_run_len = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "-filter" && i + 2 < argc)
		{
			params.filter_file_name = argv[i + 1];
			params.filter_thr = atof(argv[i + 2]);
			i += 3;
		}
		else if (par == "--verbose"s && i + 1 < argc)
		{
			params.verbosity_level = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "--output-type"s && i + 1 < argc)
		{
			string par_type = argv[i + 1];
			if (par_type == "single-file"s)
				params.output_type = output_type_t::single_file;
			else if (par_type == "split-files"s)
				params.output_type = output_type_t::split_files;
			else
			{
				cerr << "Unknown output-type: " << par_type << endl;
				usage();
				exit(0);
			}
			i += 2;
		}
		else if (par == "--multisample-fasta"s)
		{
			params.multisample_fasta = true;
			++i;
		}
		else if (par == "--store-total-ani"s)
		{
			params.store_total_ani = true;
			++i;
		}
		else if (par == "--store-global-ani"s)
		{
			params.store_global_ani = true;
			++i;
		}
		else if (par == "--store-local-ani"s)
		{
			params.store_local_ani = true;
			++i;
		}
		else if (par == "--store-coverage"s)
		{
			params.store_coverage = true;
			++i;
		}
		else if (par == "--store-regions"s)
		{
			params.store_regions = true;
			++i;
		}
		else if (par == "--store-full-seq-ids"s)
		{
			params.store_full_seq_ids = true;
			++i;
		}
		else if (par == "--store-condensed"s)
		{
			params.store_condensed = true;
			++i;
		}
		else
		{
			cerr << "Unknown parameter: " << string(argv[i]) << endl;
			usage();
			exit(0);
		}
	}

	if (working_mode == working_mode_t::all2all && params.input_file_names.empty())
	{
		cerr << "Input file names not provided\n";
		return false;
	}

	return true;
}

// ****************************************************************************
void split(const std::string& str, std::vector<std::string>& parts, char sep)
{
	parts.clear();

	std::string s;

	for (auto c : str)
	{
		if (c == sep)
		{
			parts.emplace_back(s);
			s.clear();
		}
		else
			s.push_back(c);
	}

	if (!s.empty())
		parts.emplace_back(s);
}

// ****************************************************************************
int main(int argc, char **argv)
{
	if (!parse_params(argc, argv))
		return 0;

	params.adjust_threads();
	params.multisample_fasta = true;

	CLZMatcher lzm(params);

	if (!lzm.load_sequences())
		return 0;
	
	if (!lzm.load_filter())
		return 0;

	if (!lzm.compare_sequences())
		return 0;

	lzm.reorder_sequences();

	lzm.run_all2all();

	return 0;
}

// EOF

// ani-entropy.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

// all2all --in-file-names fl/ds50_0.txt -out aaa_new -t 32 -bs -cd 128 -reg 48 -mml 6 -mdl 11 -aw 16 -am 6 -ar 2  -filter kmer-db.db/ds50_0.a2a 3 --verbose 2

//all2all -t 1 --verbose 2 -cd 128 -reg 80 -mml 7 -mdl 11 -mlrim 32 -aw 16 -am 6 -ar 2 -out err.lzani  --in-fasta last_err.fna -filter err.filter 0.9 --output-type split-files --store-total-ani --store-local-ani --store-global-ani --store-coverage --store-condensed
// ../../lz-ani-0.2 all2all --multisample-fasta --in-fasta pair.fna -out aa -t 16 -filter pair.a2a.ani-shorter 0.7

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
	cerr << LZ_ANI_INFO << endl;
	cerr << "Tool for rapid determination of similarities among sets of DNA sequences\n";
	cerr << "Usage:\n";
	cerr << "lz-ani <mode> [options]\n";
	cerr << "Modes:\n";
	cerr << "   all2all                - all to all\n";
	cerr << "Options:\n";
	cerr << "   --in-fasta <file_name> - FASTA file (for multisample-fasta mode)\n";
	cerr << "   --in-txt <file_name>   - text file with FASTA file names\n";
	cerr << "   --in-dir <path>        - directory with FASTA files\n";

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
	cerr << "   -filter-kmerdb <file_name> <min_val> - filtering file (kmer-db output) and threshold\n";
//	cerr << "   -filter-pairs <file_name> - filtering file (tsv with pairs)\n";
	cerr << "   --verbose <int>        - verbosity level (default: " << params.verbosity_level << ")\n";
	cerr << "   --multisample-fasta <true|false> - multi sample FASTA input (default: " << boolalpha << params.multisample_fasta << noboolalpha << ")\n";

	cerr << "   --output-type <type>   - one of: 'single-file', 'split-files' (default: " << "single-file" << ")\n";
	cerr << "   --output-format <type> - comma-separated list of values: " << CParams::list_component_types() << " (default: " << params.output_format << "), you can include also meta-names:" << endl;
	
	for (const auto& x : CParams::list_component_metas())
		cerr << "                          - " << x << endl;
//	cerr << "   --store-condensed      - \n";												
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
		if (par == "--in-dir"s && i + 1 < argc)
		{
			filesystem::path fp(argv[i+1]);

			try
			{
				params.input_file_names.clear();

				filesystem::directory_iterator fsdi(fp);

				for (const auto& fs : fsdi)
					params.input_file_names.push_back(fs.path().string());

			}
			catch (...)
			{
				cerr << "Non-existing directory: " << argv[i+1] << endl;
				return false;
			}

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
		else if (par == "-filter-kmerdb" && i + 2 < argc)
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
		else if (par == "--output-format"s && i + 1 < argc)
		{
			auto ret = params.parse_output_format(argv[i + 1]);
			if (!ret.empty())
			{
				cerr << "Unknown output-format component: " << ret;
				return false;
			}
			i += 2;
		}
		else if (par == "--multisample-fasta"s && i + 1 < argc)
		{
			if(argv[i + 1] == "true"s)
				params.multisample_fasta = true;
			else if(argv[i + 1] == "false"s)
				params.multisample_fasta = false;
			else
			{
				cerr << "Unknown value for --multisample-fasta: " << argv[1] << endl;
				return false;
			}
			i += 2;
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
/*void split(const std::string& str, std::vector<std::string>& parts, char sep)
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
}*/

// ****************************************************************************
int main(int argc, char **argv)
{
	if (!parse_params(argc, argv))
		return 0;

	params.adjust_threads();
	params.multisample_fasta = true;

	CLZMatcher lzm(params);

	lzm.run_all2all();

	return 0;
}

// EOF

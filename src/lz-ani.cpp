// *******************************************************************************************
// This file is a part of LZ-ANI software distributed under GNU GPL 3 license.
// The homepage of the LZ-ANI project is https://github.com/refresh-bio/LZ-ANI
//
// Copyright(C) 2024-2024, S.Deorowicz, A.Gudys
//
// Version: 1.1.0
// Date   : 2024-09-05
// *******************************************************************************************

#include <iostream>
#include <fstream>
#include <iterator>
#include <cstdio>
#include <algorithm>

#include "defs.h"
#include "lz_matcher.h"

#ifdef _MSC_VER 
#include <mimalloc.h>
#include <mimalloc-new-delete.h>
#endif

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
	cerr << "  all2all                        - all to all\n";
	
	cerr << "Options - input specification:\n";
	cerr << "      --in-fasta <file_name>     - FASTA file (for multisample-fasta mode)\n";
	cerr << "      --in-txt <file_name>       - text file with FASTA file names\n";
	cerr << "      --in-dir <path>            - directory with FASTA files\n";
	cerr << "      --multisample-fasta <bool> - multi sample FASTA input (default: " << boolalpha << params.multisample_fasta << noboolalpha << ")\n";
	cerr << "      --flt-kmerdb <fn> <float>  - filtering file (kmer-db output) and threshold\n";
	//	cerr << "      --flt-pairs <file_name>    - filtering file (tsv with pairs)\n";

	cerr << "Options - output specification:\n";
	cerr << "  -o, --out <file_name>          - output file name\n";
	cerr << "      --out-ids <file_name>      - output file name for ids file (optional)\n";
	cerr << "      --out-alignment <file_name>- output file name for ids file (optional)\n";
	cerr << "      --out-in-percent <bool>    - output in percent (default: " << boolalpha << params.output_in_percent << noboolalpha << ")\n";
	cerr << "      --out-type <type>          - one of:\n";
	cerr << "                                   tsv - two tsv files with: results defined by --out-format and sequence ids (default)\n";
	cerr << "                                   single-txt - combined results in single txt file\n";
	cerr << "      --out-format <type>        - comma-separated list of values: " << endl;
	cerr << "                                   " << CParams::list_component_types() << endl;
	cerr << "                                   you can include also meta-names:" << endl;
	for (const auto& x : CParams::list_component_metas())
		cerr << "                                   " << x << endl;
	cerr << "                                   (default: " << params.output_format << ")" << endl;
	cerr << "      --out-filter <par> <float> - store only results with <par> (can be: tani, gani, ani, cov) at least <float>; can be used multiple times" << endl;

	cerr << "Options - LZ-parsing-related:\n";
	cerr << "  -a, --mal <int>                - min. anchor length (default: " << params.min_anchor_len << ")\n";
	cerr << "  -s, --msl <int>                - min. seed length (default: " << params.min_seed_len << ")\n";
	cerr << "  -r, --mrd <int>                - max. dist. between approx. matches in reference (default: " << params.max_dist_in_ref << ")\n";
	cerr << "  -q, --mqd <int>                - max. dist. between approx. matches in query (default: " << params.max_dist_in_query << ")\n";
	cerr << "  -g, --reg <int>                - min. considered region length (default: " << params.min_region_len << ")\n";
	cerr <<	"      --aw <int>                 - approx. window length (default: " << params.approx_window << ")\n";
	cerr << "      --am <int>                 - max. no. of mismatches in approx. window (default: " << params.approx_mismatches << ")\n";
	cerr << "      --ar <int>                 - min. length of run ending approx. extension (default: " << params.approx_run_len << ")\n";

	cerr << "Options - other:\n";
	cerr << "  -t, --threads <int>            - no of threads; 0 means auto-detect (default: " << params.no_threads << ")\n";
	cerr << "  -V, --verbose <int>            - verbosity level (default: " << params.verbosity_level << ")\n";
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
	if (argc == 2 && argv[1] == "--version"s)
	{
		cerr << LZ_ANI_VERSION << endl;
		return true;
	}

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
		else if (par == "--in-dir"s && i + 1 < argc)
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
		else if ((par == "-o"s || par == "--out"s) && i + 1 < argc)
		{
			if (i + 1 >= argc)
			{
				cerr << "Unknown out name\n";
				return false;
			}

			params.output_file_name = argv[i + 1];
			i += 2;
		}
		else if (par == "--out-ids"s && i + 1 < argc)
		{
			if (i + 1 >= argc)
			{
				cerr << "Unknown out name\n";
				return false;
			}

			params.output_ids_file_name = argv[i + 1];
			i += 2;
		}
		else if (par == "--out-alignment"s && i + 1 < argc)
		{
			if (i + 1 >= argc)
			{
				cerr << "Unknown out name\n";
				return false;
			}

			params.output_alignment_file_name = argv[i + 1];
			i += 2;
		}
		else if ((par == "-t"s || par == "--threads"s) && i + 1 < argc)
		{
			params.no_threads = atoi(argv[i + 1]);
			i += 2;
		}
		else if ((par == "-s"s || par == "--msl"s) && i + 1 < argc)
		{
			params.min_seed_len = atoi(argv[i + 1]);
			i += 2;
		}
		else if ((par == "-a" || par == "--mal"s) && i + 1 < argc)
		{
			params.min_anchor_len = atoi(argv[i + 1]);
			i += 2;
		}
		else if ((par == "-r"s || par == "--mrd"s) && i + 1 < argc)
		{
			params.max_dist_in_ref = atoi(argv[i + 1]);
			i += 2;
		}
		else if ((par == "-q" || par == "--mqd"s) && i + 1 < argc)
		{
			params.max_dist_in_query = atoi(argv[i + 1]);
			i += 2;
		}
		else if ((par == "-g"s || par == "--reg"s) && i + 1 < argc)
		{
			params.min_region_len = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "--aw"s && i + 1 < argc)
		{
			params.approx_window = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "--am"s && i + 1 < argc)
		{
			params.approx_mismatches = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "--ar"s && i + 1 < argc)
		{
			params.approx_run_len = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "--flt-kmerdb"s && i + 2 < argc)
		{
			params.filter_file_name = argv[i + 1];
			params.filter_thr = atof(argv[i + 2]);
			i += 3;
		}
		else if ((par == "-V"s || par == "--verbose"s) && i + 1 < argc)
		{
			params.verbosity_level = atoi(argv[i + 1]);
			i += 2;
		}
		else if (par == "--out-type"s && i + 1 < argc)
		{
			string par_type = argv[i + 1];
			if (par_type == "single-txt"s)
				params.output_type = output_type_t::single_txt;
			else if (par_type == "tsv"s)
				params.output_type = output_type_t::two_tsv;
			else
			{
				cerr << "Unknown output-type: " << par_type << endl;
				usage();
				exit(0);
			}
			i += 2;
		}
		else if (par == "--out-format"s && i + 1 < argc)
		{
			auto ret = params.parse_output_format(argv[i + 1]);
			if (!ret.empty())
			{
				cerr << "Unknown output-format component: " << ret;
				return false;
			}
			i += 2;
		}
		else if (par == "--out-filter"s && i + 2 < argc)
		{
			if (!params.set_output_filter(argv[i + 1], argv[i + 2]))
			{
				cerr << "Unknown output-filter component: " << argv[i + 1] << " " << argv[i + 2] << endl;
				return false;
			}
			i += 3;
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
		else if (par == "--out-in-percent"s && i + 1 < argc)
		{
			if(argv[i + 1] == "true"s)
				params.output_in_percent = true;
			else if(argv[i + 1] == "false"s)
				params.output_in_percent = false;
			else
			{
				cerr << "Unknown value for --out-in-percent: " << argv[1] << endl;
				return false;
			}
			i += 2;
		}
		else
		{
			cerr << "Unknown parameter: " << string(argv[i]) << endl;
			usage();
			exit(1);
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
int main(int argc, char **argv)
{
	if (!parse_params(argc, argv))
		return 0;

	if (working_mode == working_mode_t::none)
		return 0;

	params.adjust_threads();

	CLZMatcher lzm(params);

	if (!lzm.run_all2all())
		exit(1);

	return 0;
}

// EOF

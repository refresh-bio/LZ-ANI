// *******************************************************************************************
// This file is a part of LZ-ANI software distributed under GNU GPL 3 license.
// The homepage of the LZ-ANI project is https://github.com/refresh-bio/LZ-ANI
//
// Copyright(C) 2024-2024, S.Deorowicz, A.Gudys
//
// Version: 1.0.0
// Date   : 2024-06-26
// *******************************************************************************************

#pragma once

#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <map>
#include <array>

#include "utils.h"

using namespace std;

enum class output_type_t {single_txt, two_tsv};
enum class output_component_t {seq_id1,seq_id2,seq_idx1,seq_idx2,len1,len2,total_ani,global_ani,local_ani,cov,len_ratio,nt_match,nt_mismatch,no_reg};

struct CParams
{
public:
	uint32_t verbosity_level = 1;

	uint32_t no_threads = 0;

	int min_anchor_len = 11;
	int min_seed_len = 7;
	int max_dist_in_ref = 40;
	int max_dist_in_query = 40;
	int min_region_len = 35;
	int approx_window = 15;
	int approx_mismatches = 7;
	int approx_run_len = 3;
	bool multisample_fasta = true;
	double filter_thr = 0.0;
	bool output_in_percent = false;

	output_type_t output_type = output_type_t::two_tsv;

	vector<string> input_file_names;
	string output_file_name;
	string output_ids_file_name;
	string filter_file_name;

	vector<output_component_t> output_components;
	string output_format = "standard";

	static inline std::map<std::string, std::string> std_comp = {
	{"complete", "idx1,idx2,id1,id2,tani,gani,ani,cov,num_alns,len_ratio,len1,len2,nt_match,nt_mismatch"},
	{"standard", "idx1,idx2,id1,id2,tani,gani,ani,cov,num_alns,len_ratio"},
	{"lite", "idx1,idx2,tani,gani,ani,cov,num_alns,len_ratio"}
	};

	const map<string, output_component_t> comp_name_id = {
		{"id1", output_component_t::seq_id1},
		{"id2", output_component_t::seq_id2},
		{"idx1", output_component_t::seq_idx1},
		{"idx2", output_component_t::seq_idx2},
		{"len1", output_component_t::len1},
		{"len2", output_component_t::len2},
		{"tani", output_component_t::total_ani},
		{"gani", output_component_t::global_ani},
		{"ani", output_component_t::local_ani},
		{"cov", output_component_t::cov},
		{"len_ratio", output_component_t::len_ratio},
		{"nt_match", output_component_t::nt_match},
		{"nt_mismatch", output_component_t::nt_mismatch},
		{"num_alns", output_component_t::no_reg}
	};

	const map<string, output_component_t> comp_flt_id = {
		{"tani", output_component_t::total_ani},
		{"gani", output_component_t::global_ani},
		{"ani", output_component_t::local_ani},
		{"cov", output_component_t::cov}
	};

	map<output_component_t, string> comp_id_name;
	uint64_t output_filter_mask = 0;
	vector<double> output_filter_vals;

	CParams()
	{
		for (const auto& x : comp_name_id)
			comp_id_name.emplace(x.second, x.first);

		output_filter_vals.resize(comp_name_id.size(), 0);
		output_filter_mask = 0;

		parse_output_format(output_format);
	}

	static string list_component_types()
	{
		return "id1,id2,idx1,idx2,len1,len2,tani,gani,ani,cov,len_ratio,nt_match,nt_mismatch,num_alns";
	}

	static vector<string> list_component_metas()
	{
		vector<string> ret;
		for (const auto& x : std_comp)
			ret.emplace_back(x.first + "=" + x.second);

		return ret;
	}

	string str()
	{
		stringstream ss;

		ss << "[params]" << endl;
		ss << "min_anchor_len        : " << min_anchor_len << endl;
		ss << "min_seed_len          : " << min_seed_len << endl;
		ss << "max_dist_in_ref       : " << max_dist_in_ref << endl;
		ss << "max_dist_in_query     : " << max_dist_in_query << endl;
		ss << "min_region_len        : " << min_region_len << endl;
		ss << "approx_window         : " << approx_window << endl;
		ss << "approx_mismatches     : " << approx_mismatches << endl;
		ss << "approx_run_len        : " << approx_run_len << endl;
		ss << "multisample_fasta     : " << boolalpha << multisample_fasta << noboolalpha << endl;
		ss << "filter_thr            : " << filter_thr << endl;
		ss << "output_format         : " << output_format << endl;
		ss << "output_in_percent     : " << boolalpha << output_in_percent << noboolalpha << endl;
		ss << noboolalpha;

		ss << "no_threads            : " << no_threads << endl;

		ss << "output_file_name      : " << output_file_name << endl;
		ss << "output_ids_file_name  : " << output_ids_file_name << endl;
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

	string parse_output_format(const string &of)
	{
		output_components.clear();

		vector<string> tmp = split(of, ',');
		vector<string> comp_str;

		// Replace meta-names
		for (const auto& x : tmp)
		{
			auto p = std_comp.find(x);
			if (p == std_comp.end())
				comp_str.emplace_back(x);
			else
			{
				auto tmp2 = split(p->second, ',');
				comp_str.insert(comp_str.end(), tmp2.begin(), tmp2.end());
			}
		}

		for (const auto& x : comp_str)
		{
			auto p = comp_name_id.find(x);
			if (p == comp_name_id.end())
				return x;							// Error: unknown component type
			output_components.emplace_back(p->second);
		}

		return "";									// Ok
	}

	bool set_output_filter(const string& par, const string& val)
	{
		auto p = comp_flt_id.find(par);

		if (p == comp_flt_id.end())
			return false;

		output_filter_mask |= 1ull << (uint32_t)(p->second);
		output_filter_vals[(uint32_t)(p->second)] = stod(val);

		return true;
	}
};

// EOF
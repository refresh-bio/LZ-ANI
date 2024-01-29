#pragma once

#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <map>
#include <array>

#include "utils.h"

using namespace std;

enum class output_type_t {single_file, split_files};
enum class output_component_t {seq_id1,seq_id2,seq_idx1,seq_idx2,len1,len2,total_ani,global_ani,local_ani,cov,sim,len_ratio,sym_mat,sym_lit,no_reg};

struct CParams
{
public:
	uint32_t verbosity_level = 1;

	uint32_t no_threads = 0;

	int min_match_len = 8;
	int min_close_match_len = 7;
	int min_anchor_len = 11;
	int close_dist = 128;
	int max_lit_run_in_match = 32;
	double min_coverage = 0.10;
	int min_region_len = 80;
	int approx_window = 16;
	int approx_mismatches = 6;
	int approx_run_len = 2;
	bool multisample_fasta = true;
	double filter_thr = 0.0;

	output_type_t output_type = output_type_t::split_files;

	vector<string> input_file_names;
	string output_file_name;
	string filter_file_name;

	vector<output_component_t> output_components;
	string output_format = "FULL";

//	static const std::map<std::string, std::string> std_comp;

	static inline std::map<std::string, std::string> std_comp = {
	{"ALL", "seq_id1,seq_id2,seq_idx1,seq_idx2,len1,len2,total_ani,global_ani,local_ani,cov,sim,len_ratio,sym_mat,sym_lit,no_reg"},
	{"FULL", "seq_id1,seq_id2,total_ani,global_ani,local_ani,cov,len_ratio"},
	{"META", "seq_idx1,seq_idx2,sim,local_ani,cov"},
	{"TAXA", "seq_idx1,seq_idx2,total_ani"}
	};

	const map<string, output_component_t> comp_name_id = {
		{"seq_id1", output_component_t::seq_id1},
		{"seq_id2", output_component_t::seq_id2},
		{"seq_idx1", output_component_t::seq_idx1},
		{"seq_idx2", output_component_t::seq_idx2},
		{"len1", output_component_t::len1},
		{"len2", output_component_t::len2},
		{"total_ani", output_component_t::total_ani},
		{"global_ani", output_component_t::global_ani},
		{"local_ani", output_component_t::local_ani},
		{"cov", output_component_t::cov},
		{"sim", output_component_t::sim},
		{"len_ratio", output_component_t::len_ratio},
		{"sym_mat", output_component_t::sym_mat},
		{"sym_lit", output_component_t::sym_lit},
		{"no_reg", output_component_t::no_reg}
	};

	const map<string, output_component_t> comp_flt_id = {
		{"total_ani", output_component_t::total_ani},
		{"global_ani", output_component_t::global_ani},
		{"local_ani", output_component_t::local_ani},
		{"cov", output_component_t::cov},
		{"sim", output_component_t::sim}
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
		return "seq_id1,seq_id2,seq_idx1,seq_idx2,len1,len2,total_ani,global_ani,local_ani,cov,sim,len_ratio,sym_mat,sym_lit,no_reg";
	}

	static vector<string> list_component_metas()
	{
		vector<string> ret;
		for (const auto& x : std_comp)
			ret.emplace_back(x.first + "=" + x.second);

		return ret;
	}

	bool store_condensed = false;

	string str()
	{
		stringstream ss;

		ss << "[params]" << endl;
		ss << "min_match_len         : " << min_match_len << endl;
		ss << "min_close_match_len   : " << min_close_match_len << endl;
		ss << "min_anchor_len        : " << min_anchor_len << endl;
		ss << "close_dist            : " << close_dist << endl;
		ss << "max_lit_run_in_match  : " << max_lit_run_in_match << endl;
		ss << "min_coverage          : " << min_coverage << endl;
		ss << "min_region_len        : " << min_region_len << endl;
		ss << "approx_window         : " << approx_window << endl;
		ss << "approx_mismatches     : " << approx_mismatches << endl;
		ss << "approx_run_len        : " << approx_run_len << endl;
		ss << "multisample_fasta     : " << boolalpha << multisample_fasta << noboolalpha << endl;
		ss << "filter_thr            : " << filter_thr << endl;
		ss << "output_format         : " << output_format << endl;
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
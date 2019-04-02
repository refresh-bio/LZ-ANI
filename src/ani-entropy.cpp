// ani-entropy.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"

#include "defs.h"
#include "files.h"

#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std::chrono;

seq_t s_reference;
seq_t s_data;

uint32_t n_reference;
uint32_t n_data;

uint32_t ht_size;
uint32_t ht_mask;
double ht_max_fill_factor = 0.5;
vector<int> ht;
vector<CFactor> v_parsing;

bool load_data(int argc, char **argv);
void parse();
void parsing_postprocess();
void export_parsing();
void prepare_ht();
int my_hash(seq_t::iterator p, int len);
int equal_len(int ref_pos, int data_pos);

void calc_ani();

// ****************************************************************************
int equal_len(int ref_pos, int data_pos)
{
	int r;
	int ref_len = (int) s_reference.size();
	int data_len = (int) s_data.size();

	for (r = 0; ref_pos + r < ref_len && data_pos + r < data_len; ++r)
		if (s_reference[ref_pos + r] != s_data[data_pos + r])
			break;

	return r;
}

// ****************************************************************************
int my_hash(seq_t::iterator p, int len)
{
	uint32_t r = 0;

	for (int i = 0; i < len; ++i)
	{
		if (*(p + i) == sym_N)
			return HT_FAIL;

		r = r * 127 + *(p + i);
	}

	return r & ht_mask;
}

// ****************************************************************************
bool load_data(int argc, char **argv)
{
	if (argc < 3)
	{
		cerr << "Usage: ani-entropy <ref_name> <data_name>\n";
		return false;
	}

	if (!load_file(string(argv[1]), s_reference, n_reference))
	{
		cerr << "Cannot load " << argv[1] << endl;
		return false;
	}

	if (!load_file(string(argv[2]), s_data, n_data))
	{
		cerr << "Cannot load " << argv[2] << endl;
		return false;
	}

	duplicate_rev_comp(s_reference);

	return true;
}

// ****************************************************************************
void parse()
{
	v_parsing.clear();

	int data_size = (int) s_data.size();
	int ref_pred_pos = -data_size;
	int cur_lit_run_len = 0;

	for (int i = 0; i + MIN_MATCH_LEN < data_size;)
	{
		auto h = my_hash(s_data.begin() + i, MIN_MATCH_LEN);
		uint32_t best_pos = 0;
		uint32_t best_len = 0;

		for (; ht[h] != HT_EMPTY; h = (h + 1) & ht_mask)
		{
			uint32_t matching_len = equal_len(ht[h], i);
			
			if (matching_len < MIN_MATCH_LEN)
				continue;
			if (matching_len < MIN_DISTANT_MATCH_LEN && abs(ht[h] - ref_pred_pos) > CLOSE_DIST)
				continue;

			if (matching_len > best_len)
			{
				best_len = matching_len;
				best_pos = ht[h];
			}
		}

		if (best_len >= MIN_MATCH_LEN)
		{
			v_parsing.push_back(CFactor(flag_t::match, best_pos, best_len, 0));
			i += best_len;
			ref_pred_pos = best_pos + best_len;
			cur_lit_run_len = 0;
		}
		else
		{
			v_parsing.push_back(CFactor(flag_t::literal, 0, 0, s_data[i]));
			++i;
			++ref_pred_pos;
			++cur_lit_run_len;
		}

		if (cur_lit_run_len > MAX_LIT_RUN_IN_MATCH)
			ref_pred_pos = -data_size;
	}
}

// ****************************************************************************
void parsing_postprocess()
{
	vector<CFactor> new_parsing;

	int lit_run_len = 0;
	int ref_pred_pos = -(int)s_data.size();

	for (auto &x : v_parsing)
	{
		if (x.flag == flag_t::literal)
			++lit_run_len;
		else if (x.flag == flag_t::match)
		{
			if (lit_run_len)
			{
				ref_pred_pos += lit_run_len;
				new_parsing.push_back(CFactor(flag_t::run_literals, 0, lit_run_len, 0));
				lit_run_len = 0;
			}

			if (abs(ref_pred_pos - x.offset) <= CLOSE_DIST)
			{
				x.flag = flag_t::match_close;
				ref_pred_pos = x.offset;
			}
			else
			{
				x.flag = flag_t::match_distant;
				ref_pred_pos = x.offset;
			}

			new_parsing.push_back(x);
			ref_pred_pos += x.len;
		}
	}

	swap(v_parsing, new_parsing);
}

// ****************************************************************************
void export_parsing()
{
	FILE *f = fopen("parsing.log", "wb");

	if (!f)
	{
		cerr << "Cannot open log file\n";
		exit(0);
	}

	for (auto &x : v_parsing)
	{
		if (x.flag == flag_t::literal)
			fprintf(f, "Literal: %c\n", x.symbol);
		else if (x.flag == flag_t::run_literals)
			fprintf(f, "Run-lit: %d\n", x.len);
		else if (x.flag == flag_t::match)
			fprintf(f, "Match  : Off:%8d  Len: %8d\n", x.offset, x.len);
		else if (x.flag == flag_t::match_close)
			fprintf(f, "Match-c: Off:%8d  Len: %8d\n", x.offset, x.len);
		else if (x.flag == flag_t::match_distant)
			fprintf(f, "Match-d: Off:%8d  Len: %8d\n", x.offset, x.len);
	}

	fclose(f);
}

// ****************************************************************************
void prepare_ht()
{
	uint32_t x = (uint32_t) (s_reference.size() / ht_max_fill_factor);

	while (x & (x - 1))
		x &= x - 1;

	ht_size = 2 * x;
	ht_mask = ht_size - 1;

	ht.clear();
	ht.resize(ht_size, HT_EMPTY);

	for (size_t i = 0; i + MIN_MATCH_LEN < s_reference.size(); ++i)
	{
		auto ht_idx = my_hash(s_reference.begin() + i, MIN_MATCH_LEN);

		if (ht_idx != HT_FAIL)
		{
			while (ht[ht_idx] != HT_EMPTY)
				ht_idx = (ht_idx + 1) & ht_mask;

			ht[ht_idx] = (uint32_t) i;
		}
	}
}

// ****************************************************************************
void calc_ani()
{
	int ref_len = (int)s_reference.size() - n_reference * CLOSE_DIST;
	int data_len = (int)s_data.size() - n_data * CLOSE_DIST;
	int n_sym_in_matches = 0;
	int n_sym_in_literals = 0;
	int last_run_len = 0;

	for (auto &x : v_parsing)
	{
		if (x.flag == flag_t::match_distant)
		{
			n_sym_in_matches += x.len;
			last_run_len = 0;
		}
		else if (x.flag == flag_t::match_close)
		{
			n_sym_in_matches += x.len;
			n_sym_in_literals += last_run_len;
			last_run_len = 0;
		}
		else if (x.flag == flag_t::run_literals)
		{
			if (x.len <= LONG_LITERAL_RUN_LEN)
				last_run_len = x.len;
			else
				last_run_len = 0;
		}
	}
	
	cout << "Ref size          : " << ref_len << endl;
	cout << "Data size         : " << data_len << endl;
	cout << "No_sym_in_matches : " << n_sym_in_matches << endl;
	cout << "No_sym_in_literals: " << n_sym_in_literals << endl;

	cout << "ANI est: " << std::setw(7) << std::setprecision(5) << 100.0 * (double) n_sym_in_matches / (n_sym_in_matches + n_sym_in_literals) << endl;
}

// ****************************************************************************
int main(int argc, char **argv)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	cout << "*** Loading...";
	if (!load_data(argc, argv))
	{
		return 0;
	}
	cout << " ok\n";

	cout << "*** HT preparation...";
	prepare_ht();
	cout << " ok\n";

	cout << "*** Parsing...";
	parse();
	parsing_postprocess();
	export_parsing();

	cout << " ok\n";

	calc_ani();

	cout << "*** End of work\n";

	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "Running time: " << time_span.count() << " seconds.\n";

//	getchar();

	return 0;
}

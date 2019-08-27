#include "base_worker.h"
#include <algorithm>
#include <iostream>
#include <functional>

extern int MIN_MATCH_LEN;
extern int MIN_CLOSE_MATCH_LEN;
extern int MIN_DISTANT_MATCH_LEN;
extern int CLOSE_DIST;
extern int MAX_LIT_RUN_IN_MATCH;
extern double MIN_COVERAGE;
extern int MIN_REGION_LEN;
extern int APPROX_WINDOW;
extern int APPROX_MISMATCHES;
extern int APPROX_RUNLEN;


using namespace std;

// ****************************************************************************
BaseWorker::BaseWorker() {
	fill_n(codes, 256, 4);
	codes['A'] = 0;
	codes['C'] = 1;
	codes['G'] = 2;
	codes['T'] = 3;

	hts_mask = (int)(1u << (2 * (MIN_DISTANT_MATCH_LEN - MIN_MATCH_LEN))) - 1;
}


// ****************************************************************************
void BaseWorker::export_parsing()
{
	FILE *f = fopen("parsing.log", "wb");

	if (!f)
	{
		cerr << "Cannot open log file\n";
		exit(0);
	}

	setvbuf(f, nullptr, _IOFBF, 32 << 20);

	int pred_data_pos = 0;

	for (auto &x : v_parsing)
	{
		if (pred_data_pos != x.query_pos)
			fprintf(f, "*******\n");
		fprintf(f, "Data pos: %8d   ", x.query_pos);
		if (x.flag == flag_t::literal)
			fprintf(f, "Literal    : %c\n", x.symbol);
		else if (x.flag == flag_t::run_literals)
			fprintf(f, "Run-lit    : %d\n", x.len);
		else if (x.flag == flag_t::match)
			fprintf(f, "Match      : Off:%8d  Len: %8d\n", x.ref_pos, x.len);
		else if (x.flag == flag_t::match_close)
			fprintf(f, "Match-close: Off:%8d  Len: %8d\n", x.ref_pos, x.len);
		else if (x.flag == flag_t::match_distant)
			fprintf(f, "Match-dist : Off:%8d  Len: %8d\n", x.ref_pos, x.len);
		else if (x.flag == flag_t::match_literal)
			fprintf(f, "Match-lit  : Off:%8d  Len: %8d\n", x.ref_pos, x.len);

		pred_data_pos += x.len;
	}

	f = fopen("parsing.log2", "wb");

	if (!f)
	{
		cerr << "Cannot open log file\n";
		exit(0);
	}

	setvbuf(f, nullptr, _IOFBF, 32 << 20);

	vector<pair<int, int>> v_matches;
	int cur_match_len = 0;
	int cur_match_lit = 0;
	int n_lit = 0;

	for (auto x : v_parsing)
	{
		if (x.flag == flag_t::match_distant)
		{
			if (cur_match_len)
				v_matches.emplace_back(make_pair(cur_match_len, cur_match_lit));

			cur_match_len = x.len;
			cur_match_lit = 0;
			n_lit = 0;
		}
		else if (x.flag == flag_t::match_close)
		{
			cur_match_len += x.len;
			cur_match_lit += n_lit;
			n_lit = 0;
		}
		else if (x.flag == flag_t::run_literals)
		{
			n_lit += x.len;
		}
	}

	if (cur_match_len)
		v_matches.emplace_back(make_pair(cur_match_len, cur_match_lit));

	std::sort(v_matches.begin(), v_matches.end(), std::greater<std::pair<int, int>>());

	for (auto x : v_matches)
		fprintf(f, "%d : %d\n", x.first, x.second);

	fclose(f);
}

// ****************************************************************************
bool BaseWorker::load_file(const string &file_name, seq_t &seq, uint32_t &n_parts, int separator)
{
	seq.clear();

	FILE *f = fopen(file_name.c_str(), "rb");
	if (!f)
		return false;

	setvbuf(f, nullptr, _IOFBF, 32 << 20);

	int c;
	bool is_comment = false;

	n_parts = 0;

	while ((c = getc(f)) != EOF)
	{
		if (c == '>')
			is_comment = true;
		else
		{
			if (c == '\n' || c == '\r')
			{
				if (is_comment)
				{
					is_comment = false;
					if (!seq.empty())
						for (int i = 0; i < CLOSE_DIST; ++i)
							seq.emplace_back(separator);
					++n_parts;
				}
			}
			else if (!is_comment)
				seq.emplace_back((uint8_t)c);
		}
	}

	fclose(f);

	return true;
}


// ****************************************************************************
void BaseWorker::duplicate_rev_comp(seq_t &seq)
{
	int size = (int)seq.size();

	seq.reserve(2 * size);

	for (int i = size - 1; i >= 0; --i)
	{
		if (seq[i] == sym_A)
			seq.emplace_back(sym_T);
		else if (seq[i] == sym_C)
			seq.emplace_back(sym_G);
		else if (seq[i] == sym_G)
			seq.emplace_back(sym_C);
		else if (seq[i] == sym_T)
			seq.emplace_back(sym_A);
		else
			seq.emplace_back(seq[i]);
	}
}

// ****************************************************************************
void BaseWorker::prepare_kmers(vector<pair<int64_t, int>> &v_kmers, const seq_t &seq, int len, bool store_all)
{
	v_kmers.clear();
	v_kmers.reserve(seq.size());

	uint64_t k = 0u;
	uint64_t mask = (~0ull) >> (64 - 2 * len);
	int k_len = 0;
	int seq_size = (int)seq.size();

	for (int i = 0; i < seq_size; ++i)
	{
		if (codes[seq[i]] == 4)
		{
			k = 0ull;
			k_len = 0;
		}
		else
		{
			k <<= 2;
			++k_len;
			k += codes[seq[i]];
			k &= mask;
		}

		if (i >= len - 1)
		{
			if (k_len >= len)
				v_kmers.emplace_back(make_pair(k, i + 1 - len));
			else if (store_all)
				v_kmers.emplace_back(make_pair(-1, i + 1 - len));
		}
	}

	if (store_all)
		for (int i = 0; i < len - 1; ++i)
			v_kmers.emplace_back(make_pair(-1, 0));
}

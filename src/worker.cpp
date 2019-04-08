#include "worker.h"
#include <iostream>
#include <iomanip>
#include "xmmintrin.h"
#include <nmmintrin.h>
#include <algorithm>

// ****************************************************************************
int CWorker::equal_len(int ref_pos, int data_pos, int starting_pos)
{
	int r;
	int ref_len = (int)s_reference.size();
	int data_len = (int)s_data.size();

	for (r = starting_pos; ref_pos + r < ref_len && data_pos + r < data_len; ++r)
		if (s_reference[ref_pos + r] != s_data[data_pos + r])
			break;

	return r;
}

// ****************************************************************************
int CWorker::my_hash_l(seq_t::iterator p, int len)
{
	uint64_t r = 0;

	for (int i = 0; i < len; ++i)
	{
		r <<= 2;
		switch (*(p + i))
		{
		case 'A': r += 0ull;	break;
		case 'C': r += 1ull;	break;
		case 'G': r += 2ull;	break;
		case 'T': r += 3ull;	break;
		default: return HT_FAIL;
		}
	}

	r ^= r >> 33;
	r *= 0xff51afd7ed558ccdLL;
	r ^= r >> 33;
	r *= 0xc4ceb9fe1a85ec53LL;
	r ^= r >> 33;

	return (int) (r & htl_mask);
}

// ****************************************************************************
int CWorker::hash_mm(uint64_t x, int mask)
{
	x ^= x >> 33;
	x *= 0xff51afd7ed558ccdLL;
	x ^= x >> 33;
	x *= 0xc4ceb9fe1a85ec53LL;
	x ^= x >> 33;

	return (int)(x & mask);
}

// ****************************************************************************
int CWorker::my_hash_s(seq_t::iterator p, int len)
{
	uint64_t r = 0;

	for (int i = 0; i < len; ++i)
	{
		r <<= 2;
		switch (*(p + i))
		{
		case 'A': r += 0;	break;
		case 'C': r += 1;	break;
		case 'G': r += 2;	break;
		case 'T': r += 3;	break;
		default: return HT_FAIL;
		}
	}

	return r;
}

// ****************************************************************************
bool CWorker::load_data(string fn_ref, string fn_data)
{
	if (!load_file(fn_ref, s_reference, n_reference))
	{
		cerr << "Error: Cannot load " + fn_ref + "\n";
		return false;
	}

	if (!load_file(fn_data, s_data, n_data))
	{
		cerr << "Error: Cannot load " + fn_data + "\n";
		return false;
	}

	duplicate_rev_comp(s_reference);

	return true;
}

// ****************************************************************************
void CWorker::prefetch(int pos)
{
#ifdef _WIN32
	_mm_prefetch((const char*)(s_reference.data() + pos), _MM_HINT_T0);
#else
	__builtin_prefetch(s_reference.data() + pos);
#endif
}

// ****************************************************************************
void CWorker::parse()
{
	v_parsing.clear();

	int data_size = (int)s_data.size();
	int ref_pred_pos = -data_size;
	int cur_lit_run_len = 0;

	for (int i = 0; i + MIN_MATCH_LEN < data_size;)
	{
		uint32_t best_pos = 0;
		uint32_t best_len = 0;

		if (ref_pred_pos < 0)	// Look for long match
		{
			auto h = my_hash_l(s_data.begin() + i, MIN_DISTANT_MATCH_LEN);

			for (; htl[h] != HT_EMPTY; h = (h + 1) & htl_mask)
			{
				uint32_t matching_len = equal_len(htl[h], i);

				if (matching_len < MIN_DISTANT_MATCH_LEN)
					continue;

				if (matching_len > best_len)
				{
					best_len = matching_len;
					best_pos = htl[h];
				}
			}
		}
		else
		{
			auto h = my_hash_s(s_data.begin() + i, MIN_MATCH_LEN);

			if (h != HT_FAIL)
			{
				int bucket_size = (int)hts[h].size();
				auto &bucket = hts[h];

				for (int j = 0; j < min(3, bucket_size); ++j)
					prefetch(bucket[j]);

				for (int j = 0; j < bucket_size; ++j)
				{
					if (j + 3 < bucket_size)
						prefetch(bucket[j + 3]);

					auto pos = bucket[j];
					uint32_t matching_len = equal_len(pos, i, MIN_MATCH_LEN);

					if (matching_len < MIN_MATCH_LEN)
						continue;
					if (matching_len < MIN_DISTANT_MATCH_LEN && abs(pos - ref_pred_pos) > CLOSE_DIST)
						continue;

					if (matching_len > best_len)
					{
						best_len = matching_len;
						best_pos = pos;
					}
				}
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
void CWorker::parsing_postprocess()
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
void CWorker::export_parsing()
{
	FILE *f = fopen("parsing.log", "wb");

	if (!f)
	{
		cerr << "Cannot open log file\n";
		exit(0);
	}

	setvbuf(f, nullptr, _IOFBF, 32 << 20);

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
void CWorker::prepare_ht_long()
{
	prepare_kmers(v_kmers_l, s_reference, MIN_DISTANT_MATCH_LEN);

	uint32_t x = (uint32_t)(v_kmers_l.size() / htl_max_fill_factor);

	while (x & (x - 1))
		x &= x - 1;

	htl_size = 2 * x;
	htl_mask = htl_size - 1;

	htl.clear();
	htl.resize(htl_size, HT_EMPTY);

	const int pf_dist = 16;

	for (int i = 0; i + pf_dist < (int)v_kmers_l.size(); ++i)
	{
		prefetch_htl(hash_mm(v_kmers_l[i + pf_dist].first, htl_mask));

		auto ht_idx = hash_mm(v_kmers_l[i].first, htl_mask);

		while (htl[ht_idx] != HT_EMPTY)
			ht_idx = (ht_idx + 1) & htl_mask;

		htl[ht_idx] = v_kmers_l[i].second;
	}

	
	for (int i = max((int)v_kmers_l.size() - pf_dist, 0); i < (int)v_kmers_l.size(); ++i)
	{
		auto ht_idx = hash_mm(v_kmers_l[i].first, htl_mask);

		while (htl[ht_idx] != HT_EMPTY)
			ht_idx = (ht_idx + 1) & htl_mask;

		htl[ht_idx] = v_kmers_l[i].second;
	}
	

/*

	for (size_t i = 0; i + MIN_DISTANT_MATCH_LEN < s_reference.size(); ++i)
	{
		auto ht_idx = my_hash_l(s_reference.begin() + i, MIN_DISTANT_MATCH_LEN);

		if (ht_idx != HT_FAIL)
		{
			while (htl[ht_idx] != HT_EMPTY)
				ht_idx = (ht_idx + 1) & htl_mask;

			htl[ht_idx] = (uint32_t)i;
		}
	}
	*/
}

// ****************************************************************************
void CWorker::prefetch_hts(int pos)
{
#ifdef _WIN32
	_mm_prefetch((const char*)(hts.data() + pos), _MM_HINT_T0);
#else
	__builtin_prefetch(hts.data() + pos);
#endif
}

// ****************************************************************************
void CWorker::prefetch_htl(int pos)
{
#ifdef _WIN32
	_mm_prefetch((const char*)(htl.data() + pos), _MM_HINT_T0);
#else
	__builtin_prefetch(htl.data() + pos);
#endif
}

// ****************************************************************************
// !!! Opt: mozna hasze buforowac w malym wektorze i zapisywac do HT po zrobieniu prefetcha
void CWorker::prepare_ht_short()
{
	uint32_t ht_size = 1u << (2 * MIN_MATCH_LEN);
	uint32_t ht_mask = ht_size - 1u;

	hts.clear();
	hts.resize(ht_size);

	prepare_kmers(v_kmers_s, s_reference, MIN_MATCH_LEN);

	const int pf_dist = 32;
	
	for (int i = 0; i + pf_dist < (int)v_kmers_s.size(); ++i)
	{
		prefetch_hts(v_kmers_s[i + pf_dist].first);
		hts[v_kmers_s[i].first].push_back(v_kmers_s[i].second);
	}

	for (int i = max((int)v_kmers_s.size() - pf_dist, 0); i < (int)v_kmers_s.size(); ++i)
		hts[v_kmers_s[i].first].push_back(v_kmers_s[i].second);
	

/*	for (size_t i = 0; i + MIN_MATCH_LEN < s_reference.size(); ++i)
	{
		auto ht_idx = my_hash_s(s_reference.begin() + i, MIN_MATCH_LEN);

		if (ht_idx != HT_FAIL)
			hts[ht_idx].push_back(i);
	}*/
}

// ****************************************************************************
void CWorker::prepare_pf()
{
	prepare_kmers(v_kmers_s, s_data, MIN_MATCH_LEN);
	prepare_kmers(v_kmers_l, s_data, MIN_DISTANT_MATCH_LEN);
}

// ****************************************************************************
void CWorker::calc_ani(CResults &res)
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

	res.ref_size = ref_len / 2;
	res.query_size = data_len;
	res.sym_in_literals = n_sym_in_literals;
	res.sym_in_matches = n_sym_in_matches;
	res.coverage = (double) (n_sym_in_literals + n_sym_in_matches) / data_len;
	res.ani = (double)n_sym_in_matches / (n_sym_in_matches + n_sym_in_literals);
}

// ****************************************************************************
bool CWorker::load_file(const string &file_name, seq_t &seq, uint32_t &n_parts)
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
							seq.push_back(sym_N);
					++n_parts;
				}
			}
			else if (!is_comment)
				seq.push_back((uint8_t)c);
		}
	}

	fclose(f);

	return true;
}

// ****************************************************************************
void CWorker::duplicate_rev_comp(seq_t &seq)
{
	int size = (int)seq.size();

	seq.reserve(2 * size);

	for (int i = size - 1; i >= 0; --i)
	{
		if (seq[i] == sym_N)
			seq.push_back(sym_N);
		else if (seq[i] == sym_A)
			seq.push_back(sym_T);
		else if (seq[i] == sym_C)
			seq.push_back(sym_G);
		else if (seq[i] == sym_G)
			seq.push_back(sym_C);
		else if (seq[i] == sym_T)
			seq.push_back(sym_A);
	}
}

// ****************************************************************************
void CWorker::clear()
{
	s_reference.clear();
	s_data.clear();

	n_reference = 0;
	n_data = 0;

	htl_size = 0;
	htl_mask = 0;

	htl.clear();
	htl.shrink_to_fit();

	hts.clear();
	hts.shrink_to_fit();

	v_parsing.clear();
	v_parsing.shrink_to_fit();
}

// ****************************************************************************
void CWorker::prepare_kmers(vector<pair<uint64_t, int>> &v_kmers, const seq_t &seq, int len)
{
	v_kmers.clear();
	v_kmers.reserve(seq.size());

	uint64_t k = 0u;
	uint64_t mask = (~0ull) >> (64 - 2 * len);
	int k_len = 0;

	for (int i = 0; i < seq.size(); ++i)
	{
		k <<= 2;
		++k_len;
		switch (seq[i])
		{
		case 'A': 
			k += 0ull;		break;
		case 'C':
			k += 1ull;		break;
		case 'G':
			k += 2ull;		break;
		case 'T':
			k += 3ull;		break;
		default:
			k = 0ull;
			k_len = 0;
		}

		k &= mask;

		if (k_len >= len)
			v_kmers.emplace_back(make_pair(k, i + 1 - len));
	}
}

// EOF
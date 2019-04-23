#include "worker.h"
#include <iostream>
#include <iomanip>
#include "xmmintrin.h"
#include <nmmintrin.h>
#include <immintrin.h>
#include <algorithm>

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


// ****************************************************************************
CWorker::CWorker()
{
	fill_n(codes, 256, 4);
	codes['A'] = 0;
	codes['C'] = 1;
	codes['G'] = 2;
	codes['T'] = 3;
}

// ****************************************************************************
int CWorker::equal_len(int ref_pos, int data_pos, int starting_pos)
{
	int r;
	int ref_len = (int)s_reference.size();
	int data_len = (int)s_data.size();
	int max_r = min(ref_len - ref_pos, data_len - data_pos);

	uint8_t *p_ref = s_reference.data() + ref_pos + starting_pos;
	uint8_t *p_data = s_data.data() + data_pos + starting_pos;

	for (r = starting_pos; r < max_r; ++r)
		if (*p_ref++ != *p_data++)
			break;

	return r;
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
bool CWorker::load_data(string fn_ref, string fn_data)
{
	if (!load_file(fn_ref, s_reference, n_reference, sym_N1))
	{
		cerr << "Error: Cannot load " + fn_ref + "\n";
		return false;
	}

	if (!load_file(fn_data, s_data, n_data, sym_N2))
	{
		cerr << "Error: Cannot load " + fn_data + "\n";
		return false;
	}

	duplicate_rev_comp(s_reference);

	return true;
}

// ****************************************************************************
void CWorker::swap_data()
{
	s_reference.resize(s_reference.size() / 2);
	swap(s_data, s_reference);
	swap(n_data, n_reference);

	duplicate_rev_comp(s_reference);
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
void CWorker::compare_ranges(int data_start_pos, int ref_start_pos, int len, bool backward = false)
{
	int r_len = 0;
	bool is_matching = false;
	flag_t flag = backward ? flag_t::match_distant : flag_t::match_close;

	for (int j = 0; j < len; ++j)
	{
		if (s_reference[ref_start_pos + j] == s_data[data_start_pos + j])
		{
			if (is_matching)
				r_len++;
			else
			{
				if(r_len)
					v_parsing.emplace_back(CFactor(data_start_pos+j-r_len, flag_t::run_literals, 0, r_len, 0));
				r_len = 1;
				is_matching = true;
			}
		}
		else
		{
			if (is_matching)
			{
				v_parsing.emplace_back(CFactor(data_start_pos + j-r_len, flag, ref_start_pos + j - r_len, r_len, 0));
				r_len = 1;
				is_matching = false;
				flag = flag_t::match_close;
			}
			else
				++r_len;
		}
	}

	if (is_matching) 
		v_parsing.emplace_back(CFactor(data_start_pos + len - r_len, flag, ref_start_pos + len - r_len, r_len, 0));
	else if(r_len)
		v_parsing.emplace_back(CFactor(data_start_pos + len - r_len, flag_t::run_literals, 0, r_len, 0));
}

// ****************************************************************************
int CWorker::try_extend_forward(int data_start_pos, int ref_start_pos)
{
	int data_size = (int)s_data.size();
	int ref_size = (int)s_reference.size();

	int approx_ext;
	int no_missmatches = 0;
	int last_match = 0;
	vector<int> window(APPROX_WINDOW, 0);

	for (approx_ext = 0; data_start_pos + approx_ext < data_size && ref_start_pos + approx_ext < ref_size; ++approx_ext)
	{
		bool is_missmatch = s_data[data_start_pos + approx_ext] != s_reference[ref_start_pos + approx_ext];
		no_missmatches -= window[approx_ext % APPROX_WINDOW];
		window[approx_ext % APPROX_WINDOW] = is_missmatch;
		no_missmatches += is_missmatch;

		if (!is_missmatch)
			last_match = approx_ext + 1;

		if (no_missmatches > APPROX_MISMATCHES)
			break;
	}

	return last_match;
}

// ****************************************************************************
int CWorker::try_extend_forward2(int data_start_pos, int ref_start_pos)
{
	int data_size = (int)s_data.size();
	int ref_size = (int)s_reference.size();

	int approx_ext;
	int no_missmatches = 0;
	int last_run_match = 0;
	vector<int> window(APPROX_WINDOW, 0);
	int match_run_len = APPROX_RUNLEN;

	for (approx_ext = 0; data_start_pos + approx_ext < data_size && ref_start_pos + approx_ext < ref_size; ++approx_ext)
	{
		bool is_missmatch = s_data[data_start_pos + approx_ext] != s_reference[ref_start_pos + approx_ext];
		no_missmatches -= window[approx_ext % APPROX_WINDOW];
		window[approx_ext % APPROX_WINDOW] = is_missmatch;
		no_missmatches += is_missmatch;

		if (!is_missmatch)
		{
			if(++match_run_len >= APPROX_RUNLEN)
				last_run_match = approx_ext + 1;
		}
		else
			match_run_len = 0;

		if (no_missmatches > APPROX_MISMATCHES)
			break;
	}

	return last_run_match;
}

// ****************************************************************************
int CWorker::try_extend_backward(int data_start_pos, int ref_start_pos, int max_len)
{
	int approx_ext;
	int no_missmatches = 0;
	int last_match = 0;
	vector<int> window(APPROX_WINDOW, 0);
	
	for (approx_ext = 0; data_start_pos - approx_ext > 0 && ref_start_pos - approx_ext > 0 && approx_ext < max_len; ++approx_ext)
	{
		bool is_missmatch = s_data[data_start_pos - approx_ext - 1] != s_reference[ref_start_pos - approx_ext - 1];
		no_missmatches -= window[approx_ext % APPROX_WINDOW];
		window[approx_ext % APPROX_WINDOW] = is_missmatch;
		no_missmatches += is_missmatch;

		if (!is_missmatch)
			last_match = approx_ext + 1;

		if (no_missmatches > APPROX_MISMATCHES)
			break;
	}

	return last_match;
}

// ****************************************************************************
int CWorker::try_extend_backward2(int data_start_pos, int ref_start_pos, int max_len)
{
	int approx_ext;
	int no_missmatches = 0;
	int last_run_match = 0;
	vector<int> window(APPROX_WINDOW, 0);
	int match_run_len = APPROX_RUNLEN;

	for (approx_ext = 0; data_start_pos - approx_ext > 0 && ref_start_pos - approx_ext > 0 && approx_ext < max_len; ++approx_ext)
	{
		bool is_missmatch = s_data[data_start_pos - approx_ext - 1] != s_reference[ref_start_pos - approx_ext - 1];
		no_missmatches -= window[approx_ext % APPROX_WINDOW];
		window[approx_ext % APPROX_WINDOW] = is_missmatch;
		no_missmatches += is_missmatch;

		if (!is_missmatch)
		{
			if (++match_run_len >= APPROX_RUNLEN)
				last_run_match = approx_ext + 1;
		}
		else
			match_run_len = 0;

		if (no_missmatches > APPROX_MISMATCHES)
			break;
	}

	return last_run_match;
}

// ****************************************************************************
void CWorker::parse()
{
	v_parsing.clear();

	int data_size = (int)s_data.size();
	int ref_pred_pos = -data_size;
	int cur_lit_run_len = 0;

	const int pf_dist_l = 8;
	const int pf_dist_s = 12;

	int i;

	int prev_region_start = -1;
	int prev_region_end = 0;

	for (i = 0; i + MIN_MATCH_LEN < data_size;)
	{
		int best_pos = 0;
		int best_len = 0;
		int h;

		if (ref_pred_pos < 0)
		{
			// Look for long match
			if (i + pf_dist_l < data_size && v_kmers_l[i + pf_dist_l].first >= 0)
				prefetch_htl(hash_mm(v_kmers_l[i + pf_dist_l].first, htl_mask));

			if (v_kmers_l[i].first >= 0)
			{
				h = hash_mm(v_kmers_l[i].first, htl_mask);

				for (; htl[h] != HT_EMPTY; h = (h + 1) & htl_mask)
				{
					int matching_len = equal_len(htl[h], i);

					if (matching_len < MIN_DISTANT_MATCH_LEN)
						continue;

					if (matching_len > best_len)
					{
						best_len = matching_len;
						best_pos = htl[h];
					}
				}
			}
		}
		else
		{
			// Look for short but close match
			if (i + pf_dist_s < data_size && v_kmers_s[i + pf_dist_s].first >= 0)
				prefetch_hts((int)v_kmers_s[i + pf_dist_s].first);

			auto h = v_kmers_s[i].first;

			if (h != HT_FAIL)
			{
				int bucket_size = (int)hts[h].size();
				auto &bucket = hts[h];
				const int pf_dist = 4;

				for (int j = 0; j < min(pf_dist, bucket_size); ++j)
					prefetch(bucket[j]);

				int best_close_len = 0;
				int best_close_pos = 0;

				for (int j = 0; j < bucket_size; ++j)
				{
					if (j + pf_dist < bucket_size)
						prefetch(bucket[j + pf_dist]);

					auto pos = bucket[j];
					int matching_len = equal_len(pos, i, MIN_MATCH_LEN);

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
			if (cur_lit_run_len)
			{
				if (ref_pred_pos >= 0)
					compare_ranges(i - cur_lit_run_len, ref_pred_pos - cur_lit_run_len, cur_lit_run_len);
				else
					v_parsing.emplace_back(CFactor(i- cur_lit_run_len, flag_t::run_literals, 0, cur_lit_run_len, 0));
			}

			flag_t flag = flag_t::match_distant;

			if (abs(best_pos - ref_pred_pos) <= CLOSE_DIST)
			{
				v_parsing.emplace_back(CFactor(i, flag_t::match_close, best_pos, best_len, 0));
			}
			else
			{
				// Remove previous region if too short
				if (prev_region_start >= 0 && prev_region_end - prev_region_start < MIN_REGION_LEN)
				{
					while (!v_parsing.empty() && v_parsing.back().data_pos >= prev_region_start)
						v_parsing.pop_back();
					int run_len = i - prev_region_start;

					while (!v_parsing.empty() && v_parsing.back().flag == flag_t::run_literals)
					{
						run_len += v_parsing.back().len;
						v_parsing.pop_back();
					}

					v_parsing.emplace_back(CFactor(i - run_len, flag_t::run_literals, 0, run_len, 0));
					prev_region_start = -1;
				}

				if (!v_parsing.empty() && v_parsing.back().flag == flag_t::run_literals)
				{
					int approx_pred = try_extend_backward2(i, best_pos, v_parsing.back().len);
					if (approx_pred)
					{
						v_parsing.back().len -= approx_pred;
						if (v_parsing.back().len == 0)
							v_parsing.pop_back();
						compare_ranges(i - approx_pred, best_pos - approx_pred, approx_pred, true);
						flag = flag_t::match_close;
						prev_region_start = -1;
					}
				}

				v_parsing.emplace_back(CFactor(i, flag, best_pos, best_len, 0));
				if (flag == flag_t::match_distant)
					prev_region_start = i;

				if(prev_region_start < 0)
					for(int j = (int) v_parsing.size() - 1; j >= 0; --j)
						if (v_parsing[j].flag == flag_t::match_distant)
						{
							prev_region_start = v_parsing[j].data_pos;
							break;
						}
			}

			i += best_len;
			ref_pred_pos = best_pos + best_len;
			cur_lit_run_len = 0;

			int approx_ext = try_extend_forward2(i, ref_pred_pos);
			compare_ranges(i, ref_pred_pos, approx_ext);

			i += approx_ext;
			ref_pred_pos += approx_ext;

			prev_region_end = i;
		}
		else
		{
			++i;
			++ref_pred_pos;
			++cur_lit_run_len;
		}

		if (cur_lit_run_len > MAX_LIT_RUN_IN_MATCH)
			ref_pred_pos = -data_size;
	}

	if(ref_pred_pos < 0)
		v_parsing.emplace_back(CFactor(i - cur_lit_run_len, flag_t::run_literals, 0, cur_lit_run_len + (data_size - i), 0));
	else
		compare_ranges(i - cur_lit_run_len, ref_pred_pos - cur_lit_run_len - MIN_MATCH_LEN, cur_lit_run_len + (data_size - i));
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

	int pred_data_pos = 0;

	for (auto &x : v_parsing)
	{
		if (pred_data_pos != x.data_pos)
			fprintf(f, "*******\n");
		fprintf(f, "Data pos: %8d   ", x.data_pos);
		if (x.flag == flag_t::literal)
			fprintf(f, "Literal    : %c\n", x.symbol);
		else if (x.flag == flag_t::run_literals)
			fprintf(f, "Run-lit    : %d\n", x.len);
		else if (x.flag == flag_t::match)
			fprintf(f, "Match      : Off:%8d  Len: %8d\n", x.offset, x.len);
		else if (x.flag == flag_t::match_close)
			fprintf(f, "Match-close: Off:%8d  Len: %8d\n", x.offset, x.len);
		else if (x.flag == flag_t::match_distant)
			fprintf(f, "Match-dist : Off:%8d  Len: %8d\n", x.offset, x.len);
		else if(x.flag == flag_t::match_literal)
			fprintf(f, "Match-lit  : Off:%8d  Len: %8d\n", x.offset, x.len);

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

	sort(v_matches.begin(), v_matches.end(), greater<pair<int, int>>());

	for (auto x : v_matches)
		fprintf(f, "%d : %d\n", x.first, x.second);

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
void CWorker::prepare_ht_short()
{
	uint32_t ht_size = 1u << (2 * MIN_MATCH_LEN);
//	uint32_t ht_mask = ht_size - 1u;

	hts.clear();
	hts.resize(ht_size);

	prepare_kmers(v_kmers_s, s_reference, MIN_MATCH_LEN, true);
	prepare_kmers(v_kmers_l, s_reference, MIN_DISTANT_MATCH_LEN, true);

	const int pf_dist = 32;
	
	for (int i = 0; i + pf_dist < (int)v_kmers_s.size(); ++i)
	{
		if(v_kmers_s[i + pf_dist].first >= 0)
			prefetch_hts((int) v_kmers_s[i + pf_dist].first);
		if(v_kmers_s[i].first >= 0)
			hts[v_kmers_s[i].first].emplace_back(v_kmers_s[i].second);
	}

	for (int i = max((int)v_kmers_s.size() - pf_dist, 0); i < (int)v_kmers_s.size(); ++i)
		if (v_kmers_s[i].first >= 0)
			hts[v_kmers_s[i].first].emplace_back(v_kmers_s[i].second);
}

// ****************************************************************************
void CWorker::prepare_pf()
{
	prepare_kmers(v_kmers_s, s_data, MIN_MATCH_LEN, true);
	prepare_kmers(v_kmers_l, s_data, MIN_DISTANT_MATCH_LEN, true);
}

// ****************************************************************************
void CWorker::calc_ani(CResults &res, int mode)
{
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

	sort(v_matches.begin(), v_matches.end(), greater<pair<int, int>>());

	int ref_len = (int)s_reference.size() - n_reference * CLOSE_DIST;
	int data_len = (int)s_data.size() - n_data * CLOSE_DIST;
	int n_sym_in_matches = 0;
	int n_sym_in_literals = 0;
	
	for (auto x : v_matches)
		if (x.first + x.second >= MIN_REGION_LEN) // && (double) x.first / (x.first + x.second) > 0.5
		{
			n_sym_in_matches += x.first;
			n_sym_in_literals += x.second;
		}

	if (mode == 1)
	{
		res.ref_size = ref_len / 2;
		res.query_size = data_len;
	}
	res.sym_in_literals[mode] = n_sym_in_literals;
	res.sym_in_matches[mode] = n_sym_in_matches;
	res.coverage[mode] = (double)(n_sym_in_literals + n_sym_in_matches) / data_len;
	res.ani[mode] = (double)n_sym_in_matches / (n_sym_in_matches + n_sym_in_literals);

	if (res.coverage[mode] < MIN_COVERAGE)
		res.ani[mode] -= 0.4;
}

// ****************************************************************************
bool CWorker::load_file(const string &file_name, seq_t &seq, uint32_t &n_parts, int separator)
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
void CWorker::duplicate_rev_comp(seq_t &seq)
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
void CWorker::clear()
{
	s_reference.clear();
	s_data.clear();

	n_reference = 0;
	n_data = 0;

	htl_size = 0;
	htl_mask = 0;

	htl.clear();
	hts.clear();
	v_parsing.clear();
}

// ****************************************************************************
void CWorker::prepare_kmers(vector<pair<int64_t, int>> &v_kmers, const seq_t &seq, int len, bool store_all)
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

	if(store_all)
		for(int i = 0; i < len-1; ++i)
			v_kmers.emplace_back(make_pair(-1, 0));
}

// EOF
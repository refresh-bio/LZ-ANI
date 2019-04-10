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
void CWorker::parse()
{
	v_parsing.clear();

	int data_size = (int)s_data.size();
	int ref_pred_pos = -data_size;
	int cur_lit_run_len = 0;

	const int pf_dist_l = 8;
	const int pf_dist_s = 12;

	for (int i = 0; i + MIN_MATCH_LEN < data_size;)
	{
		uint32_t best_pos = 0;
		uint32_t best_len = 0;
		int h;

		if (ref_pred_pos < 0)	
		{
			// Look for long match
			if (i + pf_dist_l < data_size && v_kmers_l[i + pf_dist_l].first >= 0)
				prefetch_htl(hash_mm(v_kmers_l[i+pf_dist_l].first, htl_mask));

			if (v_kmers_l[i].first >= 0)
			{
				h = hash_mm(v_kmers_l[i].first, htl_mask);

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
		}
		else 
		{
			// Look for short but close match
			if (i + pf_dist_s < data_size && v_kmers_s[i + pf_dist_s].first >= 0)
				prefetch_hts((int) v_kmers_s[i + pf_dist_s].first);

			auto h = v_kmers_s[i].first;

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
			if (cur_lit_run_len)
			{
				// Tu spróbowaæ zmieniæ run litera³ów na matche/mismatche
				v_parsing.emplace_back(CFactor(flag_t::run_literals, 0, cur_lit_run_len, 0));
			}

			// !!! Ju¿ tu zdecydowaæ czy to close czy distant match
			v_parsing.emplace_back(CFactor(flag_t::match, best_pos, best_len, 0));
			i += best_len;
			ref_pred_pos = best_pos + best_len;
			cur_lit_run_len = 0;
		}
		else
		{
//			v_parsing.emplace_back(CFactor(flag_t::literal, 0, 0, s_data[i]));
			++i;
			++ref_pred_pos;
			++cur_lit_run_len;
		}

		if (cur_lit_run_len > MAX_LIT_RUN_IN_MATCH)
			ref_pred_pos = -data_size;
	}

//	if (cur_lit_run_len)
	v_parsing.emplace_back(CFactor(flag_t::run_literals, 0, cur_lit_run_len + MIN_MATCH_LEN, 0));
}

// ****************************************************************************
void CWorker::parse_new()
{
	v_parsing.clear();

	int data_size = (int)s_data.size();
	int ref_pred_pos = -data_size;
	int cur_lit_run_len = 0;

	const int pf_dist_l = 8;
	const int pf_dist_s = 12;

	for (int i = 0; i + MIN_MATCH_LEN < data_size;)
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
				prefetch_hts((int) v_kmers_s[i + pf_dist_s].first);

			auto h = v_kmers_s[i].first;

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
				{
					int r_len = 0;
					bool is_matching = false;
					int ref_start_pos = ref_pred_pos - cur_lit_run_len;
					int data_start_pos = i - cur_lit_run_len;
					for (int j = 0; j < cur_lit_run_len; ++j)
					{
						if (s_reference[ref_start_pos + j] == s_data[data_start_pos + j])
						{
							if (is_matching)
								r_len++;
							else
							{
								v_parsing.emplace_back(CFactor(flag_t::run_literals, 0, r_len, 0));
								r_len = 1;
								is_matching = true;
							}
						}
						else
						{
							if (is_matching)
							{
								v_parsing.emplace_back(CFactor(flag_t::match_close, ref_start_pos + j - r_len, r_len, 0));
								r_len = 1;
								is_matching = false;
							}
							else
								++r_len;
						}
					}

					if (is_matching)
						v_parsing.emplace_back(CFactor(flag_t::match_close, ref_start_pos + cur_lit_run_len - r_len, r_len, 0));
					else
						v_parsing.emplace_back(CFactor(flag_t::run_literals, 0, r_len, 0));
				}
				else
					v_parsing.emplace_back(CFactor(flag_t::run_literals, 0, cur_lit_run_len, 0));
			}

			if(abs(best_pos - ref_pred_pos) <= CLOSE_DIST)
				v_parsing.emplace_back(CFactor(flag_t::match_close, best_pos, best_len, 0));
			else
				v_parsing.emplace_back(CFactor(flag_t::match_distant, best_pos, best_len, 0));

			i += best_len;
			ref_pred_pos = best_pos + best_len;
			cur_lit_run_len = 0;
		}
		else
		{
			//			v_parsing.emplace_back(CFactor(flag_t::literal, 0, 0, s_data[i]));
			++i;
			++ref_pred_pos;
			++cur_lit_run_len;
		}

		if (cur_lit_run_len > MAX_LIT_RUN_IN_MATCH)
			ref_pred_pos = -data_size;
	}

	//	if (cur_lit_run_len)
	v_parsing.emplace_back(CFactor(flag_t::run_literals, 0, cur_lit_run_len + MIN_MATCH_LEN, 0));
}

// ****************************************************************************
void CWorker::parsing_postprocess()
{
	vector<CFactor> new_parsing;

	int lit_run_len = 0;
	int ref_pred_pos = -(int)s_data.size();
	int data_pos = 0;

	for (auto &x : v_parsing)
	{
		if (x.flag == flag_t::literal)
		{
			++lit_run_len;
			++data_pos;
		}
		else if (x.flag == flag_t::match)
		{
			if (lit_run_len)
			{
				ref_pred_pos += lit_run_len;
				new_parsing.emplace_back(CFactor(flag_t::run_literals, 0, lit_run_len, 0));
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

			new_parsing.emplace_back(x);
			ref_pred_pos += x.len;
			data_pos += x.len;
		}
		else if (x.flag == flag_t::run_literals)
		{
			new_parsing.emplace_back(x);
			if(ref_pred_pos + x.len >= (int) s_reference.size() || ref_pred_pos < 0)
				new_parsing.emplace_back(x);
			else
				for (int i = 0; i < x.len; ++i)
				{
					if (s_data[data_pos + i] == s_reference[ref_pred_pos + i])
						new_parsing.emplace_back(CFactor(flag_t::match_literal, ref_pred_pos + i, 1, 0));
					else
						new_parsing.emplace_back(CFactor(flag_t::run_literals, 0, 1, 0));
				}
			data_pos += x.len;
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
	uint32_t ht_mask = ht_size - 1u;

	hts.clear();
	hts.resize(ht_size);

	prepare_kmers(v_kmers_s, s_reference, MIN_MATCH_LEN);

	const int pf_dist = 32;
	
	for (int i = 0; i + pf_dist < (int)v_kmers_s.size(); ++i)
	{
		prefetch_hts((int) v_kmers_s[i + pf_dist].first);
		hts[v_kmers_s[i].first].emplace_back(v_kmers_s[i].second);
	}

	for (int i = max((int)v_kmers_s.size() - pf_dist, 0); i < (int)v_kmers_s.size(); ++i)
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
		else if (x.flag == flag_t::match_literal)
		{
			n_sym_in_matches += x.len;
			n_sym_in_literals += last_run_len;
			last_run_len = 0;
		}
	}

	if (mode == 1)
	{
		res.ref_size = ref_len / 2;
		res.query_size = data_len;
	}
	res.sym_in_literals[mode] = n_sym_in_literals;
	res.sym_in_matches[mode] = n_sym_in_matches;
	res.coverage[mode] = (double) (n_sym_in_literals + n_sym_in_matches) / data_len;
	res.ani[mode] = (double)n_sym_in_matches / (n_sym_in_matches + n_sym_in_literals);

	if (res.coverage[mode] < MIN_COVERAGE)
		res.ani[mode] -= 0.4;
}

// ****************************************************************************
void CWorker::calc_ani_thr(CResults &res, int mode)
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
		if (x.first + x.second >= MIN_REGION_LEN)
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
							seq.emplace_back(sym_N);
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
		if (seq[i] == sym_N)
			seq.emplace_back(sym_N);
		else if (seq[i] == sym_A)
			seq.emplace_back(sym_T);
		else if (seq[i] == sym_C)
			seq.emplace_back(sym_G);
		else if (seq[i] == sym_G)
			seq.emplace_back(sym_C);
		else if (seq[i] == sym_T)
			seq.emplace_back(sym_A);
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
//	htl.shrink_to_fit();

	hts.clear();
//	hts.shrink_to_fit();

	v_parsing.clear();
//	v_parsing.shrink_to_fit();
}

// ****************************************************************************
void CWorker::prepare_kmers(vector<pair<uint64_t, int>> &v_kmers, const seq_t &seq, int len, bool store_all)
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

		if (i >= len - 1)
		{
			if (k_len >= len)
				v_kmers.emplace_back(make_pair(k, i + 1 - len));
			else if (store_all)
				v_kmers.emplace_back(make_pair(-1, i + 1 - len));
		}
	}
}

// EOF
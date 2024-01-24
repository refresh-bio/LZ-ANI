#include "parser.h"

// ****************************************************************************
bool CParser::prepare_reference(const seq_view ref_view, uint32_t n_seqs)
{
	seq_ref.clear();

//	append(seq_ref, params.close_dist, code_N_ref);
	append(seq_ref, ref_view, code_N_ref, code_N_seq);
	append(seq_ref, params.close_dist, code_N_ref);
	append(seq_ref, params.close_dist, code_N_ref);
	append_rc(seq_ref, ref_view, code_N_ref, code_N_seq);
	append(seq_ref, params.close_dist, code_N_ref);

	n_ref_seqs = n_seqs;

	prepare_kmers(v_kmers_ref_long, seq_ref, params.min_distant_match_len, true);
	prepare_kmers(v_kmers_ref_short, seq_ref, params.min_match_len, true);
	prepare_ht_long();
	prepare_ht_short();

	return true;
}

// ****************************************************************************
bool CParser::prepare_data(const seq_view data_view, uint32_t n_seqs)
{
	seq_data.clear();

//	append(seq_data, params.close_dist, code_N_seq);
	append(seq_data, data_view, code_N_seq, code_N_ref);
	append(seq_data, params.close_dist, code_N_seq);

	n_data_seqs = n_seqs;

	prepare_kmers(v_kmers_data_short, seq_data, params.min_match_len, true);
	prepare_kmers(v_kmers_data_long, seq_data, params.min_distant_match_len, true);

	return true;
}

// ****************************************************************************
void CParser::prepare_kmers(vector<pair<int64_t, int64_t>>& v_kmers, const seq_t& seq, int len, bool store_all)
{
	v_kmers.clear();
	v_kmers.resize(seq.size());

	uint64_t k = 0u;
	uint64_t mask = (~0ull) >> (64 - 2 * len);
	int k_len = 0;
	int seq_size = (int)seq.size();

	int i;
	int i_kmer = 0;

	for (i = 0; i < seq_size && i + 1 < len; ++i)
	{
		if (seq[i] >= code_N)
		{
			k = 0ull;
			k_len = 0;
		}
		else
		{
			k <<= 2;
			++k_len;
			k += seq[i];
			k &= mask;
		}
	}

	for (; i < seq_size; ++i, ++i_kmer)
	{
		uint64_t c = seq[i];

		k <<= 2;
		++k_len;
		k += c;
		k &= mask;

		if (c >= code_N)
			k_len = 0;

		v_kmers[i_kmer].first = (k_len >= len) ? (int64_t)k : -1;
		v_kmers[i_kmer].second = i + 1 - len;
	}

	if (store_all)
		for (int i = 0; i < len - 1; ++i)
			v_kmers[i_kmer++] = make_pair(-1, 0);

	v_kmers.resize(i_kmer);
}

// ****************************************************************************
void CParser::prepare_ht_short()
{
	uint32_t ht_size = 1u << (2 * params.min_match_len);

	ht_short_desc.clear();
	ht_short_desc.resize(ht_size, make_pair(0, 0));
	ht_short.resize(v_kmers_ref_short.size());

	const int pf_dist = 16;
	int n_kmers = (int)v_kmers_ref_short.size();

	for (int i = 0; i < n_kmers; ++i)
	{
		if (i + pf_dist < n_kmers && v_kmers_ref_short[i + pf_dist].first >= 0)
			prefetch_hts1((int)v_kmers_ref_short[i + pf_dist].first);

		if (v_kmers_ref_short[i].first >= 0)
			++ht_short_desc[v_kmers_ref_short[i].first].second;
	}

	for (int i = 1; i < (int)ht_size; ++i)
		ht_short_desc[i].first = ht_short_desc[i - 1].first + ht_short_desc[i - 1].second;

	for (int i = 0; i < n_kmers; ++i)
	{
		if (i + pf_dist < n_kmers && v_kmers_ref_short[i + pf_dist].first >= 0)
			prefetch_hts1((int)v_kmers_ref_short[i + pf_dist].first);

		if (v_kmers_ref_short[i].first >= 0)
		{
			ht_short[ht_short_desc[v_kmers_ref_short[i].first].first].first = v_kmers_ref_short[i].second;

			if (v_kmers_ref_long[i].first >= 0)
				ht_short[ht_short_desc[v_kmers_ref_short[i].first].first].second = ((int)v_kmers_ref_long[i].first) & ht_short_mask;
			else
				ht_short[ht_short_desc[v_kmers_ref_short[i].first].first].second = -1;

			++ht_short_desc[v_kmers_ref_short[i].first].first;
		}
	}

	for (int i = 0; i < (int)ht_size; ++i)
		ht_short_desc[i].first -= ht_short_desc[i].second;
}

// ****************************************************************************
void CParser::prepare_ht_long()
{
	uint32_t x = (uint32_t)(v_kmers_ref_long.size() / ht_long_max_fill_factor);

	while (x & (x - 1))
		x &= x - 1;

	ht_long_size = 2 * x;
	ht_long_mask = ht_long_size - 1;

	ht_long.clear();
	ht_long.resize(ht_long_size, HT_EMPTY);

	const int pf_dist = 16;

	for (int i = 0; i + pf_dist < (int)v_kmers_ref_long.size(); ++i)
	{
		if (v_kmers_ref_long[i + pf_dist].first >= 0)
			prefetch_htl(hash_mm(v_kmers_ref_long[i + pf_dist].first) & ht_long_mask);

		if (v_kmers_ref_long[i].first < 0)
			continue;

		auto ht_idx = hash_mm(v_kmers_ref_long[i].first) & ht_long_mask;

		while (ht_long[ht_idx] != HT_EMPTY)
			ht_idx = (ht_idx + 1) & ht_long_mask;

		ht_long[ht_idx] = (int)v_kmers_ref_long[i].second;
	}

	for (int i = max((int)v_kmers_ref_long.size() - pf_dist, 0); i < (int)v_kmers_ref_long.size(); ++i)
	{
		if (v_kmers_ref_long[i].first < 0)
			continue;

		auto ht_idx = hash_mm(v_kmers_ref_long[i].first) & ht_long_mask;

		while (ht_long[ht_idx] != HT_EMPTY)
			ht_idx = (ht_idx + 1) & ht_long_mask;

		ht_long[ht_idx] = (int)v_kmers_ref_long[i].second;
	}
}

// ****************************************************************************
int CParser::equal_len(int ref_pos, int data_pos, int starting_pos)
{
	int r;
	int ref_len = (int)seq_ref.size();
	int data_len = (int)seq_data.size();
	int max_r = min(ref_len - ref_pos, data_len - data_pos);

	uint8_t* p_ref = seq_ref.data() + ref_pos + starting_pos;
	uint8_t* p_data = seq_data.data() + data_pos + starting_pos;

	for (r = starting_pos; r < max_r; ++r)
		if (*p_ref++ != *p_data++)
			break;

	return r;
}

// ****************************************************************************
int CParser::est_equal_len(int64_t x, int64_t y)
{
	if (x < 0 || y < 0)
		return params.min_distant_match_len;

	//	return MIN_MATCH_LEN + lzcnt32((uint32_t) ((int) x & hts_mask) ^ (uint32_t)(y)) / 2 - (16 - (MIN_DISTANT_MATCH_LEN - MIN_MATCH_LEN));
	return est_len_correction + lzcnt32((uint32_t)((int)x & ht_short_mask) ^ (uint32_t)(y)) / 2;
}

// ****************************************************************************
void CParser::compare_ranges(int data_start_pos, int ref_start_pos, int len, bool backward = false)
{
	int r_len = 0;
	bool is_matching = false;
	flag_t flag = backward ? flag_t::match_distant : flag_t::match_close;

	for (int j = 0; j < len; ++j)
	{
		if (seq_ref[ref_start_pos + j] == seq_data[data_start_pos + j])
		{
			if (is_matching)
				r_len++;
			else
			{
				if (r_len)
					v_parsing.emplace_back(factor_t(data_start_pos + j - r_len, flag_t::run_literals, 0, r_len, 0));
				r_len = 1;
				is_matching = true;
			}
		}
		else
		{
			if (is_matching)
			{
				v_parsing.emplace_back(factor_t(data_start_pos + j - r_len, flag, ref_start_pos + j - r_len, r_len, 0));
				r_len = 1;
				is_matching = false;
				flag = flag_t::match_close;
			}
			else
				++r_len;
		}
	}

	if (is_matching)
		v_parsing.emplace_back(factor_t(data_start_pos + len - r_len, flag, ref_start_pos + len - r_len, r_len, 0));
	else if (r_len)
		v_parsing.emplace_back(factor_t(data_start_pos + len - r_len, flag_t::run_literals, 0, r_len, 0));
}

// ****************************************************************************
/*int CParser::try_extend_forward(int data_start_pos, int ref_start_pos)
{
	int data_size = (int)seq_data.size();
	int ref_size = (int)seq_ref.size();

	int approx_ext;
	int no_missmatches = 0;
	int last_match = 0;
	vector<int> window(params.approx_window, 0);

	for (approx_ext = 0; data_start_pos + approx_ext < data_size && ref_start_pos + approx_ext < ref_size; ++approx_ext)
	{
		bool is_missmatch = seq_data[data_start_pos + approx_ext] != seq_ref[ref_start_pos + approx_ext];
		no_missmatches -= window[approx_ext % params.approx_window];
		window[approx_ext % params.approx_window] = is_missmatch;
		no_missmatches += is_missmatch;

		if (!is_missmatch)
			last_match = approx_ext + 1;

		if (no_missmatches > params.approx_mismatches)
			break;
	}

	return last_match;
}*/

// ****************************************************************************
int CParser::try_extend_forward2(int data_start_pos, int ref_start_pos)
{
	int data_size = (int)seq_data.size();
	int ref_size = (int)seq_ref.size();

	int approx_ext;
	int no_missmatches = 0;
	int last_run_match = 0;
	vector<int> window(params.approx_window, 0);
	int match_run_len = params.approx_run_len;

	for (approx_ext = 0; data_start_pos + approx_ext < data_size && ref_start_pos + approx_ext < ref_size; ++approx_ext)
	{
		bool is_missmatch = seq_data[data_start_pos + approx_ext] != seq_ref[ref_start_pos + approx_ext];
		no_missmatches -= window[approx_ext % params.approx_window];
		window[approx_ext % params.approx_window] = is_missmatch;
		no_missmatches += is_missmatch;

		if (!is_missmatch)
		{
			if (++match_run_len >= params.approx_run_len)
				last_run_match = approx_ext + 1;
		}
		else
			match_run_len = 0;

		if (no_missmatches > params.approx_mismatches)
			break;
	}

	return last_run_match;
}

// ****************************************************************************
/*int CParser::try_extend_backward(int data_start_pos, int ref_start_pos, int max_len)
{
	int approx_ext;
	int no_missmatches = 0;
	int last_match = 0;
	vector<int> window(params.approx_window, 0);

	for (approx_ext = 0; data_start_pos - approx_ext > 0 && ref_start_pos - approx_ext > 0 && approx_ext < max_len; ++approx_ext)
	{
		bool is_missmatch = seq_data[data_start_pos - approx_ext - 1] != seq_ref[ref_start_pos - approx_ext - 1];
		no_missmatches -= window[approx_ext % params.approx_window];
		window[approx_ext % params.approx_window] = is_missmatch;
		no_missmatches += is_missmatch;

		if (!is_missmatch)
			last_match = approx_ext + 1;

		if (no_missmatches > params.approx_mismatches)
			break;
	}

	return last_match;
}*/

// ****************************************************************************
int CParser::try_extend_backward2(int data_start_pos, int ref_start_pos, int max_len)
{
	int approx_ext;
	int no_missmatches = 0;
	int last_run_match = 0;
	vector<int> window(params.approx_window, 0);
	int match_run_len = params.approx_run_len;

	for (approx_ext = 0; data_start_pos - approx_ext > 0 && ref_start_pos - approx_ext > 0 && approx_ext < max_len; ++approx_ext)
	{
		bool is_missmatch = seq_data[data_start_pos - approx_ext - 1] != seq_ref[ref_start_pos - approx_ext - 1];
		no_missmatches -= window[approx_ext % params.approx_window];
		window[approx_ext % params.approx_window] = is_missmatch;
		no_missmatches += is_missmatch;

		if (!is_missmatch)
		{
			if (++match_run_len >= params.approx_run_len)
				last_run_match = approx_ext + 1;
		}
		else
			match_run_len = 0;

		if (no_missmatches > params.approx_mismatches)
			break;
	}

	return last_run_match;
}

// ****************************************************************************
void CParser::parse()
{
	v_parsing.clear();

	int data_size = (int)seq_data.size();
	int ref_pred_pos = -data_size;
	int cur_lit_run_len = 0;

	const int pf_dist_l = 8;
	const int pf_dist_s1 = 24;
	const int pf_dist_s2 = 12;

	int i;

	int prev_region_start = -1;
	int prev_region_end = 0;

	for (i = 0; i + params.min_match_len < data_size;)
	{
		int best_pos = 0;
		int best_len = 0;
		int h;

		if (ref_pred_pos < 0)
		{
			// Look for long match
			if (i + pf_dist_l < data_size && v_kmers_data_long[i + pf_dist_l].first >= 0)
				if (v_kmers_data_long[i + pf_dist_l].first >= 0)
					prefetch_htl(hash_mm(v_kmers_data_long[i + pf_dist_l].first) & ht_long_mask);

			if (v_kmers_data_long[i].first >= 0)
			{
				h = hash_mm(v_kmers_data_long[i].first) & ht_long_mask;

				for (; ht_long[h] != HT_EMPTY; h = (h + 1) & ht_long_mask)
				{
					int matching_len = equal_len(ht_long[h], i);

					if (matching_len < params.min_distant_match_len)
						continue;

					if (matching_len > best_len)
					{
						best_len = matching_len;
						best_pos = ht_long[h];
					}
				}
			}
		}
		else
		{
			// Look for short but close match
			if (i + pf_dist_s1 < data_size && v_kmers_data_short[i + pf_dist_s1].first >= 0)
				prefetch_hts1((int)v_kmers_data_short[i + pf_dist_s1].first);
			if (i + pf_dist_s2 < data_size && v_kmers_data_short[i + pf_dist_s2].first >= 0)
				prefetch_hts2((int)v_kmers_data_short[i + pf_dist_s2].first);

			auto h = v_kmers_data_short[i].first;

			if (h != HT_FAIL)
			{
				int bucket_size = ht_short_desc[h].second;
				auto* bucket = ht_short.data() + ht_short_desc[h].first;

				for (int j = 0; j < bucket_size; ++j)
				{
					auto pos = bucket[j].first;
//#if 0
					int est_matching_len = est_equal_len(v_kmers_data_long[i].first, bucket[j].second);

					int matching_len;

					if (est_matching_len >= params.min_distant_match_len)
						matching_len = equal_len(pos, i, params.min_match_len);
					else
						matching_len = est_matching_len;
//#endif

//					int matching_len = equal_len(pos, i, params.min_match_len);
					if (matching_len < params.min_distant_match_len)			
					{
						int dist = pos - ref_pred_pos;
						if (dist > params.close_dist || dist <= -matching_len)
							continue;
					}

					if (matching_len > best_len)
					{
						best_len = matching_len;
						best_pos = pos;
					}
				}
			}
		}

		if (best_len >= params.min_match_len)
		{
			if (cur_lit_run_len)
			{
				if (ref_pred_pos >= 0)
					compare_ranges(i - cur_lit_run_len, ref_pred_pos - cur_lit_run_len, cur_lit_run_len);
				else
					v_parsing.emplace_back(factor_t(i - cur_lit_run_len, flag_t::run_literals, 0, cur_lit_run_len, 0));
			}

			flag_t flag = flag_t::match_distant;

			if (abs(best_pos - ref_pred_pos) <= params.close_dist)
			{
				v_parsing.emplace_back(factor_t(i, flag_t::match_close, best_pos, best_len, 0));
			}
			else
			{
				// Remove previous region if too short
				if (prev_region_start >= 0 && prev_region_end - prev_region_start < params.min_region_len)
				{
					while (!v_parsing.empty() && v_parsing.back().data_pos >= prev_region_start)
						v_parsing.pop_back();
					int run_len = i - prev_region_start;

					while (!v_parsing.empty() && v_parsing.back().flag == flag_t::run_literals)
					{
						run_len += v_parsing.back().len;
						v_parsing.pop_back();
					}

					v_parsing.emplace_back(factor_t(i - run_len, flag_t::run_literals, 0, run_len, 0));
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

				v_parsing.emplace_back(factor_t(i, flag, best_pos, best_len, 0));
				if (flag == flag_t::match_distant)
					prev_region_start = i;

				if (prev_region_start < 0)
					for (int j = (int)v_parsing.size() - 1; j >= 0; --j)
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

		if (cur_lit_run_len > params.max_lit_run_in_match)
			ref_pred_pos = -data_size;
	}

	if (ref_pred_pos < 0)
		v_parsing.emplace_back(factor_t(i - cur_lit_run_len, flag_t::run_literals, 0, cur_lit_run_len + (data_size - i), 0));
	else
		compare_ranges(i - cur_lit_run_len, ref_pred_pos - cur_lit_run_len - params.min_match_len, cur_lit_run_len + (data_size - i));
}

// ****************************************************************************
results_t CParser::calc_stats()
{
	vector<pair<int, int>> v_matches;
	int cur_match_len = 0;
	int cur_match_lit = 0;
	int n_lit = 0;

	for (const auto& x : v_parsing)
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

//	int ref_len = (int)seq_ref.size() - n_ref_seqs * params.close_dist;
//	int data_len = (int)seq_data.size() - n_data_seqs * params.close_dist;
	int n_sym_in_matches = 0;
	int n_sym_in_literals = 0;

	int n_components = 0;

	for (auto x : v_matches)
		if (x.first + x.second >= params.min_region_len) // && (double) x.first / (x.first + x.second) > 0.5
		{
			n_sym_in_matches += x.first;
			n_sym_in_literals += x.second;
			++n_components;
		}

	return results_t(n_sym_in_matches, n_sym_in_literals, n_components);
}


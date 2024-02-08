#include "parser.h"

#include <iostream>

// ****************************************************************************
bool CParser::prepare_reference(const seq_view ref_view, uint32_t n_seqs)
{
	seq_ref.clear();

	append(seq_ref, ref_view, code_N_ref, code_N_seq);
	append(seq_ref, params.close_dist, code_N_ref);
	append(seq_ref, params.close_dist, code_N_ref);
	append_rc(seq_ref, ref_view, code_N_ref, code_N_seq);
	append(seq_ref, params.close_dist, code_N_ref);

	n_ref_seqs = n_seqs;

	prepare_kmers(v_kmers_ref_long, seq_ref, params.min_anchor_len, true);
	prepare_kmers(v_kmers_ref_short, seq_ref, params.min_match_len, true);
	prepare_ht_long();
	prepare_ht_short();

	return true;
}

// ****************************************************************************
bool CParser::prepare_data(const seq_view data_view, uint32_t n_seqs)
{
	seq_data.clear();

	append(seq_data, data_view, code_N_seq, code_N_ref);
	append(seq_data, params.close_dist, code_N_seq);

	n_data_seqs = n_seqs;

	prepare_kmers(v_kmers_data_short, seq_data, params.min_match_len, true);
	prepare_kmers(v_kmers_data_long, seq_data, params.min_anchor_len, true);

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
			ht_short[ht_short_desc[v_kmers_ref_short[i].first].first] = v_kmers_ref_short[i].second;
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
int CParser::equal_len(const int ref_pos, const int data_pos, const int starting_pos)  const
{
	int r;
	int ref_len = (int)seq_ref.size();
	int data_len = (int)seq_data.size();
	int max_r = min(ref_len - ref_pos, data_len - data_pos);

	const uint8_t* p_ref = seq_ref.data() + ref_pos + starting_pos;
	const uint8_t* p_data = seq_data.data() + data_pos + starting_pos;

	for (r = starting_pos; r < max_r; ++r)
		if (*p_ref++ != *p_data++)
			break;

	return r;
}

// ****************************************************************************
void CParser::compare_ranges(const int data_start_pos, const int ref_start_pos, const int len, const bool backward = false)
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
					v_parsing.emplace_back(data_start_pos + j - r_len, flag_t::run_literals, 0, r_len);
				r_len = 1;
				is_matching = true;
			}
		}
		else
		{
			if (is_matching)
			{
				v_parsing.emplace_back(data_start_pos + j - r_len, flag, ref_start_pos + j - r_len, r_len);
				r_len = 1;									
				is_matching = false;
				flag = flag_t::match_close;
			}
			else
				++r_len;
		}
	}

	if (is_matching)
		v_parsing.emplace_back(data_start_pos + len - r_len, flag, ref_start_pos + len - r_len, r_len);
	else if (r_len)
		v_parsing.emplace_back(data_start_pos + len - r_len, flag_t::run_literals, 0, r_len);
}

// ****************************************************************************
void CParser::compare_ranges_both_ways(const int data_start_pos, const int ref_start_pos_left, const int ref_end_pos_right, const int len, 
	vector<pair<int, bool>>& left_side, vector<pair<int, bool>>& right_side)
{
	left_side.clear();
	right_side.clear();

	int to_scan;
	
	if (ref_end_pos_right < ref_start_pos_left)
		to_scan = len;
	else
		to_scan = min(ref_end_pos_right - ref_start_pos_left, len);

#if _DEBUG
	string data, ref_left, ref_right;

	for (int i = 0; i < len; ++i)
	{
		data.push_back("ACGTNN"[seq_data[data_start_pos + i]]);
	}

	for (int i = 0; i < to_scan; ++i)
	{
		ref_left.push_back("ACGTNN"[seq_ref[ref_start_pos_left + i]]);
		ref_right.push_back("ACGTNN"[seq_ref[ref_end_pos_right - len + i]]);
	}

	vector<factor_t> loc_parsing;
#else
	vector<factor_t> &loc_parsing = v_parsing;
#endif

	int no_matches = 0;

	left_side.emplace_back(0, false);
	for (int i = 0; i < to_scan; ++i)
	{
		bool is_match = seq_ref[ref_start_pos_left + i] == seq_data[data_start_pos + i];
		left_side.emplace_back(no_matches += (int) is_match, is_match);
	}

	no_matches = 0;

	right_side.emplace_back(0, false);
	for (int i = 1; i <= min(to_scan, ref_end_pos_right); ++i)
	{
		bool is_match = seq_ref[ref_end_pos_right - i] == seq_data[data_start_pos + len - i];
		right_side.emplace_back(no_matches += (int) is_match, is_match);
	}
	right_side.resize(to_scan + 1, make_pair(0, false));

	int max_no_matches = 0;
	int best_split = 0;

	for (int i = 0; i <= to_scan; ++i)
	{
		no_matches = left_side[i].first + right_side[to_scan - i].first;
		if (no_matches >= max_no_matches)
		{
			max_no_matches = no_matches;
			best_split = i;
		}
	}

	flag_t item_flags[2] = { flag_t::run_literals, flag_t::match_close };

	// Store left
	if (best_split > 0)
	{
		bool is_match = left_side[1].second;
		auto cf = item_flags[is_match];
		int data_p = data_start_pos;
		loc_parsing.emplace_back(data_p++, cf, left_side[1].second ? ref_start_pos_left : 0, 1);

		for (int i = 2; i <= best_split; ++i, ++data_p)
		{
			is_match = left_side[i].second;
			cf = item_flags[is_match];
			if (cf == loc_parsing.back().flag)
				loc_parsing.back().len++;
			else
				loc_parsing.emplace_back(data_p, cf, is_match ? ref_start_pos_left + i - 1 : 0, 1);
		}
	}

	// Store middle
	if (to_scan < len)
	{
		if(best_split > 0 && loc_parsing.back().flag == flag_t::run_literals)
			loc_parsing.back().len += len - to_scan;
		else
			loc_parsing.emplace_back(data_start_pos + best_split, flag_t::run_literals, 0, len - to_scan);
	}

	// Store right
	if (best_split < to_scan)
	{
		int shift = len - to_scan;
		int from_right = to_scan - best_split;

		bool is_match = right_side[from_right].second;
		auto cf = item_flags[is_match];
		int data_p = data_start_pos + best_split + shift;

		if(!is_match && (best_split > 0 || shift > 0) && loc_parsing.back().flag == flag_t::run_literals)
			loc_parsing.back().len++;
		else
			loc_parsing.emplace_back(data_p++, cf, is_match ? ref_end_pos_right - from_right : 0, 1);

		for (int i = from_right - 1; i > 0; --i, ++data_p)
		{
			is_match = right_side[i].second;
			cf = item_flags[is_match];
			if (cf == loc_parsing.back().flag)
				loc_parsing.back().len++;
			else
				loc_parsing.emplace_back(data_p, cf, is_match ? ref_end_pos_right - i : 0, 1);
		}
	}

#if _DEBUG
	v_parsing.insert(v_parsing.end(), loc_parsing.begin(), loc_parsing.end());
#endif

//	compare_ranges(data_start_pos, ref_start_pos_left, len);
}

// ****************************************************************************
int CParser::try_extend_forward(const int data_start_pos, const int ref_start_pos, vector<int> &window)
{
	int data_size = (int)seq_data.size();
	int ref_size = (int)seq_ref.size();

	int approx_ext;
	int no_missmatches = 0;
	int last_run_match = 0;
	window.clear();
	window.resize(params.approx_window, 0);
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
int CParser::try_extend_backward(const int data_start_pos, const int ref_start_pos, const int max_len, vector<int> &window)
{
	int approx_ext;
	int no_missmatches = 0;
	int last_run_match = 0;
	window.clear();
	window.resize(params.approx_window, 0);
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
// Check if the region should be preserved
// True means preserve
bool CParser::eval_region(int region_start, int region_end)
{
	// Simple test: just check region length
	return region_end - region_start >= params.min_region_len;

/*	int no_miss = 0;

	if (v_parsing.empty())
		return false;

	auto p = v_parsing.rbegin();
	if (p->flag == flag_t::run_literals)
		++p;

	for (; p != v_parsing.rend() && p->data_pos >= region_start; ++p)
		if (p->flag == flag_t::run_literals)
			no_miss += p->len;

	double region_len = region_end - region_start;

	if (region_len == 0)
		return true;					// Never should be here

	if ((region_len - no_miss) / region_len < 0.9)
		return region_len >= params.min_region_len;

	// Preserve shorter but high-quality regions
	return region_len >= params.min_region_len / 2;*/
}

// ****************************************************************************
void CParser::parse()
{
	v_parsing.clear();

	int data_size = (int)seq_data.size();
	int ref_pred_pos = -data_size;
	int cur_lit_run_len = 0;

	const int pf_dist_l = 8;
//	const int pf_dist_s1 = 24;
//	const int pf_dist_s2 = 12;

	int i;

	int prev_region_start = -1;
	int prev_region_end = 0;

	vector<pair<int, bool>> left_side;
	vector<pair<int, bool>> right_side;
	vector<int> window;

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

					if (matching_len < params.min_anchor_len)
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
			// Look for long match
/*			if (i + pf_dist_l < data_size && v_kmers_data_long[i + pf_dist_l].first >= 0)
				if (v_kmers_data_long[i + pf_dist_l].first >= 0)
					prefetch_htl(hash_mm(v_kmers_data_long[i + pf_dist_l].first) & ht_long_mask);
					*/
			prefetch_htl(hash_mm(v_kmers_data_long[i].first) & ht_long_mask);

			// Look for short but close match
/*				if (i + pf_dist_s1 < data_size && v_kmers_data_short[i + pf_dist_s1].first >= 0)
				prefetch_hts1((int)v_kmers_data_short[i + pf_dist_s1].first);
			if (i + pf_dist_s2 < data_size && v_kmers_data_short[i + pf_dist_s2].first >= 0)
				prefetch_hts2((int)v_kmers_data_short[i + pf_dist_s2].first);
				*/
			auto h = v_kmers_data_short[i].first;

			if (h != HT_FAIL)
			{
				int bucket_size = ht_short_desc[h].second;
				auto* bucket = ht_short.data() + ht_short_desc[h].first;

				int j_start = lower_bound(bucket, bucket + bucket_size, ref_pred_pos - cur_lit_run_len) - bucket;

				for (int j = j_start; j < bucket_size && bucket[j] < ref_pred_pos + params.close_dist; ++j)
				{
					auto pos = bucket[j];
					int matching_len = equal_len(pos, i, params.min_match_len);

//					if (matching_len == params.min_match_len && abs(pos - ref_pred_pos) >= params.close_dist / 16)
					if (matching_len == params.min_match_len && cur_lit_run_len > params.max_lit_run_in_match / 8)
						matching_len = 0;

					if (matching_len >= best_len)
					{
						if (matching_len == best_len)
						{
							if (abs(pos - ref_pred_pos) < abs(best_pos - ref_pred_pos))
								best_pos = pos;
						}
						else
						{
							best_len = matching_len;
							best_pos = pos;
						}
					}
				}
			}

			int best_anchor_len = 0;
			int best_anchor_pos = 0;

			if (v_kmers_data_long[i].first >= 0)
			{
				h = hash_mm(v_kmers_data_long[i].first) & ht_long_mask;

				for (; ht_long[h] != HT_EMPTY; h = (h + 1) & ht_long_mask)
				{
					int matching_len = equal_len(ht_long[h], i);

					if (matching_len < params.min_anchor_len)
						continue;

					if (matching_len > best_anchor_len)
					{
						best_anchor_len = matching_len;
						best_anchor_pos = ht_long[h];
					}
				}
			}

			if (best_anchor_pos)
			{
				if (!best_pos)
				{
					best_pos = best_anchor_pos;
					best_len = best_anchor_len;
				}
				else
				{
					// Approximate probabilities that match is by a chance
					double anchor_prob = ipow(1 - prob_len(best_anchor_len), 2 * (seq_ref.size() + 1 - best_anchor_len));
					double close_prob = ipow(1 - prob_len(best_len), cur_lit_run_len + params.close_dist + 1 - best_len);

					if (anchor_prob > close_prob)
					{
						best_pos = best_anchor_pos;
						best_len = best_anchor_len;
					}
				}
			}
		}

		if (best_len >= params.min_match_len)
		{
			flag_t flag = flag_t::match_distant;

			if (ref_pred_pos >= 0 && abs(best_pos - ref_pred_pos) <= params.close_dist)
			{
				compare_ranges_both_ways(i - cur_lit_run_len, ref_pred_pos - cur_lit_run_len, best_pos + best_len, cur_lit_run_len, left_side, right_side);

				v_parsing.emplace_back(i, flag_t::match_close, best_pos, best_len);
			}
			else
			{
				if (cur_lit_run_len)
					v_parsing.emplace_back(i - cur_lit_run_len, flag_t::run_literals, 0, cur_lit_run_len);

				// Remove previous region if too short
//				if (prev_region_start >= 0 && prev_region_end - prev_region_start < params.min_region_len)
				if (prev_region_start >= 0 && !eval_region(prev_region_start, prev_region_end))
				{
					while (!v_parsing.empty() && v_parsing.back().data_pos >= prev_region_start)
						v_parsing.pop_back();
					int run_len = i - prev_region_start;

					while (!v_parsing.empty() && v_parsing.back().flag == flag_t::run_literals)
					{
						run_len += v_parsing.back().len;
						v_parsing.pop_back();
					}

					v_parsing.emplace_back(i - run_len, flag_t::run_literals, 0, run_len);
					prev_region_start = -1;
				}

				if (!v_parsing.empty() && v_parsing.back().flag == flag_t::run_literals)
				{
					int approx_pred = try_extend_backward(i, best_pos, v_parsing.back().len, window);
					if (approx_pred)
					{
						v_parsing.back().len -= approx_pred;
						if (v_parsing.back().len == 0)
							v_parsing.pop_back();
						compare_ranges(i - approx_pred, best_pos - approx_pred, approx_pred, true);
						flag = flag_t::match_close;
//						prev_region_start = -1;
						prev_region_start = i - approx_pred;
					}
				}

				v_parsing.emplace_back(i, flag, best_pos, best_len);
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

			int approx_ext = try_extend_forward(i, ref_pred_pos, window);
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
		v_parsing.emplace_back(i - cur_lit_run_len, flag_t::run_literals, 0, cur_lit_run_len + (data_size - i));
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

// EOF

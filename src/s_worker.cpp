#include "s_worker.h"
#include <iostream>
#include <iomanip>
#include "xmmintrin.h"
#include <nmmintrin.h>
#include <immintrin.h>
//#include <intrin.h>
#include <algorithm>
#include <utility>
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


// ****************************************************************************
CSharedWorker::CSharedWorker() : BaseWorker()
{
	
	s_reference = nullptr;
	htl = nullptr;
	hts = nullptr;
	hts2 = nullptr;
	hts3 = nullptr;
	hts3_desc = nullptr;

	est_len_correction = MIN_MATCH_LEN - (16 - (MIN_DISTANT_MATCH_LEN - MIN_MATCH_LEN));
}

// ****************************************************************************
CSharedWorker::~CSharedWorker()
{
	if (s_reference)
		delete s_reference;

	if(htl)			delete htl;
	if(hts)			delete hts;
	if(hts2)		delete hts2;
	if(hts3)		delete hts3;
	if(hts3_desc)	delete hts3_desc;
}

// ****************************************************************************
int CSharedWorker::equal_len(int ref_pos, int data_pos, int starting_pos)
{
	int r;
	int ref_len = (int)s_reference->size();
	int data_len = (int)s_data.size();
	int max_r = min(ref_len - ref_pos, data_len - data_pos);

	uint8_t *p_ref = s_reference->data() + ref_pos + starting_pos;
	uint8_t *p_data = s_data.data() + data_pos + starting_pos;

	for (r = starting_pos; r < max_r; ++r)
		if (*p_ref++ != *p_data++)
			break;

	return r;
}

// ****************************************************************************
int CSharedWorker::est_equal_len(int64_t x, int64_t y)
{
	if (x < 0 || y < 0)
		return MIN_DISTANT_MATCH_LEN;
	
//	return MIN_MATCH_LEN + lzcnt32((uint32_t) ((int) x & hts_mask) ^ (uint32_t)(y)) / 2 - (16 - (MIN_DISTANT_MATCH_LEN - MIN_MATCH_LEN));
	return est_len_correction + lzcnt32((uint32_t)((int)x & hts_mask) ^ (uint32_t)(y)) / 2;
}

// ****************************************************************************
bool CSharedWorker::load_reference(string fn_ref, Genome* buffered_data)
{
	if (s_reference)
		delete s_reference;
	s_reference = new seq_t;

	if (buffered_data)
	{
		if (buffered_data->n_seqs() == 0)
		{
			if (!load_file(fn_ref, *buffered_data, sym_N1))
			{
				cerr << "Error: Cannot load " + fn_ref + "\n";
				return false;
			}
		}

		*s_reference = buffered_data->seq;
		n_reference = buffered_data->n_seqs();
	}
	else 
	{
		Genome genome;
		
		if (!load_file(fn_ref, genome, sym_N1))
		{
			cerr << "Error: Cannot load " + fn_ref + "\n";
			return false;
		}

		*s_reference = genome.seq;
		n_reference = genome.n_seqs();
	}

	duplicate_rev_comp(*s_reference);

	return true;
}

// ****************************************************************************
bool CSharedWorker::load_data(string fn_data, Genome *buffered_data)
{
	if (buffered_data && buffered_data->n_seqs() > 0)
	{
		s_data = buffered_data->seq;
		n_data = buffered_data->n_seqs();
		
		return true;
	}

	Genome genome;

	if (!load_file(fn_data, genome, sym_N2))
	{
		cerr << "Error: Cannot load " + fn_data + "\n";
		return false;
	}

	s_data = genome.seq;
	n_data = genome.n_seqs();

	if (buffered_data)
	{
		*buffered_data = genome;
	}

	return true;
}

// ****************************************************************************
bool CSharedWorker::share_from(CSharedWorker* base)
{
	if (base == nullptr)
	{
		s_reference = nullptr;
		n_reference = 0;
		htl = nullptr;
		hts = nullptr;
		hts2 = nullptr;
		hts3 = nullptr;
		hts3_desc = nullptr;
	}
	else
	{
		s_reference = base->s_reference;
		n_reference = base->n_reference;
		htl = base->htl;
		hts = base->hts;
		hts2 = base->hts2;
		hts3 = base->hts3;
		hts3_desc = base->hts3_desc;

		htl_size = base->htl_size;
		htl_mask = base->htl_mask;
		hts_mask = base->hts_mask;
	}

	return true;
}

// ****************************************************************************
void CSharedWorker::parse()
{
	v_parsing.clear();

	int data_size = (int)s_data.size();
	int ref_pred_pos = -data_size;
	int cur_lit_run_len = 0;

	const int pf_dist_l = 8;
	const int pf_dist_s1 = 24;
	const int pf_dist_s2 = 12;

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
			if (i + pf_dist_l < data_size && v_kmers_dl[i + pf_dist_l].first >= 0)
				if(v_kmers_dl[i + pf_dist_l].first >= 0)
					prefetch_htl(hash_mm(v_kmers_dl[i + pf_dist_l].first, htl_mask));

			if (v_kmers_dl[i].first >= 0)
			{
				h = hash_mm(v_kmers_dl[i].first, htl_mask);

				for (; (*htl)[h] != HT_EMPTY; h = (h + 1) & htl_mask)
				{
					int matching_len = equal_len((*htl)[h], i);

					if (matching_len < MIN_DISTANT_MATCH_LEN)
						continue;

					if (matching_len > best_len)
					{
						best_len = matching_len;
						best_pos = (*htl)[h];
					}
				}
			}
		}
		else
		{
			// Look for short but close match
			if (i + pf_dist_s1 < data_size && v_kmers_ds[i + pf_dist_s1].first >= 0)
				prefetch_hts1((int)v_kmers_ds[i + pf_dist_s1].first);
			if (i + pf_dist_s2 < data_size && v_kmers_ds[i + pf_dist_s2].first >= 0)
				prefetch_hts2((int)v_kmers_ds[i + pf_dist_s2].first);

			auto h = v_kmers_ds[i].first;

			if (h != HT_FAIL)
			{
//				int bucket_size = (int)hts[h].size();
//				int bucket_size = (int)hts2[h].size();
//				auto &bucket = hts[h];
//				auto& bucket = hts2[h];
				int bucket_size = (*hts3_desc)[h].second;
				auto* bucket = hts3->data() + (*hts3_desc)[h].first;

//				const int pf_dist = 4;

/*				for (int j = 0; j < min(pf_dist, bucket_size); ++j)
//					prefetch(bucket[j]);
					prefetch(bucket[j].first);
					*/
				int best_close_len = 0;
				int best_close_pos = 0;

				for (int j = 0; j < bucket_size; ++j)
				{
/*					if (j + pf_dist < bucket_size)
//						prefetch(bucket[j + pf_dist]);
						prefetch(bucket[j + pf_dist].first);*/

//					auto pos = bucket[j];
					auto pos = bucket[j].first;
					int est_matching_len = est_equal_len(v_kmers_dl[i].first, bucket[j].second);

					/*if (est_matching_len < best_len)
						continue;*/
					
					int matching_len;
					
/*					if (est_matching_len != equal_len(pos, i, MIN_MATCH_LEN) && est_matching_len != 32 && est_matching_len != MIN_DISTANT_MATCH_LEN)
					{
						int aa = equal_len(pos, i, MIN_MATCH_LEN);

						cout << est_matching_len << ":" << aa << endl;
						cout << hex << v_kmers_dl[i].first << " : " << bucket[j].second << dec << endl;
					}*/

					if (est_matching_len >= MIN_DISTANT_MATCH_LEN)
						matching_len = equal_len(pos, i, MIN_MATCH_LEN);
					else
						matching_len = est_matching_len;

//					if (matching_len < MIN_MATCH_LEN)
//						continue;
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
					while (!v_parsing.empty() && v_parsing.back().query_pos >= prev_region_start)
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
							prev_region_start = v_parsing[j].query_pos;
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
void CSharedWorker::prepare_ht_long()
{
	uint32_t x = (uint32_t)(v_kmers_rl.size() / htl_max_fill_factor);

	while (x & (x - 1))
		x &= x - 1;

	htl_size = 2 * x;
	htl_mask = htl_size - 1;

	if (!htl)
		htl = new vector<int>;

	htl->clear();
	htl->resize(htl_size, HT_EMPTY);

	const int pf_dist = 16;

	for (int i = 0; i + pf_dist < (int)v_kmers_rl.size(); ++i)
	{
		if(v_kmers_rl[i + pf_dist].first >= 0)
			prefetch_htl(hash_mm(v_kmers_rl[i + pf_dist].first, htl_mask));

		if (v_kmers_rl[i].first < 0)
			continue;

		auto ht_idx = hash_mm(v_kmers_rl[i].first, htl_mask);

		while ((*htl)[ht_idx] != HT_EMPTY)
			ht_idx = (ht_idx + 1) & htl_mask;

		(*htl)[ht_idx] = v_kmers_rl[i].second;
	}
	
	for (int i = max((int)v_kmers_rl.size() - pf_dist, 0); i < (int)v_kmers_rl.size(); ++i)
	{
		if (v_kmers_rl[i].first < 0)
			continue;

		auto ht_idx = hash_mm(v_kmers_rl[i].first, htl_mask);

		while ((*htl)[ht_idx] != HT_EMPTY)
			ht_idx = (ht_idx + 1) & htl_mask;

		(*htl)[ht_idx] = v_kmers_rl[i].second;
	}
}

// ****************************************************************************
void CSharedWorker::prefetch_hts1(int pos)
{
#ifdef _WIN32
//	_mm_prefetch((const char*)(hts.data() + pos), _MM_HINT_T0);
//	_mm_prefetch((const char*)(hts2.data() + pos), _MM_HINT_T0);
	_mm_prefetch((const char*)(hts3_desc->data() + pos), _MM_HINT_T0);
#else
//	__builtin_prefetch(hts->data() + pos);
//	__builtin_prefetch(hts2->data() + pos);
	__builtin_prefetch(hts3_desc->data() + pos);
#endif
}

// ****************************************************************************
void CSharedWorker::prefetch_hts2(int pos)
{
#ifdef _WIN32
	//	_mm_prefetch((const char*)(hts.data() + pos), _MM_HINT_T0);
	//	_mm_prefetch((const char*)(hts2.data() + pos), _MM_HINT_T0);
	_mm_prefetch((const char*)(hts3->data() + (*hts3_desc)[pos].first), _MM_HINT_T0);
#else
	//	__builtin_prefetch(hts->data() + pos);
//	__builtin_prefetch(hts2->data() + pos);
	__builtin_prefetch((const char*)(hts3->data() + (*hts3_desc)[pos].first), _MM_HINT_T0);
#endif
}

// ****************************************************************************
void CSharedWorker::prefetch_htl(int pos)
{
#ifdef _WIN32
	_mm_prefetch((const char*)(htl->data() + pos), _MM_HINT_T0);
#else
	__builtin_prefetch(htl->data() + pos);
#endif
}

// ****************************************************************************
void CSharedWorker::prepare_ht_short()
{
	uint32_t ht_size = 1u << (2 * MIN_MATCH_LEN);
//	uint32_t ht_mask = ht_size - 1u;

	if(!hts3_desc)
		hts3_desc = new vector<pair<int, int>>;
	if(!hts3)
		hts3 = new vector<pair<int, int>>;

//	hts.clear();
//	hts2.clear();
	hts3_desc->clear();
//	hts.resize(ht_size);
//	hts2.resize(ht_size);
	hts3_desc->resize(ht_size, make_pair(0, 0));
	hts3->resize(v_kmers_rs.size());

	const int pf_dist = 16;
	int n_kmers = (int)v_kmers_rs.size();

	for (int i = 0; i < n_kmers; ++i)
	{
		if (i + pf_dist < n_kmers && v_kmers_rs[i + pf_dist].first >= 0)
			prefetch_hts1((int)v_kmers_rs[i + pf_dist].first);

		if (v_kmers_rs[i].first >= 0)
			++(*hts3_desc)[v_kmers_rs[i].first].second;
	}

	for (int i = 1; i < ht_size; ++i)
		(*hts3_desc)[i].first = (*hts3_desc)[i - 1].first + (*hts3_desc)[i - 1].second;

	for (int i = 0; i < n_kmers; ++i)
	{
		if (i + pf_dist < n_kmers && v_kmers_rs[i + pf_dist].first >= 0)
			prefetch_hts1((int)v_kmers_rs[i + pf_dist].first);

		if (v_kmers_rs[i].first >= 0)
		{
			(*hts3)[(*hts3_desc)[v_kmers_rs[i].first].first].first = v_kmers_rs[i].second;

			if (v_kmers_rl[i].first >= 0)
				(*hts3)[(*hts3_desc)[v_kmers_rs[i].first].first].second = ((int)v_kmers_rl[i].first) & hts_mask;
			else
				(*hts3)[(*hts3_desc)[v_kmers_rs[i].first].first].second = -1;

			++(*hts3_desc)[v_kmers_rs[i].first].first;
		}
	}

	for (int i = 0; i < ht_size; ++i)
		(*hts3_desc)[i].first -= (*hts3_desc)[i].second;


/*	const int pf_dist = 16;

	int est_hts_entry_len = (int) (2 * v_kmers_rs.size() / ht_size);
	
	for (auto& x : hts2)
		x.reserve(est_hts_entry_len);

	for (int i = 0; i + pf_dist < (int)v_kmers_rs.size(); ++i)
	{
		if(v_kmers_rs[i + pf_dist].first >= 0)
			prefetch_hts((int) v_kmers_rs[i + pf_dist].first);
		if (v_kmers_rs[i].first >= 0)
		{
			//			hts[v_kmers_rs[i].first].emplace_back(v_kmers_rs[i].second);
			if (v_kmers_rl[i].first >= 0)
				hts2[v_kmers_rs[i].first].emplace_back(make_pair(v_kmers_rs[i].second, ((int)v_kmers_rl[i].first) & hts_mask));
			else
				hts2[v_kmers_rs[i].first].emplace_back(make_pair(v_kmers_rs[i].second, -1));
		}
	}

	for (int i = max((int)v_kmers_rs.size() - pf_dist, 0); i < (int)v_kmers_rs.size(); ++i)
		if (v_kmers_rs[i].first >= 0)
		{
			//			hts[v_kmers_rs[i].first].emplace_back(v_kmers_rs[i].second);
			if (v_kmers_rl[i].first >= 0)
				hts2[v_kmers_rs[i].first].emplace_back(make_pair(v_kmers_rs[i].second, ((int)v_kmers_rl[i].first) & hts_mask));
			else
				hts2[v_kmers_rs[i].first].emplace_back(make_pair(v_kmers_rs[i].second, -1));
		}*/
}

// ****************************************************************************
void CSharedWorker::prepare_kmers_ref_short()
{
	prepare_kmers(v_kmers_rs, *s_reference, MIN_MATCH_LEN, true);
}

// ****************************************************************************
void CSharedWorker::prepare_kmers_ref_long()
{
	prepare_kmers(v_kmers_rl, *s_reference, MIN_DISTANT_MATCH_LEN, true);
}

// ****************************************************************************
void CSharedWorker::prepare_kmers_data()
{
	prepare_kmers(v_kmers_ds, s_data, MIN_MATCH_LEN, true);
	prepare_kmers(v_kmers_dl, s_data, MIN_DISTANT_MATCH_LEN, true);
}

// ****************************************************************************
void CSharedWorker::clear_ref()
{
	if(s_reference)
		s_reference->clear();

	n_reference = 0;

	htl_size = 0;
	htl_mask = 0;

	if(htl)
		htl->clear();
	if(hts)
		hts->clear();
	if(hts2)
		hts2->clear();

	if(hts3)
		hts3->clear();
	if(hts3_desc)	
		hts3_desc->clear();
}

// ****************************************************************************
void CSharedWorker::clear_data()
{
	s_data.clear();

	n_data = 0;

	v_parsing.clear();
}

// EOF
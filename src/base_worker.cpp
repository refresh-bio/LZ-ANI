#include "base_worker.h"
#include "distributions.h"

#include <algorithm>
#include <iostream>
#include <functional>
#include <fstream>
#include <string>
#include <iterator>
#include <cstring>

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
BaseWorker::BaseWorker() : bufferSize(INITIAL_BUFFER_SIZE) {
	fill_n(codes, 256, 4);
	codes['A'] = 0;
	codes['C'] = 1;
	codes['G'] = 2;
	codes['T'] = 3;

	hts_mask = (int)(1u << (2 * (MIN_DISTANT_MATCH_LEN - MIN_MATCH_LEN))) - 1;

	loadBuffer = new char[bufferSize];
}

// ****************************************************************************
BaseWorker::~BaseWorker() {
	delete[] loadBuffer;
}


// ****************************************************************************
void BaseWorker::prefetch(int pos)
{
#ifdef _WIN32
	_mm_prefetch((const char*)(s_reference->data() + pos), _MM_HINT_T0);
#else
	__builtin_prefetch(s_reference->data() + pos);
#endif
}

// ****************************************************************************
void BaseWorker::compare_ranges(int data_start_pos, int ref_start_pos, int len, bool backward)
{
	int r_len = 0;
	bool is_matching = false;
	flag_t flag = backward ? flag_t::match_distant : flag_t::match_close;

	for (int j = 0; j < len; ++j)
	{
		if ((*s_reference)[ref_start_pos + j] == s_data[data_start_pos + j])
		{
			if (is_matching)
				r_len++;
			else
			{
				if (r_len)
					v_parsing.emplace_back(CFactor(data_start_pos + j - r_len, flag_t::run_literals, 0, r_len, 0));
				r_len = 1;
				is_matching = true;
			}
		}
		else
		{
			if (is_matching)
			{
				v_parsing.emplace_back(CFactor(data_start_pos + j - r_len, flag, ref_start_pos + j - r_len, r_len, 0));
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
	else if (r_len)
		v_parsing.emplace_back(CFactor(data_start_pos + len - r_len, flag_t::run_literals, 0, r_len, 0));
}

// ****************************************************************************
int BaseWorker::try_extend_forward(int data_start_pos, int ref_start_pos)
{
	int data_size = (int)s_data.size();
	int ref_size = (int)s_reference->size();

	int approx_ext;
	int no_missmatches = 0;
	int last_match = 0;
	vector<int> window(APPROX_WINDOW, 0);

	for (approx_ext = 0; data_start_pos + approx_ext < data_size && ref_start_pos + approx_ext < ref_size; ++approx_ext)
	{
		bool is_missmatch = s_data[data_start_pos + approx_ext] != (*s_reference)[ref_start_pos + approx_ext];
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
int BaseWorker::try_extend_forward2(int data_start_pos, int ref_start_pos)
{
	int data_size = (int)s_data.size();
	int ref_size = (int)s_reference->size();

	int approx_ext;
	int no_missmatches = 0;
	int last_run_match = 0;
	vector<int> window(APPROX_WINDOW, 0);
	int match_run_len = APPROX_RUNLEN;

	for (approx_ext = 0; data_start_pos + approx_ext < data_size && ref_start_pos + approx_ext < ref_size; ++approx_ext)
	{
		bool is_missmatch = s_data[data_start_pos + approx_ext] != (*s_reference)[ref_start_pos + approx_ext];
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
int BaseWorker::try_extend_backward(int data_start_pos, int ref_start_pos, int max_len)
{
	int approx_ext;
	int no_missmatches = 0;
	int last_match = 0;
	vector<int> window(APPROX_WINDOW, 0);

	for (approx_ext = 0; data_start_pos - approx_ext > 0 && ref_start_pos - approx_ext > 0 && approx_ext < max_len; ++approx_ext)
	{
		bool is_missmatch = s_data[data_start_pos - approx_ext - 1] != (*s_reference)[ref_start_pos - approx_ext - 1];
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
int BaseWorker::try_extend_backward2(int data_start_pos, int ref_start_pos, int max_len)
{
	int approx_ext;
	int no_missmatches = 0;
	int last_run_match = 0;
	vector<int> window(APPROX_WINDOW, 0);
	int match_run_len = APPROX_RUNLEN;

	for (approx_ext = 0; data_start_pos - approx_ext > 0 && ref_start_pos - approx_ext > 0 && approx_ext < max_len; ++approx_ext)
	{
		bool is_missmatch = s_data[data_start_pos - approx_ext - 1] != (*s_reference)[ref_start_pos - approx_ext - 1];
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
bool BaseWorker::load_file(const string &file_name, Genome& genome, int separator)
{
//	cout << "LOAD: " << file_name << endl;
	
	genome.seq.clear();

	std::ifstream f(file_name.c_str());
	if (!f)
		return false;

	// adjust size of buffer
	f.seekg(0, std::ios::end);
	size_t size = f.tellg();
	f.seekg(0, std::ios::beg);

	if (size + 1 > bufferSize) {
		bufferSize = std::max(bufferSize * 2, size + 1);
		delete[] loadBuffer;
		loadBuffer = new char[bufferSize];
	}
	
	f.read(loadBuffer, size);

	loadBuffer[size] = 0; // add null termination 
	
	f.close();
	
	std::vector<char*> subsequences;
	std::vector<char*> headers;
	std::vector<size_t> lengths;
	
	extractSubsequences(loadBuffer, size, subsequences, lengths, headers);

	genome.totalLen = size;
	genome.seq.resize(size + subsequences.size() * Genome::SEPARATOR_LENGTH);

	auto out = genome.seq.begin();
	for (int i = 0; i < subsequences.size(); ++i) {
		std::copy(subsequences[i], subsequences[i] + lengths[i], out);
		out += lengths[i];

		std::fill(out, out + Genome::SEPARATOR_LENGTH, separator);
		out += Genome::SEPARATOR_LENGTH;
	}

	genome.seq.erase(out, genome.seq.end());
	
	// fill in genome structure 
	genome.lengths = lengths;
	for (auto h : headers) {
		genome.headers.emplace_back(h);
	}
	
	return true;
}


// ****************************************************************************
bool BaseWorker::extractSubsequences(
	char* data,
	size_t& totalLen,
	std::vector<char*>& subsequences,
	std::vector<size_t>& lengths,
	std::vector<char*>& headers) {

	// extract contigs
	char * header = nullptr;
	char * ptr = data;

	while (header = strchr(ptr, '>')) { // find begining of header
		*header = 0; // put 0 as separator (end of previous chromosome)
		if (subsequences.size()) {
			lengths.push_back(header - subsequences.back());
		}

		++header; // omit '<'
		headers.push_back(header);

		ptr = strchr(header, '\n'); // find end of header
		if (*(ptr - 1) == '\r') { // on Windows
			*(ptr - 1) = 0;
		}
		*ptr = 0; // put 0 as separator
		++ptr; // move to next character (begin of chromosome)
		subsequences.push_back(ptr); // store chromosome
	}

	lengths.push_back(data + totalLen - subsequences.back());

	// remove newline characters from chromosome
	totalLen = 0;
	for (int i = 0; i < subsequences.size(); ++i) {
		// determine chromosome end
		char* newend = std::remove_if(subsequences[i], subsequences[i] + lengths[i], [](char c) -> bool { return c == '\n' || c == '\r';  });
		*newend = 0;
		lengths[i] = newend - subsequences[i];
		totalLen += lengths[i];
		//	assert(lengths[i] == strlen(subsequences[i]));
	}

	return true;
}

// ****************************************************************************
void BaseWorker::duplicate_rev_comp(seq_t &seq)
{
	int size = (int)seq.size();
	int separator = seq.back();

	seq.resize(2 * size);

	auto complement = [](int8_t in) -> int8_t {
		if (in == sym_A)		return sym_T;
		else if (in == sym_C)	return sym_G;
		else if (in == sym_G)	return sym_C;
		else if (in == sym_T)	return sym_A;
		else return in;
	};

	std::transform(seq.begin(), seq.begin() + size, seq.rbegin(), complement);

	// add separators at the end of the duplicated sequence
	seq.resize(seq.size() + Genome::SEPARATOR_LENGTH);
	std::fill(seq.rbegin(), seq.rbegin() + Genome::SEPARATOR_LENGTH, separator);
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


// ****************************************************************************
void BaseWorker::calc_ani(CResults &res, int mode, std::vector<Region>& v_matches)
{

	int cur_match_len = 0;
	int cur_match_lit = 0;
	int n_lit = 0;

	CFactor *firstInRegion = nullptr;
	CFactor *lastInRegion = nullptr;

	for (auto& x : v_parsing)
	{
		if (x.flag == flag_t::match_distant)
		{
			// store a
			if (cur_match_len) {
				v_matches.emplace_back(Region(
					cur_match_len, cur_match_lit,
					firstInRegion->query_pos,
					firstInRegion->ref_pos, lastInRegion->ref_pos + lastInRegion->len - firstInRegion->ref_pos));
			}

			cur_match_len = x.len;
			cur_match_lit = 0;
			n_lit = 0;

			firstInRegion = lastInRegion = &x; // update last as well in the case there is only one match
		}
		else if (x.flag == flag_t::match_close)
		{
			cur_match_len += x.len;
			cur_match_lit += n_lit;
			n_lit = 0;
			lastInRegion = &x;
		}
		else if (x.flag == flag_t::run_literals)
		{
			n_lit += x.len;
		}
	}

	if (cur_match_len) 
		v_matches.emplace_back(Region(
			cur_match_len, cur_match_lit,
			firstInRegion->query_pos,
			firstInRegion->ref_pos, lastInRegion->ref_pos + lastInRegion->len - firstInRegion->ref_pos));


	sort(v_matches.begin(), v_matches.end());
	
	int ref_len = (int)s_reference->size() - n_reference * Genome::SEPARATOR_LENGTH;
	int data_len = (int)s_data.size() - n_data * Genome::SEPARATOR_LENGTH;
	int n_sym_in_matches = 0;
	int n_sym_in_literals = 0;

	for (auto& x : v_matches) {
		if (x.num_matches + x.num_literals >= MIN_REGION_LEN) // && (double) x.first / (x.first + x.second) > 0.5
		{
			n_sym_in_matches += x.num_matches;
			n_sym_in_literals += x.num_literals;
		}
	}

	if (mode == 1)
	{
		res.ref_size = ref_len / 2;
		res.query_size = data_len;
	}
	res.sym_in_literals[mode] = n_sym_in_literals;
	res.sym_in_matches[mode] = n_sym_in_matches;
	res.coverage[mode] = (double)(n_sym_in_literals + n_sym_in_matches) / data_len;
	if (n_sym_in_matches + n_sym_in_literals)
		res.ani[mode] = (double)n_sym_in_matches / (n_sym_in_matches + n_sym_in_literals);
	else
		res.ani[mode] = 0.0;

	/*	if (res.coverage[mode] < MIN_COVERAGE)
	res.ani[mode] -= 0.4;*/

	if (res.ani[mode] > 1.0)
		res.ani[mode] = 1.0;
	if (res.coverage[mode] > 1.0)
		res.coverage[mode] = 1.0;

	// calculate p-values for all matches
	double p_success = res.ani[mode];

	for (auto& match : v_matches) {
		// take into account only high-conserved regions
		int span = match.num_matches + match.num_literals;

		if ((double)match.num_matches / span > p_success) {
			// use binomial distribution with global ANI as the probability of success (symbol conservation) 
			BinomialDistributionApproximation dist(span, p_success);

			// calculate p-value as the probability of obtaining same and more matching symbols then observed for particular region
			match.p_value = 1 - dist.cdf((double)match.num_matches - 0.5);
		}
	}

	std::sort(v_matches.begin(), v_matches.end(), [](const Region& lhs, const Region& rhs)->bool {
		return lhs.p_value < rhs.p_value;
	});

}

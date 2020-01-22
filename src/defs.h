#pragma once

#include <cstdint>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;

const uint8_t sym_A = 'A';
const uint8_t sym_C = 'C';
const uint8_t sym_G = 'G';
const uint8_t sym_T = 'T';
const uint8_t sym_N1 = 'N';
const uint8_t sym_N2 = 'n';

const int DEF_MIN_MATCH_LEN = 8;
const int DEF_MIN_CLOSE_MATCH_LEN = 8;
const int DEF_MIN_DISTANT_MATCH_LEN = 20;
const int DEF_CLOSE_DIST = 256;
//const int DEF_LONG_LITERAL_RUN_LEN = 16;
const int DEF_MAX_LIT_RUN_IN_MATCH = 32;
const double DEF_MIN_COVERAGE = 0.10;
const int DEF_MIN_REGION_LEN = 128;
const int DEF_APPROX_WINDOW = 16;
const int DEF_APPROX_MISMATCHES = 5;
const int DEF_APPROX_RUNLEN = 3;
const int DEF_RANGE_FROM = 0;
const int DEF_RANGE_TO = 1 << 30;

enum class flag_t {match, match_close, match_distant, literal, run_literals, match_literal};

const int HT_EMPTY = -1;
const int HT_FAIL = -1;

struct CFactor {
	int query_pos;
	flag_t flag;
	int ref_pos;
	int len;
	uint8_t symbol;

	CFactor(int _data_pos, flag_t _flag, uint32_t _offset, uint32_t _len, uint8_t _symbol) :
		query_pos(_data_pos), flag(_flag), ref_pos(_offset), len(_len), symbol(_symbol)
	{}
};

typedef vector<uint8_t> seq_t;

struct CResults {
	int ref_size;
	int query_size;
	int sym_in_matches[3];
	int sym_in_literals[3];
	double coverage[3];
	double ani[3];
	double time;

	CResults() = default;
};

struct Region {
	int num_matches;
	int num_literals;

	int query_pos;

	int ref_raw_pos;
	int ref_len;
	
	double p_value;
	
	Region(int num_matches, int num_literals, 
		int query_pos,
		int ref_raw_pos, int ref_len) : 
			num_matches(num_matches), num_literals(num_literals), 
			query_pos(query_pos),
			ref_raw_pos(ref_raw_pos), ref_len(ref_len), p_value(1)
	{}

	bool operator<(const Region& rhs) {
		return this->num_matches < rhs.num_matches || ((this->num_matches == rhs.num_matches) && (this->num_literals < rhs.num_literals));
	}
};


struct Genome {
	static const size_t SEPARATOR_LENGTH = 1000; 

	seq_t seq;
	std::vector<size_t> lengths;
	std::vector<string> headers;
	size_t totalLen;

	const size_t n_seqs() const { return lengths.size();  }

	void recodeNs(uint8_t from_sym, uint8_t to_sym)
	{
		std::replace(seq.begin(), seq.end(), from_sym, to_sym);
	}

	bool translateRawPosition(size_t rawPos, size_t span, size_t &subsequenceId, size_t& startPos, size_t& endPos) {
		
		size_t half = totalLen + (lengths.size() ) * SEPARATOR_LENGTH;
	
		if (rawPos < half) {
		
			startPos = rawPos;
			// start from the first sequence
			for (subsequenceId = 0; startPos > lengths[subsequenceId]; ++subsequenceId) {
				startPos -= lengths[subsequenceId] + SEPARATOR_LENGTH;
			}

			endPos = startPos + span;
		}
		else if (rawPos < 2 * half ) {
			
			startPos = rawPos - (half + SEPARATOR_LENGTH); // seek after forward-revcompl separator
			for (subsequenceId = n_seqs() - 1; startPos > lengths[subsequenceId]; --subsequenceId) {
				startPos -= lengths[subsequenceId] + SEPARATOR_LENGTH;
			}
			startPos = lengths[subsequenceId] - startPos;
			endPos = startPos - span;		
		}
		else {
			return false;
		}

		// position inside separator area
		if (startPos > lengths[subsequenceId] || endPos > lengths[subsequenceId]) {
			return false;
		}

		return true;
	}

	std::string printDebugInfo() {
		std::ostringstream oss;
		oss << "Chomosomes length: " << totalLen << endl
			<< "Raw length: " << seq.size() << endl
			<< "#A: " << std::count(seq.begin(), seq.end(), 'A') << endl
			<< "#C: " << std::count(seq.begin(), seq.end(), 'C') << endl
			<< "#G: " << std::count(seq.begin(), seq.end(), 'G') << endl
			<< "#T: " << std::count(seq.begin(), seq.end(), 'T') << endl
			<< "#N: " << std::count(seq.begin(), seq.end(), 'N') << endl
			<< "#n: " << std::count(seq.begin(), seq.end(), 'n') << endl;
		size_t start = 0;
		for (int i = 0; i < n_seqs(); ++i) {
			size_t end = start + lengths[i];
			oss << "Chr " << i << ": " << start << "-" << end << " (" << lengths[i] << ")" << endl;
			start = end + SEPARATOR_LENGTH;
		}

		return oss.str();
		
	}

};

// EOF

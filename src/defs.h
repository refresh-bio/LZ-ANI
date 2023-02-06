#pragma once

#include <cstdint>
#include <vector>

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
	int data_pos;
	flag_t flag;
	int offset;
	int len;
	uint8_t symbol;

	CFactor(int _data_pos, flag_t _flag, uint32_t _offset, uint32_t _len, uint8_t _symbol) :
		data_pos(_data_pos), flag(_flag), offset(_offset), len(_len), symbol(_symbol)
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
	double total_ani;

	CResults() :
		ref_size(0),
		query_size(0),
		sym_in_matches{ 0,0,0 },
		sym_in_literals{ 0,0,0 },
		coverage{ 0, 0, 0 },
		ani{ 0,0,0 },
		time{ 0 },
		total_ani{ 0 }
	{};
};

// EOF

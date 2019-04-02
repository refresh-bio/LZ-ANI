#pragma once

#include <cstdint>
#include <vector>

using namespace std;

const uint8_t sym_A = 'A';
const uint8_t sym_C = 'C';
const uint8_t sym_G = 'G';
const uint8_t sym_T = 'T';
const uint8_t sym_N = 'N';

const int MIN_MATCH_LEN = 9;
const int MIN_CLOSE_MATCH_LEN = 9;
const int MIN_DISTANT_MATCH_LEN = 22;
const int CLOSE_DIST = 128;
const int LONG_LITERAL_RUN_LEN = 16;
const int MAX_LIT_RUN_IN_MATCH = 128;

enum class flag_t {match, match_close, match_distant, literal, run_literals};

const int HT_EMPTY = -1;
const int HT_FAIL = -1;

struct CFactor {
	flag_t flag;
	int offset;
	int len;
	uint8_t symbol;

	CFactor(flag_t _flag, uint32_t _offset, uint32_t _len, uint8_t _symbol) :
		flag(_flag), offset(_offset), len(_len), symbol(_symbol)
	{}
};

typedef vector<uint8_t> seq_t;


// EOF

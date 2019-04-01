#pragma once

#include <cstdint>
#include <vector>

using namespace std;


const uint32_t close_dist = 64;

const uint8_t sym_A = 'A';
const uint8_t sym_C = 'C';
const uint8_t sym_G = 'G';
const uint8_t sym_T = 'T';
const uint8_t sym_N = 'N';

const uint32_t MIN_MATCH_LEN = 8;
const uint32_t MIN_CLOSE_MATCH_LEN = 3;
const uint32_t MIN_DISTANT_MATCH_LEN = 12;

enum class flag_t {match, match_close, match_distant, literal, run_literals};

const uint32_t HT_EMPTY = ~0u;
const uint32_t HT_FAIL = ~0u;

struct CFactor {
	flag_t flag;
	uint32_t offset;
	uint32_t len;
	uint8_t symbol;

	CFactor(flag_t _flag, uint32_t _offset, uint32_t _len, uint8_t _symbol) :
		flag(_flag), offset(_offset), len(_len), symbol(_symbol)
	{}
};

typedef vector<uint8_t> seq_t;


// EOF

#pragma once

#include "params.h"
#include "defs.h"
#include "seq_reservoir.h"

class CLZMatcher2
{
	CParams2 params;

	CSeqReservoir seq_reservoir;



public:
	CLZMatcher2(CParams2& params) :
		params(params)
	{}

	bool load_sequences();
};
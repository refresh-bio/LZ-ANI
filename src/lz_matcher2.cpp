#include "lz_matcher2.h"

// ****************************************************************************
bool CLZMatcher2::load_sequences()
{
	if(params.multisample_fasta)
		return seq_reservoir.load_multifasta(params.input_file_names);
	else
		return seq_reservoir.load_fasta(params.input_file_names);
}

// EOF
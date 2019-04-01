#pragma once

#include "defs.h"
#include <vector>
#include <string>

bool load_file(const string &file_name, seq_t &seq, uint32_t &n_parts);
void duplicate_rev_comp(seq_t &seq);

// EOF

#pragma once

#include <vector>
#include <string>

#include "params.h"

using namespace std;



class CLZMatcher
{
	CParams& params;

	vector<string> input_file_names;

public:
	CLZMatcher(CParams& params) : params(params)
	{};

	bool run_all2all(vector<string>& _input_file_names);
};
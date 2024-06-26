// *******************************************************************************************
// This file is a part of LZ-ANI software distributed under GNU GPL 3 license.
// The homepage of the LZ-ANI project is https://github.com/refresh-bio/LZ-ANI
//
// Copyright(C) 2024-2024, S.Deorowicz, A.Gudys
//
// Version: 1.0.0
// Date   : 2024-06-26
// *******************************************************************************************

#include "utils.h"
#include <algorithm>

// ****************************************************************************
vector<string> split(const string& str, char sep)
{
	vector<string> parts;

	string s;

	for (auto c : str)
	{
		if (c == sep)
		{
			parts.emplace_back(s);
			s.clear();
		}
		else
			s.push_back(c);
	}

	if (!s.empty())
		parts.emplace_back(s);

	return parts;
}

// ****************************************************************************
string strip_at_space(const string& str)
{
	auto p = find(str.begin(), str.end(), ' ');

	if (p == str.end())
		return str;
	else
		return string(str.begin(), p);
}

// EOF

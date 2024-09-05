// *******************************************************************************************
// This file is a part of LZ-ANI software distributed under GNU GPL 3 license.
// The homepage of the LZ-ANI project is https://github.com/refresh-bio/LZ-ANI
//
// Copyright(C) 2024-2024, S.Deorowicz, A.Gudys
//
// Version: 1.1.0
// Date   : 2024-09-05
// *******************************************************************************************

#pragma once

#include <vector>
#include <string>
#include <cassert>

using namespace std;

vector<string> split(const string& str, char sep);
string strip_at_space(const string& str);

// EOF
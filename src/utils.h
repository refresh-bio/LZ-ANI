#pragma once

#include <vector>
#include <string>

using namespace std;

vector<string> split(const string& str, char sep);

string remove_path_from_file(const string& file_path);

string strip_at_space(const string& str);

// EOF
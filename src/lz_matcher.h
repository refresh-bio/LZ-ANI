#pragma once

#include <cinttypes>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <thread>
#include <future>
#include <chrono>

#include "params.h"
#include "defs.h"
#include "s_worker.h"
#include "data_storage.h"
#include "conversion.h"

using namespace std;
using namespace std::chrono;

using filter_dict_t = unordered_set<pair_id_t>;
//using filter_dict_t = set<pair_id_t>;

//using results_dict_t = unordered_map<pair_id_t, CResults>;
using results_dict_t = map<pair_id_t, CResults>;


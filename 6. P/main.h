#ifndef MAIN_H
#define MAIN_H

#include "helper_functions.h"
#include "model.h"
#include <functional>
#include <unordered_map>

using ModelFunction = std::function<void(Instance&, Solution&)>;

#endif
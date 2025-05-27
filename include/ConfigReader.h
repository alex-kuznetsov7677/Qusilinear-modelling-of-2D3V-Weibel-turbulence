#pragma once

#include <map>
#include <string>
#include <fstream>
#include <sstream>

std::map<std::string, double> ReadConfig(const std::string& filename);

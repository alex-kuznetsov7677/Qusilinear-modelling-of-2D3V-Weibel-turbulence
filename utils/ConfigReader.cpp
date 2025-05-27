#include <map>
#include <string>
#include <fstream>
#include <sstream>

#include "ConfigReader.h"

std::map<std::string, double> ReadConfig(const std::string& filename) {
    std::ifstream file(filename);
    std::map<std::string, double> config;
    std::string line;

    while (std::getline(file, line)) {

        if (line.empty()) continue;

        std::istringstream is_line(line);
        std::string key;
        if (std::getline(is_line, key, '=')) {
            std::string value;
            if (std::getline(is_line, value)) {
                key.erase(0, key.find_first_not_of(" \t"));
                key.erase(key.find_last_not_of(" \t") + 1);
                value.erase(0, value.find_first_not_of(" \t"));
                value.erase(value.find_last_not_of(" \t") + 1);

                config[key] = std::stod(value);
            }
        }
    }

    return config;
}
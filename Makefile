CXX = mpicxx
CXXFLAGS = -std=c++11 -O3 -march=native -ffast-math -Iinclude
LDFLAGS = -lm -lmpi

SRC_DIR = src
UTILS_DIR = utils
BUILD_DIR = build

SOURCES = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(UTILS_DIR)/*.cpp)
OBJECTS = $(addprefix $(BUILD_DIR)/, $(notdir $(SOURCES:.cpp=.o)))
TARGET = $(BUILD_DIR)/exe9_c

all: $(BUILD_DIR) output $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(UTILS_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

output:
	mkdir -p output

clean:
	rm -rf $(BUILD_DIR) output

rebuild: clean all

.PHONY: all clean rebuild output

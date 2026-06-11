# Makefile

CXX := nvcc
CXXFLAGS :=

# Directories
SRC_DIR := src
GRID_DIR := grid
BOUNDARY_DIR := boundary
CUDA_DIR := CUDA_files
TESTCASE_DIR := testcase
INCLUDE_DIR := include
BUILD_DIR := build
BIN_DIR := bin
SCHEME_DIR := schemes
TEMPORAL_DIR := temporal
SOLUTION_DIR := solution

# Source files
SRC_FILES := $(wildcard $(SRC_DIR)/*.cu) $(wildcard $(GRID_DIR)/*.cpp) $(wildcard $(BOUNDARY_DIR)/*.cpp) $(wildcard $(CUDA_DIR)/*.cu) $(wildcard $(TESTCASE_DIR)/testcase1.cpp) $(wildcard $(SCHEME_DIR)/llf.cu) $(wildcard $(TEMPORAL_DIR)/*.cpp) $(wildcard $(SOLUTION_DIR)/*.cpp)
OBJ_FILES := $(addprefix $(BUILD_DIR)/,$(notdir $(SRC_FILES:.cu=.o)))
OBJ_FILES := $(OBJ_FILES:.cpp=.o)

# Executable name
EXECUTABLE := $(BIN_DIR)/main

# Include directory
INC_DIR := $(INCLUDE_DIR)

# Targets
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) $^ -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

$(BUILD_DIR)/%.o: $(GRID_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

$(BUILD_DIR)/%.o: $(BOUNDARY_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

$(BUILD_DIR)/%.o: $(TESTCASE_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

$(BUILD_DIR)/%.o: $(SCHEME_DIR)/%.cu
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

$(BUILD_DIR)/%.o: $(TEMPORAL_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

$(BUILD_DIR)/%.o: $(SOLUTION_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

$(BUILD_DIR)/%.o: $(CUDA_DIR)/%.cu
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

clean:
	rm -f $(BUILD_DIR)/*.o $(EXECUTABLE) $(BIN_DIR)/finalSolution.vtk

.PHONY: all clean
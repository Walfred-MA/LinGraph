CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic

TARGET := build_local_graphs_parallel_hotspots
SOURCE := build_local_graphs_parallel_hotspots.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -pthread -o $@ $<

clean:
	rm -f $(TARGET)

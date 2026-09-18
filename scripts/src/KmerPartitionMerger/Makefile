CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
LDFLAGS ?=

# Conda and cluster environments often predefine CXXFLAGS/LDFLAGS.  Append
# pthread unconditionally so those environment values cannot discard the
# compile and link options required by std::thread.
override CXXFLAGS += -pthread
override LDLIBS += -pthread

TARGET := kmerpartitionmerger
SOURCES := $(wildcard src/*.cpp)
OBJECTS := $(SOURCES:.cpp=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

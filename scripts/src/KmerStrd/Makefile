CXX := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -pthread

TARGET := kmerstrd
SRCS := src/main.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -pthread

clean:
	rm -f $(TARGET)

CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS ?=
LDFLAGS ?=
LDLIBS ?= -lz
# Baseline x86-64: a host-tuned default emits AVX that older nodes reject (SIGILL).
ARCHFLAGS ?= $(if $(filter x86_64 amd64,$(shell uname -m)),-march=x86-64 -mtune=generic)

SmallVariantMerge: SmallVariantMerge.cpp support.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ARCHFLAGS) -pthread SmallVariantMerge.cpp -o $@ $(LDFLAGS) $(LDLIBS)

.PHONY: clean
clean:
	rm -f SmallVariantMerge

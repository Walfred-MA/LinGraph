# === Compiler & Flags ===
CXX      := g++
CXXFLAGS := -Wall -std=c++17 -O3 -flto -DCOUNTINT=uint8_t

# === Libraries & Linker Options ===
# Note: Keep -Wl,-rpath flags so the executable knows where Spack libraries live at runtime
LDFLAGS  := -Wl,-rpath,/apps/spack/2406/apps/linux-rocky8-x86_64_v3/gcc-13.3.0/curl-8.8.0-sqn4mlx/lib
LDLIBS   := -lhts -lcurl -lz -lpthread

# === Source Files ===
SRC_DIR  := src
SRCS     := $(wildcard $(SRC_DIR)/*.cpp)
TARGET   := kmer_searcher

# === Rules ===
all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(TARGET) *.o

.PHONY: all clean

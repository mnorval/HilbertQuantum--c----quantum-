CXX?=g++
CXXFLAGS=-std=c++17 -O2 -Wall -Wextra -Iinclude
all:
	mkdir -p build
	$(CXX) $(CXXFLAGS) src/main.cpp -o build/hilbert
clean:
	rm -rf build

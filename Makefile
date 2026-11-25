# Build the two experiment drivers.
#
# The bin/ directory is part of the repository, so no directory creation is required and the
# recipes work under both POSIX make and mingw32-make on Windows.

CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17
INCLUDES  = -Iinclude
BIN       = bin

COMMON_SRC = src/DotsException.cpp
HEADERS    = $(wildcard include/*.h) $(wildcard include/dim/*.h)

ifeq ($(OS),Windows_NT)
EXE = .exe
else
EXE =
endif

.PHONY: all clean test

all: $(BIN)/run_unified$(EXE) $(BIN)/run_bisect$(EXE)

$(BIN)/run_unified$(EXE): src/run_unified.cpp $(COMMON_SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ src/run_unified.cpp $(COMMON_SRC)

$(BIN)/run_bisect$(EXE): src/run_bisect.cpp $(COMMON_SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ src/run_bisect.cpp $(COMMON_SRC)

test: all
	python3 scripts/smoke_test.py

clean:
	-$(RM) $(BIN)/run_unified$(EXE) $(BIN)/run_bisect$(EXE)
	-$(RM) results/*.txt results/*.csv

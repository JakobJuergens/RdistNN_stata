# Makefile — RdistNN Stata plugin (Linux / macOS)
#
# Usage:
#   make           — build rdistnn.plugin
#   make clean     — remove build artifacts
#
# The compiled rdistnn.plugin must be placed on Stata's adopath so that
# `program define rdistnn_plugin, plugin` can find it.

CXX      := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -fPIC

UNAME := $(shell uname -s)
ifeq ($(UNAME), Darwin)
  # macOS: bundle (no -lc, Stata provides its own symbols)
  LDFLAGS := -bundle -undefined dynamic_lookup
else
  # Linux
  LDFLAGS := -shared
endif

TARGET := rdistnn.plugin
CORE_TEST := rdistnn_core_test
SRCS   := stplugin.c rdistnn_plugin.cpp

$(TARGET): $(SRCS) rdistnn_core.hpp stplugin.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(SRCS)

$(CORE_TEST): test_rdistnn_core.cpp rdistnn_core.hpp
	$(CXX) $(CXXFLAGS) -o $@ test_rdistnn_core.cpp

check: $(CORE_TEST)
	./$(CORE_TEST)

clean:
	rm -f $(TARGET) $(CORE_TEST)

.PHONY: check clean

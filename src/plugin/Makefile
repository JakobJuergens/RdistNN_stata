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
  CPPFLAGS := -DSYSTEM=APPLEMAC
  # macOS: bundle (no -lc, Stata provides its own symbols)
  LDFLAGS := -bundle -undefined dynamic_lookup
else
  CPPFLAGS := -DSYSTEM=OPUNIX
  # Linux
  LDFLAGS := -shared
endif

TARGET := rdistnn.plugin
CORE_TEST := rdistnn_core_test
DISTANCE_TEST := rdistnn_distance_test
PROTOCOL_TEST := rdistnn_distance_plugin_test
ROUNDING_TEST := rdistnn_distance_rounding_test
BLOCK_TEST := rdistnn_blocks_test
LOO_PROTOCOL_TEST := rdistnn_loo_plugin_test
HEADERS := rdistnn_core.hpp rdistnn_distance.hpp compensated_arithmetic.hpp
SRCS   := stplugin.c rdistnn_plugin.cpp

$(TARGET): $(SRCS) $(HEADERS) stplugin.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $(SRCS)

$(CORE_TEST): test_rdistnn_core.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ test_rdistnn_core.cpp

$(DISTANCE_TEST): test_rdistnn_distance.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ test_rdistnn_distance.cpp

$(PROTOCOL_TEST): test_rdistnn_distance_plugin.cpp rdistnn_plugin.cpp $(HEADERS) stplugin.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ test_rdistnn_distance_plugin.cpp

$(ROUNDING_TEST): test_rdistnn_distance_rounding.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ test_rdistnn_distance_rounding.cpp

$(BLOCK_TEST): test_rdistnn_blocks.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ test_rdistnn_blocks.cpp


$(LOO_PROTOCOL_TEST): test_rdistnn_loo_plugin.cpp rdistnn_plugin.cpp $(HEADERS) stplugin.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ test_rdistnn_loo_plugin.cpp

check: $(CORE_TEST) $(DISTANCE_TEST) $(PROTOCOL_TEST) $(ROUNDING_TEST) $(BLOCK_TEST) $(LOO_PROTOCOL_TEST)
	./$(CORE_TEST)
	./$(DISTANCE_TEST) 2 2 0.5790338377234229 1.4758492921057114
	./$(PROTOCOL_TEST)
	./$(ROUNDING_TEST) 2 2 0.5790338377234229 1.4758492921057114
	./$(BLOCK_TEST)
	./$(LOO_PROTOCOL_TEST)

clean:
	rm -f $(TARGET) $(CORE_TEST) $(DISTANCE_TEST) $(PROTOCOL_TEST) $(ROUNDING_TEST) $(BLOCK_TEST) $(LOO_PROTOCOL_TEST)

.PHONY: check clean

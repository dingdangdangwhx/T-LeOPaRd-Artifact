CC ?= gcc
CXX ?= g++
CFLAGS ?= -O3 -march=native -mtune=native -Wall -Wextra
PROFILE ?= partial72
HEXL_DIR := third_party/hexl
HEXL_BUILD := build/hexl
HEXL_LIB := $(HEXL_BUILD)/hexl/lib/libhexl.a
LAZER_SOURCES := $(wildcard src/lazer/*.c src/lazer/*.h)
HEXL_SOURCES := $(shell find $(HEXL_DIR) third_party/cpu_features -type f)

.PHONY: all lib-static benchmark check clean parameters proof-parameters client-parameters
all: lib-static
	$(MAKE) -C src/protocol PROFILE=$(PROFILE) all

lib-static: liblazer.a

$(HEXL_LIB): $(HEXL_SOURCES)
	cmake -S $(HEXL_DIR) -B $(HEXL_BUILD) -DCMAKE_BUILD_TYPE=Release \
	  -DHEXL_TESTING=OFF -DHEXL_BENCHMARK=OFF -DHEXL_DOCS=OFF \
	  -DCMAKE_DISABLE_FIND_PACKAGE_CpuFeatures=ON -DCMAKE_INSTALL_LIBDIR=lib
	cmake --build $(HEXL_BUILD) --target hexl --parallel $(or $(JOBS),4)

build/lazer.o: $(LAZER_SOURCES) lazer.h config.h
	@mkdir -p build
	$(CC) $(CFLAGS) -DNDEBUG -I. -Isrc/lazer -c src/lazer/lazer.c -o $@

build/hexl.o: src/lazer/hexl.cpp src/lazer/hexl.h $(HEXL_LIB)
	$(CXX) $(CFLAGS) -DNDEBUG -Isrc/lazer -I$(HEXL_DIR)/hexl/include -c $< -o $@

liblazer.a: build/lazer.o build/hexl.o
	$(AR) rcs $@ $^

benchmark:
	bash run.sh

check: lib-static
	$(MAKE) -C src/protocol PROFILE=partial72 check
	$(MAKE) -C src/protocol PROFILE=plain72 check
	$(MAKE) -C src/protocol PROFILE=partial72 check-client
	$(MAKE) -C src/protocol PROFILE=plain72 check-client

client-parameters:
	@mkdir -p results
	python3 parameters/client_screen.py > results/client-screen.json
	sage parameters/client_commitment.sage > results/client-commitment.json
	sage parameters/client_challenge.sage > results/client-challenge.json

parameters:
	sage -python parameters/compute.py --kappa 32 --output results/parameters.json

proof-parameters:
	bash parameters/regenerate_nizk2.sh

clean:
	$(MAKE) -C src/protocol clean
	rm -rf build liblazer.a
	rm -f third_party/hexl/pkgconfig/hexl.pc third_party/hexl/cmake/hexl-*/HEXLConfig.cmake
	rm -f parameters/*.sage.py
	find parameters third_party -type d -name __pycache__ -prune -exec rm -rf {} +

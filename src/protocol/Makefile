PROFILE ?= partial72
ROOT := ../..
NIZK1 := $(ROOT)/third_party/nizk1
NIZK2 := $(ROOT)/third_party/nizk2
CC ?= gcc

ifeq ($(PROFILE),semi_partial56)
LOGQ = 56
M = 32
ELL = 35
BASE = 255
KMSIS = 25
M2 = 62
SIGMA1 = 26
SEED = 2
OUTER_BASE = 32768
FG_EXTRA = 2
else ifeq ($(PROFILE),semi_plain64)
LOGQ = 64
M = 37
ELL = 40
BASE = 319
KMSIS = 23
M2 = 64
SIGMA1 = 27
SEED = 14
OUTER_BASE = 49152
FG_EXTRA = 2
else ifneq ($(filter $(PROFILE),partial72 plain72),)
LOGQ = 72
M = 42
ELL = 45
BASE = 319
KMSIS = 24
M2 = 74
SIGMA1 = 27
SEED = 2
OUTER_BASE = 40960
FG_EXTRA = 3
PROFILE_FLAGS = -include $(NIZK2)/profile_$(PROFILE).h
else
$(error PROFILE must be partial72, plain72, semi_partial56 or semi_plain64)
endif

COMMON_FLAGS = -std=c2x -O3 -march=native -mtune=native -fwrapv -fopenmp \
	-Wall -Wextra -Wshadow -Wundef -Wmissing-prototypes -Wredundant-decls \
	-Wno-unused-function -I$(ROOT) -I$(NIZK1) -I$(NIZK1)/src/core \
	-I$(NIZK1)/src/selector -DNIZK2_SMALL_COMMITMENT=0 \
	-DNIZK2_DOUBLE_COMMITMENT=1 -DLEOPARD_LOGQ=$(LOGQ) $(PROFILE_FLAGS)
CLIENT_FLAGS = $(COMMON_FLAGS) -DLOGQ=$(LOGQ) -DLEOPARD_ROUND_P=4 \
	-DDEMO_M=$(M) -DUMC_ELL=$(ELL) -DGADGET_BASE=$(BASE) \
	-DLNP22_KMSIS=$(KMSIS) -DLNP22_M2=$(M2) \
	-DLNP22_LOG2STDEV1=$(SIGMA1) -DLNP22_LOG2STDEV2=17 \
	-DLNP22_BATCH_ROUNDS=1 -DLNP22_CACHE_A1OPEN=0 \
	-DNIZK1_SEED_OFFSET=$(SEED) -DLNP22_DEFAULT_SEED_OFFSET=$(SEED) \
	-DLABRADOR_FG_EXTRA=$(FG_EXTRA) -DOUTER_BZ_DIGITS=2 \
	-DOUTER_BZ_DECOMP_BASE=$(OUTER_BASE) -DOUTER_XY_DECOMP_BASE=$(OUTER_BASE) \
	-DOUTER_Z2_DECOMP_BASE=2048
SERVER_FLAGS = $(COMMON_FLAGS) -DNIZK2_M=$(ELL) -DNIZK2_N=$(M) \
	-DNIZK2_PARAMS_HEADER=\"params_double_$(PROFILE)_agg.h\" \
	-DNIZK2_LEGACY_I64_INSTANCE=0 -DNIZK2_CHECK_LINEAR_RELATION=0 \
	-DNIZK2_AGGREGATE_MAX=2 -DLEOPARD_PARALLEL_BLINDEVAL=0 \
	-DLEOPARD_PARALLEL_VERIFY_REPLIES=1 -DLEOPARD_VERIFY_BUDGETED_THREADS=1 \
	-DLAZER_AFFINE_NESTED=1
SERVER_LIBS = $(ROOT)/liblazer.a -lmpfr -lgmp -lm \
	$(ROOT)/build/hexl/hexl/lib/libhexl.a -lstdc++ -fopenmp
CORE = $(addprefix $(NIZK1)/src/core/,pack.c greyhound.c dachshund.c \
	chihuahua.c labrador.c data.c jlproj.c polx.c poly.c polz.c sparsemat.c \
	ntt.S invntt.S aesctr.c gaussian_lnp.c fips202.c randombytes.c cpucycles.c)
FRONTEND = $(NIZK1)/client_inner.c $(NIZK1)/client_outer.c \
	$(wildcard $(NIZK1)/src/core/*.h $(NIZK1)/src/selector/*)
COMMON_DEPS = common.h wire.h $(NIZK2)/commitment.h .build-config

.PHONY: all check check-client clean FORCE
ifneq ($(filter $(PROFILE),semi_partial56 semi_plain64),)
all: semi-honest-bench
else
all: threshold-client threshold-server
endif

.build-config: FORCE
	@printf '%s\n' '$(CC)' '$(CLIENT_FLAGS)' '$(SERVER_FLAGS)' > $@.tmp
	@cmp -s $@.tmp $@ || mv $@.tmp $@
	@rm -f $@.tmp

threshold-client: client.c nizk1.c nizk1.h wire.c public_fastpath.h short_codec.h $(CORE) $(FRONTEND) $(COMMON_DEPS)
	$(CC) $(CLIENT_FLAGS) client.c wire.c nizk1.c $(CORE) -o $@ -lm -fopenmp

ifneq ($(filter $(PROFILE),partial72 plain72),)
$(ROOT)/build/request-round-$(PROFILE): $(ROOT)/tests/request-round-fixture.c client.c nizk1.c nizk1.h wire.c public_fastpath.h short_codec.h $(CORE) $(FRONTEND) $(COMMON_DEPS)
	@mkdir -p $(ROOT)/build
	$(CC) $(CLIENT_FLAGS) $< wire.c nizk1.c $(CORE) -Wl,--wrap=randombytes -o $@ -lm -fopenmp
endif

threshold-server: server.c nizk2.c nizk2.h wire.c $(NIZK2)/server_proof.c \
	$(NIZK2)/commitment_proof.h $(NIZK2)/params_double_$(PROFILE)_agg.h \
	$(NIZK2)/profile_$(PROFILE).h $(ROOT)/liblazer.a $(COMMON_DEPS)
	$(CC) $(SERVER_FLAGS) server.c wire.c nizk2.c $(NIZK1)/src/core/aesctr.c -o $@ $(SERVER_LIBS)

semi-honest-bench: semi-honest-bench.c nizk1.c nizk1.h wire.c public_fastpath.h short_codec.h $(CORE) $(FRONTEND) $(COMMON_DEPS)
	$(CC) $(CLIENT_FLAGS) semi-honest-bench.c wire.c $(CORE) -o $@ -lm -fopenmp

nizk1-check: $(ROOT)/tests/nizk1.c nizk1.c wire.c $(CORE) $(FRONTEND) $(COMMON_DEPS)
	$(CC) $(CLIENT_FLAGS) $(ROOT)/tests/nizk1.c wire.c $(CORE) -o $@ -lm -fopenmp

arithmetic-check: $(ROOT)/tests/arithmetic.c $(CORE) $(FRONTEND) $(COMMON_DEPS)
	$(CC) $(CLIENT_FLAGS) $(ROOT)/tests/arithmetic.c $(CORE) -o $@ -lm -fopenmp

check-client: threshold-client nizk1-check arithmetic-check
	./arithmetic-check
	./nizk1-check
	cd $(ROOT) && python3 tests/nizk1-wire.py --seeds 0,1,2 --patterns --output results/check-$(PROFILE)
	cd $(ROOT) && python3 tests/cached-request.py

nizk2-check: $(ROOT)/tests/nizk2.c nizk2.c nizk2.h wire.c \
	$(NIZK2)/server_proof.c $(NIZK2)/commitment_proof.h \
	$(NIZK2)/params_double_$(PROFILE)_agg.h $(NIZK2)/profile_$(PROFILE).h \
	$(ROOT)/liblazer.a $(COMMON_DEPS)
	$(CC) $(SERVER_FLAGS) $(ROOT)/tests/nizk2.c wire.c $(NIZK1)/src/core/aesctr.c -o $@ $(SERVER_LIBS)

check: nizk2-check
	NIZK2_L2_NEGATIVE_TESTS=1 ./nizk2-check

clean:
	rm -f threshold-client threshold-server semi-honest-bench nizk1-check nizk2-check arithmetic-check .build-config .build-config.tmp

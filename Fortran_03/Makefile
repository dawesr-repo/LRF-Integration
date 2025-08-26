# -------- Makefile --------
# Usage:
#   make test            # build (release) and run ./ejec_test.x
#   make release         # build only (O3, native, LTO)
#   make ofast           # build only (Ofast, native, LTO)
#   make test-ofast      # build Ofast and run
#   make debug           # debug build (O0, bounds checks)
#   make test-debug      # build debug and run
#   make pgo-gen         # instrumented build for PGO
#   make pgo-use         # use collected profiles
#   make clean
#
# To pass args to the program: make test ARGS="...your args..."

FC      ?= gfortran
TARGET  ?= ejec_test.x

# Explicit order to satisfy Fortran module deps:
SRCS := Fitting_Constant_v2.f90 Geometry_Constant_v2.f90 helper_functions.f90 Testing_v2.f90 test_v2.f90
OBJS := $(SRCS:.f90=.o)

ARCH    ?= -march=native -mtune=native
LTO     ?= -flto
VEC     ?= -funroll-loops -floop-nest-optimize
FPMATH  ?= -ffp-contract=fast -fno-math-errno -fno-trapping-math
WARN    ?=

ODBG    := -O0 -g -fcheck=bounds -fbacktrace
OREL    := -O3 $(ARCH) $(LTO) $(VEC) -fomit-frame-pointer $(FPMATH)
OFAST   := -Ofast $(ARCH) $(LTO) $(VEC) -fomit-frame-pointer $(FPMATH)

PGODIR  ?= pgodata
PGOGEN  := -fprofile-generate=$(PGODIR)
PGOUSE  := -fprofile-use=$(PGODIR) -fprofile-correction

.PHONY: all release debug ofast clean pgo-gen pgo-use test test-ofast test-debug run vec-report

all: release

%.o: %.f90
	$(FC) $(FFLAGS) $(WARN) -c $< -o $@

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(WARN) $(OBJS) -o $@

# Build flavors
release: FFLAGS := $(OREL)
release: $(TARGET)

debug:   FFLAGS := $(ODBG)
debug:   $(TARGET)

ofast:   FFLAGS := $(OFAST)
ofast:   $(TARGET)

# Test (build + run)
test: release
	./$(TARGET) $(ARGS)

test-ofast: ofast
	./$(TARGET) $(ARGS)

test-debug: debug
	./$(TARGET) $(ARGS)

# PGO
pgo-gen: FFLAGS := $(OREL) $(PGOGEN)
pgo-gen: $(TARGET)
	@echo "Run ./$(TARGET) with a representative workload, then 'make pgo-use'."

pgo-use: clean
pgo-use: FFLAGS := $(OREL) $(PGOUSE)
pgo-use: $(TARGET)

# Optional vectorization reports
vec-report: FFLAGS := $(OREL) -fopt-info-vec-optimized -fopt-info-vec-missed
vec-report: clean $(TARGET)
	@echo "Vectorization reports emitted (look for *.optrpt)."

clean:
	$(RM) $(OBJS) $(TARGET) *.optrpt
	@if [ -d "$(PGODIR)" ]; then find "$(PGODIR)" -type f -delete; fi

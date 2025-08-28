# ============================================================
# Makefile — LRF-Integration (Fortran) with `par` modifier
# ============================================================
# Usage:
#   make test                         # build+run SERIAL (no MPI, OpenMP runtime forced to 1 thread)
#   make par test                     # build parallel-capable, run serial (flags=0 → 1 OMP thread)
#   make par test MPI_flag=1 nproc_in=2
#   make par test OMP_flag=1 THREADS=8
#   make par test MPI_flag=1 nproc_in=2 OMP_flag=1 THREADS=8
#   make clean | make check-mpi
# Env/overrides:
#   FFLAGS_EXTRA=...  LDFLAGS_EXTRA=...  LDLIBS_EXTRA=...
# ============================================================

# ---------- primary executable under test ----------
TARGET ?= ejec_test.x
ARGS   ?=

# ---------- detect if user included `par` in goals ----------
.PHONY: par
par:
	@$(MAKE) BUILD_PAR=1 $(filter-out $@,$(MAKECMDGOALS))

# Default: serial build (no MPI/OpenMP)
BUILD_PAR ?= 0

# ---------- map user runtime vars (accept both lower/upper names) ----------
# CLI positional args expected by program: THREADS NPROC MPI_flag OMP_flag
THREADS ?= 0
ifneq ($(strip $(threads)),)
  THREADS := $(threads)
endif

NPROC   ?= 1
ifneq ($(strip $(nproc_in)),)
  NPROC := $(nproc_in)
endif

MPI_FLAG ?= 0
ifneq ($(strip $(MPI_flag)),)
  MPI_FLAG := $(MPI_flag)
endif

OMP_FLAG ?= 0
ifneq ($(strip $(OMP_flag)),)
  OMP_FLAG := $(OMP_flag)
endif

# ---------- toolchain + compile toggles ----------
FC         ?= gfortran
OPTFLAGS = -O3 -march=native -mtune=native -flto -funroll-loops -floop-nest-optimize -fomit-frame-pointer \
           -ffp-contract=fast -fno-math-errno -fno-trapping-math
FFLAGS   = $(OPTFLAGS) -cpp -std=f2008 -ffree-line-length-none
LDFLAGS  =
LDLIBS   =

# Allow easy CLI extension
FFLAGS  += $(FFLAGS_EXTRA)
LDFLAGS += $(LDFLAGS_EXTRA)
LDLIBS  += $(LDLIBS_EXTRA)

# Serial vs Parallel-capable build
ifeq ($(BUILD_PAR),1)
  # Enable OpenMP & MPI codepaths
  FFLAGS  += -fopenmp -DUSE_MPI
  LDFLAGS += -fopenmp
  # Prefer an MPI wrapper if available
  MPI_FC := $(shell command -v mpifort 2>/dev/null || command -v mpif90 2>/dev/null || command -v mpiifort 2>/dev/null)
  ifneq ($(MPI_FC),)
    FC := $(MPI_FC)
  else
    $(warning BUILD_PAR=1 but no MPI wrapper found; compiling with $(FC). Provide MPI include/libs if needed.)
  endif
else
  # Pure serial build (no MPI/OpenMP)
  # (Keep -cpp for conditional code; no -DUSE_MPI, no -fopenmp)
endif

# ---------- sources / objects ----------
# If a module version exists, use it and drop the legacy Testing_v2.f90
ifneq ($(wildcard testing_v2_mod.f90),)
  SRC_COMMON   = Fitting_Constant_v2.f90 Geometry_Constant_v2.f90 helper_functions.f90 testing_v2_mod.f90
  TEST_SRCS    = test_v2.f90
  NEEDS_MOD    = 1
else
  SRC_COMMON   = Fitting_Constant_v2.f90 Geometry_Constant_v2.f90 helper_functions.f90
  TEST_SRCS    = Testing_v2.f90 test_v2.f90
  NEEDS_MOD    = 0
endif

OBJ_COMMON   = $(SRC_COMMON:.f90=.o)
TEST_OBJS    = $(TEST_SRCS:.f90=.o)

EX_SOURCES   = min_example.f90 min_example_batch.f90
EX_OBJECTS   = $(EX_SOURCES:.f90=.o)

# Build these if present
EX_BINS      = $(foreach f,min_example min_example_batch,$(if $(wildcard $(f).f90),$(f).x,))
ALL_BINS     = $(TARGET) $(EX_BINS)

# ---------- dependencies between modules ----------
Geometry_Constant_v2.o: Fitting_Constant_v2.o
helper_functions.o:     Fitting_Constant_v2.o Geometry_Constant_v2.o
ifeq ($(NEEDS_MOD),1)
test_v2.o: testing_v2_mod.o
endif

# ---------- link rules ----------
$(TARGET): $(OBJ_COMMON) $(TEST_OBJS)
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

min_example.x: min_example.o $(OBJ_COMMON)
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

min_example_batch.x: min_example_batch.o $(OBJ_COMMON)
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

# ---------- compile pattern ----------
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# ---------- high-level targets ----------
.PHONY: all release test batch clean help check-mpi
all: release
release: $(TARGET)

# Test (build + run)
# Build may be serial or parallel-capable. We always ensure a sane OMP thread count at runtime.
test: release
	@echo ">> RUN: par=$(BUILD_PAR)  NPROC=$(NPROC) THREADS=$(THREADS) MPI_flag=$(MPI_FLAG) OMP_flag=$(OMP_FLAG)"
	@{ \
	  t="$(THREADS)"; \
	  if [ -z "$$t" ] || [ "$$t" -le 0 ]; then t=1; fi; \
	  export OMP_NUM_THREADS="$$t"; \
	  export OMP_DYNAMIC=FALSE; export OMP_PROC_BIND=close; export OMP_PLACES=cores; \
	  if [ "$(MPI_FLAG)" = "1" ] && [ "$(NPROC)" -gt 1 ]; then \
	    mpirun -np $(NPROC) ./$(TARGET) $$t $(NPROC) $(MPI_FLAG) $(OMP_FLAG) $(ARGS); \
	  else \
	    ./$(TARGET) $$t $(NPROC) $(MPI_FLAG) $(OMP_FLAG) $(ARGS); \
	  fi; \
	}

# Example batch runner (kept for convenience)
batch: min_example_batch.x
	@echo ">> Running batch: par=$(BUILD_PAR)  NPROC=$(NPROC) THREADS=$(THREADS) MPI_flag=$(MPI_FLAG) OMP_flag=$(OMP_FLAG)"
	@{ \
	  t="$(THREADS)"; \
	  if [ -z "$$t" ] || [ "$$t" -le 0 ]; then t=1; fi; \
	  export OMP_NUM_THREADS="$$t"; \
	  export OMP_DYNAMIC=FALSE; export OMP_PROC_BIND=close; export OMP_PLACES=cores; \
	  if [ "$(MPI_FLAG)" = "1" ] && [ "$(NPROC)" -gt 1 ]; then \
	    mpirun -np $(NPROC) ./min_example_batch.x $$t $(NPROC) $(MPI_FLAG) $(OMP_FLAG); \
	  else \
	    ./min_example_batch.x $$t $(NPROC) $(MPI_FLAG) $(OMP_FLAG); \
	  fi; \
	}

help:
	@echo "make test                   # build serial & run (forces OMP_NUM_THREADS>=1)"
	@echo "make par test               # build parallel-capable & run serial (flags=0)"
	@echo "make par test MPI_flag=1 nproc_in=2"
	@echo "make par test OMP_flag=1 THREADS=8"
	@echo "make par test MPI_flag=1 nproc_in=2 OMP_flag=1 THREADS=8"
	@echo "FFLAGS_EXTRA/LDFLAGS_EXTRA/LDLIBS_EXTRA accepted"
	@echo "make batch ...              # run min_example_batch.x with same flags"
	@echo "make check-mpi | make clean"

# ---------- MPI sanity check ----------
check-mpi:
	@echo "Fortran MPI wrapper:" ; \
	( command -v mpifort || command -v mpif90 || command -v mpiifort ) || { echo "  NOT FOUND"; exit 1; }
	@echo "mpirun:" ; ( mpirun --version 2>/dev/null | head -n1 || echo "  NOT FOUND" )
	@echo "mpi.mod (first hits):" ; ( find /usr /usr/local /opt -name 'mpi.mod' 2>/dev/null | head -n3 || true )
	@echo "Compile+run smoke test:" ; \
	fc=$$(command -v mpifort || command -v mpif90 || command -v mpiifort); \
	if [ -n "$$fc" ]; then \
	  printf '%s\n' \
'program hello' \
'  use mpi' \
'  implicit none' \
'  integer :: ierr, rank, size' \
'  call MPI_Init(ierr)' \
'  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)' \
'  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)' \
"  print *, 'Hello from rank', rank, 'of', size" \
'  call MPI_Finalize(ierr)' \
'end program hello' \
	  > /tmp/mpi_hello.f90; \
	  "$$fc" -ffree-line-length-none /tmp/mpi_hello.f90 -o /tmp/mpi_hello.x && mpirun -np 2 /tmp/mpi_hello.x; \
	else echo "  Skipped (no wrapper)"; fi

# ---------- clean ----------
PGODIR ?= ./pgo
clean:
	$(RM) $(OBJ_COMMON) $(TEST_OBJS) $(EX_OBJECTS) $(ALL_BINS) *.optrpt *.mod *.smod
	@if [ -d "$(PGODIR)" ]; then find "$(PGODIR)" -type f -delete; fi
	@find . \( -name '*.o' -o -name '*.mod' -o -name '*.smod' -o -name '*.optrpt' -o -name '*.x' \) -type f -delete

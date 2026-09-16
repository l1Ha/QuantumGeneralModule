# ==============================================================================
# GeneralModule - Modern Fortran Computational Physics Framework
# Makefile for compilation, testing, and maintenance
# ==============================================================================

FC := gfortran
export DEVELOPER_DIR ?= /Library/Developer/CommandLineTools
FFLAGS ?= -O2 -fPIC -Wall -Wextra -std=f2008 -ffree-line-length-none
AR ?= ar
ARFLAGS ?= rcs

SRC_DIR = src
TEST_DIR = tests
EX_DIR = examples
LIB_DIR = lib
MOD_DIR = include

# 48 核心模块源文件 (按拓扑依赖顺序排列)
MODULE_SRCS = \
	$(SRC_DIR)/mod_constants.f90 \
	$(SRC_DIR)/mod_special_functions.f90 \
	$(SRC_DIR)/mod_linear_algebra.f90 \
	$(SRC_DIR)/mod_dvr_grid.f90 \
	$(SRC_DIR)/mod_laser_pulse.f90 \
	$(SRC_DIR)/mod_absorbing_boundary.f90 \
	$(SRC_DIR)/mod_thermal_ensemble.f90 \
	$(SRC_DIR)/mod_wavepacket_propagator.f90 \
	$(SRC_DIR)/mod_coulomb_atomic.f90 \
	$(SRC_DIR)/mod_hhg_spectra.f90 \
	$(SRC_DIR)/mod_chebyshev_propagator.f90 \
	$(SRC_DIR)/mod_multistate_coupling.f90 \
	$(SRC_DIR)/mod_rovibrational.f90 \
	$(SRC_DIR)/mod_io_utils.f90 \
	$(SRC_DIR)/mod_interpolation.f90 \
	$(SRC_DIR)/mod_photofragment_flux.f90 \
	$(SRC_DIR)/mod_open_quantum.f90 \
	$(SRC_DIR)/mod_optimal_control.f90 \
	$(SRC_DIR)/mod_ti_scattering.f90 \
	$(SRC_DIR)/mod_td_scattering.f90 \
	$(SRC_DIR)/mod_field_scattering.f90 \
	$(SRC_DIR)/mod_dipolar_scattering.f90 \
	$(SRC_DIR)/mod_photoassociation.f90 \
	$(SRC_DIR)/mod_three_body_recombination.f90 \
	$(SRC_DIR)/mod_confined_scattering.f90 \
	$(SRC_DIR)/mod_autoionization_fano.f90 \
	$(SRC_DIR)/mod_crossed_field_scattering.f90 \
	$(SRC_DIR)/mod_triatomic_geometry.f90 \
	$(SRC_DIR)/mod_spinor_bec.f90 \
	$(SRC_DIR)/mod_hyperspherical_reactive.f90 \
	$(SRC_DIR)/mod_dipolar_droplets_lhy.f90 \
	$(SRC_DIR)/mod_strong_field_nsdi.f90 \
	$(SRC_DIR)/mod_feshbach_bound_states.f90 \
	$(SRC_DIR)/mod_attosecond_transient_absorption.f90 \
	$(SRC_DIR)/mod_bicircular_pecd.f90 \
	$(SRC_DIR)/mod_ultracold_reaction_shielding.f90 \
	$(SRC_DIR)/mod_rydberg_blockade.f90 \
	$(SRC_DIR)/mod_surface_scattering.f90 \
	$(SRC_DIR)/mod_surface_reaction_er.f90 \
	$(SRC_DIR)/mod_surface_electronic_friction.f90 \
	$(SRC_DIR)/mod_grazing_fast_atom_diffraction.f90 \
	$(SRC_DIR)/mod_ion_atom_scattering.f90 \
	$(SRC_DIR)/mod_surface_hopping_fssh.f90 \
	$(SRC_DIR)/mod_molecular_alignment.f90 \
	$(SRC_DIR)/mod_optical_lattice_hubbard.f90 \
	$(SRC_DIR)/mod_reaction_path_hamiltonian.f90 \
	$(SRC_DIR)/mod_relativistic_atomic.f90 \
	$(SRC_DIR)/mod_resonant_xray_scattering.f90 \
	$(SRC_DIR)/general_module.f90

MODULE_OBJS = $(MODULE_SRCS:.f90=.o)
STATIC_LIB = $(SRC_DIR)/libgeneral_module.a

.PHONY: all lib test examples clean check help

all: lib

lib: $(STATIC_LIB)

$(STATIC_LIB): $(MODULE_SRCS)
	@echo "==> Building GeneralModule static library..."
	@for src in $(MODULE_SRCS); do \
		echo "    Compiling $$src..."; \
		$(FC) $(FFLAGS) -I$(SRC_DIR) -J$(SRC_DIR) -c $$src -o $${src%.f90}.o || exit 1; \
	done
	@$(AR) $(ARFLAGS) $(STATIC_LIB) $(SRC_DIR)/*.o
	@echo "==> Library built: $(STATIC_LIB)"

test: lib
	@echo "==> Running full unit test suite (40 test suites)..."
	@bash $(TEST_DIR)/run_all_tests.sh

examples: lib
	@echo "==> Building and running physical examples (35 examples)..."
	@bash $(EX_DIR)/build_examples.sh

check: test
	@echo "==> Running Python tests..."
	@python3 python/test_pygenmod.py

clean:
	@echo "==> Cleaning build artifacts..."
	@rm -f $(SRC_DIR)/*.o $(SRC_DIR)/*.mod $(SRC_DIR)/*.a
	@rm -f $(TEST_DIR)/*.o $(TEST_DIR)/*.mod $(TEST_DIR)/*.dat
	@rm -f $(EX_DIR)/*.o $(EX_DIR)/*.mod $(EX_DIR)/*.dat
	@rm -f *.o *.mod *.a *.dylib *.so *.dat *.png
	@rm -rf build .fpm
	@echo "==> Clean complete."

help:
	@echo "Available targets in GeneralModule Makefile:"
	@echo "  make all      - Build the core static library (libgeneral_module.a)"
	@echo "  make test     - Run all 40 unit test suites (335 assertions)"
	@echo "  make examples - Build and execute all 35 physics examples"
	@echo "  make check    - Run Fortran tests and Python pygenmod test suite"
	@echo "  make clean    - Remove all compiled objects, modules, and data outputs"

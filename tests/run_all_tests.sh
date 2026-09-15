#!/bin/bash
# ==============================================================================
# GeneralModule 自动化测试套件运行脚本
# ==============================================================================
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SRC_DIR="$DIR/../src"

if [ -d "/Library/Developer/CommandLineTools" ]; then
    export DEVELOPER_DIR="/Library/Developer/CommandLineTools"
fi

echo "================================================================"
echo "          Running GeneralModule Test Suite Suite                "
echo "================================================================"

# 1. 编译 src 中的模块 (若尚未编译)
echo "[1/2] Compiling all modules in src/..."
cd "$SRC_DIR"
FFLAGS="-O2 -fPIC -ffree-line-length-none"
gfortran $FFLAGS -c mod_constants.f90
gfortran $FFLAGS -c mod_special_functions.f90
gfortran $FFLAGS -c mod_linear_algebra.f90
gfortran $FFLAGS -c mod_dvr_grid.f90
gfortran $FFLAGS -c mod_laser_pulse.f90
gfortran $FFLAGS -c mod_absorbing_boundary.f90
gfortran $FFLAGS -c mod_thermal_ensemble.f90
gfortran $FFLAGS -c mod_wavepacket_propagator.f90
gfortran $FFLAGS -c mod_coulomb_atomic.f90
gfortran $FFLAGS -c mod_hhg_spectra.f90
gfortran $FFLAGS -c mod_chebyshev_propagator.f90
gfortran $FFLAGS -c mod_multistate_coupling.f90
gfortran $FFLAGS -c mod_rovibrational.f90
gfortran $FFLAGS -c mod_io_utils.f90
gfortran $FFLAGS -c mod_interpolation.f90
gfortran $FFLAGS -c mod_photofragment_flux.f90
gfortran $FFLAGS -c mod_open_quantum.f90
gfortran $FFLAGS -c mod_optimal_control.f90
gfortran $FFLAGS -c mod_ti_scattering.f90
gfortran $FFLAGS -c mod_td_scattering.f90
gfortran $FFLAGS -c mod_field_scattering.f90
gfortran $FFLAGS -c mod_dipolar_scattering.f90
gfortran $FFLAGS -c mod_photoassociation.f90
gfortran $FFLAGS -c mod_three_body_recombination.f90
gfortran $FFLAGS -c mod_confined_scattering.f90
gfortran $FFLAGS -c mod_autoionization_fano.f90
gfortran $FFLAGS -c mod_crossed_field_scattering.f90
gfortran $FFLAGS -c mod_triatomic_geometry.f90
gfortran $FFLAGS -c mod_spinor_bec.f90
gfortran $FFLAGS -c mod_hyperspherical_reactive.f90
gfortran $FFLAGS -c mod_dipolar_droplets_lhy.f90
gfortran $FFLAGS -c mod_strong_field_nsdi.f90
gfortran $FFLAGS -c mod_feshbach_bound_states.f90
gfortran $FFLAGS -c mod_attosecond_transient_absorption.f90
gfortran $FFLAGS -c mod_bicircular_pecd.f90
gfortran $FFLAGS -c mod_ultracold_reaction_shielding.f90
gfortran $FFLAGS -c mod_rydberg_blockade.f90
gfortran $FFLAGS -c general_module.f90
ar rcs libgeneral_module.a *.o

cd "$DIR"

# 2. 编译并运行各项测试 (共 29 大完整测试套件)
TESTS=(
    "test_constants"
    "test_special_functions"
    "test_dvr_grid"
    "test_laser_pulse"
    "test_propagators"
    "test_atomic_hhg"
    "test_laser_rovibrational_control"
    "test_interpolation"
    "test_photofragment_flux"
    "test_open_quantum_opt"
    "test_ti_scattering"
    "test_td_scattering"
    "test_field_scattering"
    "test_dipolar_scattering"
    "test_photoassociation"
    "test_three_body_recombination"
    "test_confined_scattering"
    "test_autoionization_fano"
    "test_crossed_field_scattering"
    "test_triatomic_geometry"
    "test_spinor_bec"
    "test_hyperspherical_reactive"
    "test_dipolar_droplets_lhy"
    "test_strong_field_nsdi"
    "test_feshbach_bound_states"
    "test_attosecond_transient_absorption"
    "test_bicircular_pecd"
    "test_ultracold_reaction_shielding"
    "test_rydberg_blockade"
)

for test_name in "${TESTS[@]}"; do
    echo ""
    echo ">> Compiling and running $test_name..."
    gfortran $FFLAGS -I"$SRC_DIR" "$test_name.f90" "$SRC_DIR"/libgeneral_module.a -o "$test_name"
    "./$test_name"
    rm -f "$test_name"
done

echo ""
echo "================================================================"
echo "       ALL UNIT TESTS PASSED SUCCESSFULLY! (100% Pass)          "
echo "================================================================"

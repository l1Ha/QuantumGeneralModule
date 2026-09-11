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
gfortran -O2 -fPIC -c mod_constants.f90
gfortran -O2 -fPIC -c mod_special_functions.f90
gfortran -O2 -fPIC -c mod_linear_algebra.f90
gfortran -O2 -fPIC -c mod_dvr_grid.f90
gfortran -O2 -fPIC -c mod_laser_pulse.f90
gfortran -O2 -fPIC -c mod_absorbing_boundary.f90
gfortran -O2 -fPIC -c mod_thermal_ensemble.f90
gfortran -O2 -fPIC -c mod_wavepacket_propagator.f90
gfortran -O2 -fPIC -c mod_coulomb_atomic.f90
gfortran -O2 -fPIC -c mod_hhg_spectra.f90
gfortran -O2 -fPIC -c mod_chebyshev_propagator.f90
gfortran -O2 -fPIC -c mod_multistate_coupling.f90
gfortran -O2 -fPIC -c mod_rovibrational.f90
gfortran -O2 -fPIC -c mod_io_utils.f90
gfortran -O2 -fPIC -c mod_interpolation.f90
gfortran -O2 -fPIC -c mod_photofragment_flux.f90
gfortran -O2 -fPIC -c mod_open_quantum.f90
gfortran -O2 -fPIC -c mod_optimal_control.f90
gfortran -O2 -fPIC -c general_module.f90
ar rcs libgeneral_module.a *.o

cd "$DIR"

# 2. 编译并运行各项测试 (共 10 大完整测试套件)
TESTS=("test_constants" "test_special_functions" "test_dvr_grid" "test_laser_pulse" "test_propagators" "test_atomic_hhg" "test_laser_rovibrational_control" "test_interpolation" "test_photofragment_flux" "test_open_quantum_opt")

for test_name in "${TESTS[@]}"; do
    echo ""
    echo ">> Compiling and running $test_name..."
    gfortran -O2 -I"$SRC_DIR" "$test_name.f90" "$SRC_DIR"/libgeneral_module.a -o "$test_name"
    "./$test_name"
    rm -f "$test_name"
done

echo ""
echo "================================================================"
echo "       ALL UNIT TESTS PASSED SUCCESSFULLY! (100% Pass)          "
echo "================================================================"

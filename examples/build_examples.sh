#!/bin/bash
# ==============================================================================
# GeneralModule 示例编译与运行脚本
# ==============================================================================
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SRC_DIR="$DIR/../src"

# 仅在 macOS 且确实安装了 Xcode 命令行工具时设置 DEVELOPER_DIR；
# 其他平台无条件导出会污染子进程环境。
if [ -d "/Library/Developer/CommandLineTools" ]; then
    export DEVELOPER_DIR="/Library/Developer/CommandLineTools"
fi

echo "================================================================"
echo "          Building and Running GeneralModule Examples           "
echo "================================================================"

cd "$DIR"

# 确保静态库与本次构建标志一致: 避免与 OPENMP=1 构建的库混链
# (OpenMP 目标文件引用 GOMP_* 符号, 无 -fopenmp 链接会未定义)
NJOBS="${TEST_JOBS:-$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null) || echo 1)}"
FFLAGS="-O2 -fPIC -ffree-line-length-none"
if [ "${OPENMP:-0}" = "1" ]; then
    FFLAGS="$FFLAGS -fopenmp"
fi
export FFLAGS
echo "[0/1] Rebuilding module archive with consistent flags..."
FFLAGS="$FFLAGS" bash "$DIR/../scripts/compile_modules.sh" "$SRC_DIR" "$NJOBS"
ar rcs "$SRC_DIR"/libgeneral_module.a "$SRC_DIR"/*.o

EXAMPLES=(
    "ex01_fgh_diatomic_bound_states"
    "ex02_pulse_synthesis"
    "ex03_split_operator_1d"
    "ex04_field_free_orientation"
    "ex05_hhg_lewenstein_spectrum"
    "ex06_two_state_nonadiabatic"
    "ex07_scattering_wavefunctions_ti_td"
    "ex08_ultracold_feshbach_segmented"
    "ex09_dipolar_relaxation_scattering"
    "ex10_photoassociation_spectroscopy"
    "ex11_three_body_efimov_recombination"
    "ex12_confined_cir_scattering"
    "ex13_autoionization_fano_resonance"
    "ex14_spinor_bec_dynamics"
    "ex15_crossed_field_stark_zeeman"
    "ex16_triatomic_reaction_berry_phase"
    "ex17_hyperspherical_reaction_rates"
    "ex18_dipolar_quantum_droplets"
    "ex19_strong_field_nsdi_recollision"
    "ex20_feshbach_molecular_bound_states"
    "ex21_attosecond_transient_absorption"
    "ex22_bicircular_pecd_chiral"
    "ex23_ultracold_molecule_shielding"
    "ex24_rydberg_blockade_dynamics"
    "ex25_surface_corrugated_diffraction"
    "ex26_eley_rideal_surface_reaction"
    "ex27_surface_electronic_friction_gle"
    "ex28_grazing_fast_atom_diffraction"
    "ex29_cold_ion_atom_scattering"
    "ex30_tully_surface_hopping"
    "ex31_molecular_alignment_revival"
    "ex32_optical_lattice_bose_hubbard"
    "ex33_rph_variational_transition_state"
    "ex34_relativistic_dirac_cesium"
    "ex35_rixs_core_level_spectroscopy"
    "ex36_penning_associative_ionization"
)

for ex in "${EXAMPLES[@]}"; do
    echo ""
    echo ">> Building & Running $ex..."
    gfortran $FFLAGS -I"$SRC_DIR" "$ex.f90" "$SRC_DIR"/libgeneral_module.a -o "$ex"
    "./$ex"
    rm -f "$ex"
done

echo ""
echo "================================================================"
echo "          ALL EXAMPLES EXECUTED SUCCESSFULLY!                   "
echo "================================================================"

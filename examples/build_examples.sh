#!/bin/bash
# ==============================================================================
# GeneralModule 示例编译与运行脚本
# ==============================================================================
set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SRC_DIR="$DIR/../src"
export DEVELOPER_DIR="/Library/Developer/CommandLineTools"

echo "================================================================"
echo "          Building and Running GeneralModule Examples           "
echo "================================================================"

cd "$DIR"

FFLAGS="-O2 -fPIC -ffree-line-length-none"
EXAMPLES=("ex01_fgh_diatomic_bound_states" "ex02_pulse_synthesis" "ex03_split_operator_1d" "ex04_field_free_orientation" "ex05_hhg_lewenstein_spectrum" "ex06_two_state_nonadiabatic" "ex07_scattering_wavefunctions_ti_td" "ex08_ultracold_feshbach_segmented" "ex09_dipolar_relaxation_scattering" "ex10_photoassociation_spectroscopy")

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

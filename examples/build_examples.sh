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

EXAMPLES=("ex01_fgh_diatomic_bound_states" "ex02_pulse_synthesis" "ex03_split_operator_1d" "ex04_field_free_orientation" "ex05_hhg_lewenstein_spectrum" "ex06_two_state_nonadiabatic")

for ex in "${EXAMPLES[@]}"; do
    echo ""
    echo ">> Building & Running $ex..."
    gfortran -O2 -I"$SRC_DIR" "$ex.f90" "$SRC_DIR"/libgeneral_module.a -o "$ex"
    "./$ex"
    rm -f "$ex"
done

echo ""
echo "================================================================"
echo "          ALL EXAMPLES EXECUTED SUCCESSFULLY!                   "
echo "================================================================"

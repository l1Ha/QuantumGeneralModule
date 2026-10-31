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

# 远程诊断: CI 中通过公开注解标记阶段与失败行号 (无需登录读日志)
if [ -n "${GITHUB_ACTIONS:-}" ]; then
    trap 'rc=$?; echo "::error::run_all_tests.sh aborted at line $LINENO (exit $rc)"' ERR
fi

# ---------------------------------------------------------------------------
# worker 模式 (必须先于主流程): xargs 通过 `bash "$SELF" --run-one <name>`
# 自调用执行单个套件。不使用 export -f / bash -c —— macOS /bin/bash 3.2
# 的函数导出行为不可靠, 曾导致并行子进程整体静默失败。
# ---------------------------------------------------------------------------
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SRC_DIR="$DIR/../src"

# worker 自调用所需的脚本绝对路径: 必须在任何 cd 之前、用 $PWD 解析
# (脚本中的 `cd "$DIR"` 会使相对路径 $0 失效)
SELF="${BASH_SOURCE[0]}"
case "$SELF" in
    /*) ;;
    *) SELF="$PWD/$SELF" ;;
esac

# 可移植纳秒时间戳: GNU date 支持 %N, BSD/macOS date 不支持 (回退到秒)
t_stamp() {
    local s
    s="$(date +%s%N 2>/dev/null)"
    case "$s" in
        ''|*[!0-9]*) date +%s ;;
        *) echo "$s" ;;
    esac
}

run_one_test() {
    local name="$1"
    local rc=0
    local t0 t1
    t0="$(t_stamp)"
    rm -f "$DIR/$name"
    gfortran $FFLAGS -I"$SRC_DIR" "$DIR/$name.f90" "$SRC_DIR/libgeneral_module.a" -o "$DIR/$name" \
        > "$LOG_DIR/$name.build.log" 2>&1 || rc=1
    if [ "$rc" -eq 0 ]; then
        "$DIR/$name" > "$LOG_DIR/$name.log" 2>&1 || rc=1
    fi
    t1="$(t_stamp)"
    if [ -n "$t0" ] && [ -n "$t1" ] && [ "$t1" -ge "$t0" ] 2>/dev/null; then
        echo "$(( (t1 - t0) / 1000000 ))" > "$LOG_DIR/$name.ms"
    fi
    rm -f "$DIR/$name"
    if [ "$rc" -eq 0 ]; then
        echo "PASS" > "$LOG_DIR/$name.status"
    else
        echo "FAIL" > "$LOG_DIR/$name.status"
    fi
    return 0
}

if [ "${1:-}" = "--run-one" ]; then
    run_one_test "$2"
    exit 0
fi

# 1. 编译 src 中的模块 (依赖自愈的并行编译; TEST_JOBS 可覆盖并行度)
# 透传 OPENMP 标志: make OPENMP=1 test 时库与测试均以 -fopenmp 编译
NJOBS="${TEST_JOBS:-$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null) || echo 1)}"
echo "[1/2] Compiling all modules in src/ (parallel, up to $NJOBS jobs)..."
FFLAGS="-O2 -fPIC -ffree-line-length-none"
if [ "${OPENMP:-0}" = "1" ]; then
    FFLAGS="$FFLAGS -fopenmp"
fi
export FFLAGS
FFLAGS="$FFLAGS" bash "$DIR/../scripts/compile_modules.sh" "$SRC_DIR" "$NJOBS"
ar rcs "$SRC_DIR"/libgeneral_module.a "$SRC_DIR"/*.o
if [ -n "${GITHUB_ACTIONS:-}" ]; then
    echo "::notice::run_all_tests.sh: module build done (bash $BASH_VERSION, jobs=$NJOBS)"
fi

cd "$DIR"

# 2. 编译并运行各项测试 (共 46 大完整测试套件)
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
    "test_surface_scattering"
    "test_surface_reaction_er"
    "test_surface_electronic_friction"
    "test_grazing_fast_atom_diffraction"
    "test_ion_atom_scattering"
    "test_surface_hopping_fssh"
    "test_molecular_alignment"
    "test_optical_lattice_hubbard"
    "test_reaction_path_hamiltonian"
    "test_relativistic_atomic"
    "test_resonant_xray_scattering"
    "test_penning_associative_ionization"
    "test_coulomb_threebody"
    "test_qct_dynamics"
    "test_sop_hamiltonian"
    "test_tensor_train"
    "test_sparse_grid"
)

# 并行执行: 各测试只打印 stdout、无共享文件 I/O, 天然可安全并行。
# 每套件输出写入独立日志, 全部完成后按固定顺序回放, 保证输出确定性;
# 失败检测不依赖 xargs 退出码 (避免 set -e 提前中断汇总)。
echo ""
echo ">> Compiling and running ${#TESTS[@]} test suites with $NJOBS parallel jobs..."
LOG_DIR="$(mktemp -d)"
export FFLAGS SRC_DIR DIR LOG_DIR

printf '%s\n' "${TESTS[@]}" | xargs -P "$NJOBS" -I '{}' bash "$SELF" --run-one '{}' || true
if [ -n "${GITHUB_ACTIONS:-}" ]; then
    echo "::notice::run_all_tests.sh: parallel test phase done"
fi

n_fail=0
for test_name in "${TESTS[@]}"; do
    echo ""
    echo ">> $test_name ($(cat "$LOG_DIR/$test_name.ms" 2>/dev/null || echo '?') ms)"
    cat "$LOG_DIR/$test_name.log" 2>/dev/null || true
    if [ "$(cat "$LOG_DIR/$test_name.status" 2>/dev/null)" != "PASS" ]; then
        n_fail=$((n_fail + 1))
        echo "[BUILD/RUN FAILURE] $test_name"
        cat "$LOG_DIR/$test_name.build.log" 2>/dev/null
        # 在 CI 中把失败详情转发为公开 workflow 注解 (无需登录读日志)
        if [ -n "${GITHUB_ACTIONS:-}" ]; then
            detail="$(tail -n 15 "$LOG_DIR/$test_name.log" 2>/dev/null; echo '--- build ---'; tail -n 12 "$LOG_DIR/$test_name.build.log" 2>/dev/null)"
            detail="${detail//\%/\%25}"
            detail="${detail//$'\n'/%0A}"
            echo "::error::SUITE FAILED: $test_name :: ${detail}"
        fi
    fi
done
rm -rf "$LOG_DIR"

if [ "$n_fail" -ne 0 ]; then
    echo ""
    echo "****************************************************************"
    echo "   $n_fail TEST SUITE(S) FAILED - see [BUILD/RUN FAILURE] above "
    echo "****************************************************************"
    exit 1
fi

echo ""
echo "================================================================"
echo "       ALL UNIT TESTS PASSED SUCCESSFULLY! (100% Pass)          "
echo "================================================================"

#!/bin/bash
# ==============================================================================
# compile_modules.sh - 依赖自愈的并行模块编译器 (GeneralModule)
#
# 用法: compile_modules.sh <src_dir> [njobs]
#   FC/FFLAGS 环境变量可覆盖编译器与编译选项 (默认 gfortran / -O2 -fPIC)。
#   通过环境变量 TEST_JOBS / 第二参数控制并行度, 默认取 nproc。
#
# 原理: 不硬编码模块依赖顺序。每轮并行编译全部未完成文件, 依赖未就绪的
# 文件会因缺少 .mod 而快速失败, 下一轮重试; 依赖链深度约 5 层, 因此通常
# 5 个波次内全部完成。若某一轮无任何进展 (存在真实语法错误等), 则对剩余
# 文件逐一详细编译并如实报错退出, 不会静默吞掉失败。
# ==============================================================================
set -u

SRC_DIR="${1:?usage: compile_modules.sh <src_dir> [njobs]}"
if [ "${2:-}" != "" ]; then
    NJOBS="$2"
else
    NJOBS="$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null) || echo 1)"
fi
FC="${FC:-gfortran}"
FFLAGS="${FFLAGS:--O2 -fPIC}"

export FC FFLAGS SRC_DIR

compile_one() {
    # shellcheck disable=SC2086
    "$FC" $FFLAGS -I"$SRC_DIR" -J"$SRC_DIR" -c "$1" -o "${1%.f90}.o"
}
export -f compile_one

# 清理旧产物, 使"文件是否存在 .o"不能误判本轮编译成功
rm -f "$SRC_DIR"/*.o "$SRC_DIR"/*.mod

# 收集全部模块源 (glob 顺序无关, 由重试机制自然收敛到拓扑序)
pending=()
for f in "$SRC_DIR"/mod_*.f90 "$SRC_DIR"/general_module.f90; do
    [ -e "$f" ] && pending+=("$f")
done
n_total="${#pending[@]}"

MARK="$(mktemp -d)"
trap 'rm -rf "$MARK"' EXIT

wave=0
while [ "${#pending[@]}" -gt 0 ]; do
    wave=$((wave + 1))
    rm -f "$MARK"/*.ok
    export MARK
    # shellcheck disable=SC2016
    printf '%s\n' "${pending[@]}" | \
        xargs -P "$NJOBS" -I '{}' bash -c '
            if compile_one "{}" >/dev/null 2>&1; then
                : > "$MARK/$(basename "{}").ok"
            fi' || true

    next=()
    for f in "${pending[@]}"; do
        [ -f "$MARK/$(basename "$f").ok" ] || next+=("$f")
    done

    if [ "${#next[@]}" -eq "${#pending[@]}" ]; then
        # 无进展: 逐文件详细编译, 真实错误直接呈现
        echo "[compile_modules] no progress in wave $wave; compiling remaining files verbosely..." >&2
        for f in "${next[@]}"; do
            compile_one "$f" || {
                echo "[compile_modules] FAILED: $f" >&2
                exit 1
            }
        done
        next=()
    fi
    pending=("${next[@]}")
done

echo "[compile_modules] compiled $n_total modules in $wave wave(s) with up to $NJOBS parallel jobs"

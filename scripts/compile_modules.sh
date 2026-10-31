#!/bin/bash
# ==============================================================================
# compile_modules.sh - 依赖自愈的并行模块编译器 (GeneralModule)
#
# 用法: compile_modules.sh <src_dir> [njobs]
#   FC/FFLAGS 环境变量可覆盖编译器与编译选项 (默认 gfortran / -O2 -fPIC)。
#   并行度取第二参数, 否则 nproc/sysctl 自动侦测。
#
# 原理: 不硬编码模块依赖顺序。每轮并行编译全部未完成文件, 依赖未就绪的
# 文件会因缺少 .mod 而快速失败, 下一轮重试; 依赖链深度约 5 层, 因此通常
# 3-5 个波次内全部完成。若某一轮无任何进展 (存在真实编译错误), 则对剩余
# 文件逐一详细编译、打印编译器原始错误并退出非零, 不会静默吞掉失败。
#
# 可移植性: 并行 worker 通过"脚本自调用"(--compile-one)实现, 不依赖
# bash 的 export -f 函数导出与 bash -c 语义 (macOS /bin/bash 3.2 上函数
# 导出行为不可靠, 曾导致全部子进程静默失败)。
# ==============================================================================
set -u

FC="${FC:-gfortran}"
FFLAGS="${FFLAGS:--O2 -fPIC}"

# ---- 单文件 worker 模式 (由 xargs 调用; 必须先于任何其他逻辑) ----------------
if [ "${1:-}" = "--compile-one" ]; then
    SRC_DIR_W="${2:?}"; f_W="${3:?}"; mark_W="${4:-}"
    # shellcheck disable=SC2086
    if "$FC" $FFLAGS -I"$SRC_DIR_W" -J"$SRC_DIR_W" -c "$f_W" -o "${f_W%.f90}.o"; then
        if [ -n "$mark_W" ]; then
            : > "$mark_W/$(basename "$f_W").ok"
        fi
        exit 0
    fi
    exit 1
fi

if [ -n "${GITHUB_ACTIONS:-}" ]; then
    trap 'rc=$?; echo "::error::compile_modules.sh aborted at line $LINENO (exit $rc)"' ERR
    echo "::notice::compile_modules.sh start: bash $BASH_VERSION"
fi

SRC_DIR="${1:?usage: compile_modules.sh <src_dir> [njobs]}"
if [ "${2:-}" != "" ]; then
    NJOBS="$2"
else
    NJOBS="$( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null) || echo 1)"
fi
SELF="$0"

# 清理旧产物, 健康起点: 成功文件必然在本脚本运行期间生成
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
    printf '%s\n' "${pending[@]}" | \
        xargs -P "$NJOBS" -I '{}' bash "$SELF" --compile-one "$SRC_DIR" '{}' "$MARK" \
        >/dev/null 2>&1 || true

    next=()
    n_next=0
    n_pending=0
    for f in "${pending[@]}"; do
        n_pending=$((n_pending + 1))
        if [ ! -f "$MARK/$(basename "$f").ok" ]; then
            next+=("$f")
            n_next=$((n_next + 1))
        fi
    done

    if [ "$n_next" -eq "$n_pending" ]; then
        # 无进展: 逐文件详细编译, 编译器原始错误直接呈现并转发为公开注解
        echo "[compile_modules] no progress in wave $wave; compiling remaining files verbosely..." >&2
        for f in "${next[@]}"; do
            if ! cerr="$("$FC" $FFLAGS -I"$SRC_DIR" -J"$SRC_DIR" -c "$f" -o "${f%.f90}.o" 2>&1)"; then
                echo "[compile_modules] FAILED: $f" >&2
                echo "$cerr" >&2
                if [ -n "${GITHUB_ACTIONS:-}" ]; then
                    cerr_esc="${cerr//\%/\%25}"
                    cerr_esc="${cerr_esc//$'\n'/%0A}"
                    echo "::error::compile_modules FAILED: $(basename "$f") :: ${cerr_esc}"
                fi
                exit 1
            fi
        done
        n_next=0
    fi
    # 注意: 不用 pending=("${next[@]}"): 空数组值展开在 set -u + bash<4.4
    # (macOS /bin/bash 3.2) 下报 unbound variable, 曾导致 CI 在最后一波死亡
    if [ "$n_next" -gt 0 ]; then
        pending=("${next[@]}")
    else
        pending=()
    fi
done

echo "[compile_modules] compiled $n_total modules in $wave wave(s) with up to $NJOBS parallel jobs"

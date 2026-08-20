#!/usr/bin/env bash

# Shared module logging and lightweight resource tracking for METAHICT modules.
# Each module writes:
#   - module.log: stdout/stderr from the module run
#   - resources.txt: peak memory and elapsed wall time

metahict_detect_outdir() {
    local script_name
    script_name="$(basename "$0")"

    if [[ "$script_name" == "6_binning.sh" && $# -ge 3 ]]; then
        printf '%s\n' "$3"
        return 0
    fi

    while [[ $# -gt 0 ]]; do
        case "$1" in
            -o|--output|--out|--outdir)
                if [[ $# -ge 2 ]]; then
                    printf '%s\n' "$2"
                    return 0
                fi
                ;;
            --output=*|--out=*|--outdir=*)
                printf '%s\n' "${1#*=}"
                return 0
                ;;
            --)
                break
                ;;
        esac
        shift
    done

    return 1
}

_metahict_resource_descendants() {
    local pid="$1"
    local child
    printf '%s\n' "$pid"
    while read -r child; do
        [[ -n "$child" ]] || continue
        _metahict_resource_descendants "$child"
    done < <(pgrep -P "$pid" 2>/dev/null || true)
}

_metahict_resource_monitor() {
    local root_pid="$1"
    local peak_file="$2"
    local interval="${METAHICT_RESOURCE_INTERVAL:-5}"
    local pids rss cgroup_current_kb

    printf '0\n' > "$peak_file"
    while kill -0 "$root_pid" 2>/dev/null; do
        cgroup_current_kb=0
        if [[ -n "${METAHICT_RESOURCE_CGROUP_MEMORY_CURRENT:-}" && -r "${METAHICT_RESOURCE_CGROUP_MEMORY_CURRENT}" ]]; then
            cgroup_current_kb="$(awk '{printf "%d", $1 / 1024}' "${METAHICT_RESOURCE_CGROUP_MEMORY_CURRENT}" 2>/dev/null || printf '0')"
        fi

        if [[ "${cgroup_current_kb}" =~ ^[0-9]+$ ]] && (( cgroup_current_kb > 0 )); then
            rss="${cgroup_current_kb}"
        else
            pids="$(_metahict_resource_descendants "$root_pid" | sort -n -u | paste -sd, -)"
            if [[ -n "$pids" ]]; then
                rss="$(ps -o rss= -p "$pids" 2>/dev/null | awk 'BEGIN {max = 0} {if ($1 > max) max = $1} END {print max + 0}')"
            else
                rss=0
            fi
        fi

        if [[ "${rss}" =~ ^[0-9]+$ ]]; then
            awk -v rss="$rss" 'NR == 1 {if (rss > $1) print rss; else print $1}' "$peak_file" > "${peak_file}.tmp"
            mv "${peak_file}.tmp" "$peak_file"
        fi
        sleep "$interval"
    done
}

_metahict_format_hms() {
    local seconds="$1"
    printf '%02d:%02d:%02d\n' "$((seconds / 3600))" "$(((seconds % 3600) / 60))" "$((seconds % 60))"
}

_metahict_resource_finish() {
    local status="$?"
    local end_time elapsed peak_kb peak_mb peak_gb elapsed_hms

    trap - EXIT

    if [[ -n "${METAHICT_RESOURCE_MONITOR_PID:-}" ]]; then
        kill "$METAHICT_RESOURCE_MONITOR_PID" 2>/dev/null || true
        wait "$METAHICT_RESOURCE_MONITOR_PID" 2>/dev/null || true
    fi

    end_time="$(date +%s)"
    elapsed=$((end_time - METAHICT_RESOURCE_START_EPOCH))
    elapsed_hms="$(_metahict_format_hms "$elapsed")"
    peak_kb="$(cat "$METAHICT_RESOURCE_PEAK_FILE" 2>/dev/null || printf '0')"
    peak_mb="$(awk -v kb="$peak_kb" 'BEGIN {printf "%.2f", kb / 1024}')"
    peak_gb="$(awk -v kb="$peak_kb" 'BEGIN {printf "%.2f", kb / 1024 / 1024}')"

    {
        printf 'highest_memory\t%s KB (%s MB, %s GB)\n' "$peak_kb" "$peak_mb" "$peak_gb"
        printf 'elapsed_time\t%s\n' "$elapsed_hms"
    } > "$METAHICT_RESOURCE_FILE"

    rm -f "$METAHICT_RESOURCE_PEAK_FILE" "${METAHICT_RESOURCE_PEAK_FILE}.tmp"
    echo "[INFO] Resource report written to: $METAHICT_RESOURCE_FILE"
    echo "[INFO] Module log written to: $METAHICT_LOG_FILE"

    exit "$status"
}

metahict_resource_start() {
    local outdir="$1"
    local module_name="${2:-$(basename "$0")}"

    [[ -n "$outdir" ]] || return 0
    mkdir -p "$outdir"

    export METAHICT_RESOURCE_START_EPOCH
    export METAHICT_RESOURCE_FILE
    export METAHICT_RESOURCE_PEAK_FILE
    export METAHICT_LOG_FILE
    export METAHICT_RESOURCE_MONITOR_PID
    export METAHICT_RESOURCE_CGROUP_MEMORY_CURRENT

    METAHICT_RESOURCE_START_EPOCH="$(date +%s)"
    METAHICT_RESOURCE_FILE="${outdir}/resources.txt"
    METAHICT_RESOURCE_PEAK_FILE="${outdir}/.resources_peak.$$"
    METAHICT_LOG_FILE="${outdir}/module.log"
    METAHICT_RESOURCE_CGROUP_MEMORY_CURRENT=""
    if [[ -r /proc/self/cgroup ]]; then
        local cgroup_path
        cgroup_path="$(awk -F: '$2 == "" {print $3; exit}' /proc/self/cgroup 2>/dev/null || true)"
        if [[ "$cgroup_path" == */app.slice/*.service && -r "/sys/fs/cgroup${cgroup_path}/memory.current" ]]; then
            METAHICT_RESOURCE_CGROUP_MEMORY_CURRENT="/sys/fs/cgroup${cgroup_path}/memory.current"
        fi
    fi

    : > "$METAHICT_LOG_FILE"
    exec > >(tee -a "$METAHICT_LOG_FILE") 2>&1

    echo "[INFO] ===== METAHICT module started: ${module_name} ====="
    echo "[INFO] Start time: $(date '+%F %T')"
    echo "[INFO] Output dir: $outdir"
    echo "[INFO] Module log: $METAHICT_LOG_FILE"
    echo "[INFO] Resource report: $METAHICT_RESOURCE_FILE"

    _metahict_resource_monitor "$$" "$METAHICT_RESOURCE_PEAK_FILE" 2>/dev/null &
    METAHICT_RESOURCE_MONITOR_PID="$!"
    trap _metahict_resource_finish EXIT
}

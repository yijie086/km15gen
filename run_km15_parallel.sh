#!/bin/bash

set -u

# ===== 可配置 =====
OUTPUT_PREFIX="output_7p5"   # <<< 在这里改
N_BATCH=10
FILES_PER_BATCH=10000
MAX_CONCURRENT=40
# ==================

BASE_DIR="$(pwd)"
MAIN_PY="${BASE_DIR}/main.py"
LOG_BASE="${BASE_DIR}/logs_km15"

mkdir -p "$LOG_BASE"

run_one_file() {
    local batch=$1
    local i=$2

    local idx
    idx=$(printf "%05d" "$i")

    local outdir="${BASE_DIR}/${OUTPUT_PREFIX}_${batch}"
    local fname="km15_${idx}"
    local logfile="${LOG_BASE}/batch_${batch}_job_${idx}.log"

    cd "$outdir" || exit 1

    python "$MAIN_PY" \
        -trig 5000 \
        -fname "$fname" \
        -model km15 \
        -Ed 7.546 \
        -xBmin 0.1 \
        -xBmax 0.65 \
        -Q2min 1.0 \
        -Q2max 8.0 \
        -tmin 0 \
        -tmax 1.1 \
        -ymin 0.05 \
        -ymax 0.9 \
        -w2min 4.0 \
        > "$logfile" 2>&1

    if [ -f "${fname}.dat" ]; then
        mv "${fname}.dat" "${fname}.txt"
    fi
}

for ((batch=1; batch<=N_BATCH; batch++)); do
    OUTDIR="${BASE_DIR}/${OUTPUT_PREFIX}_${batch}"
    mkdir -p "$OUTDIR"

    echo "[INFO] Starting batch $batch -> $OUTDIR"

    for ((i=1; i<=FILES_PER_BATCH; i++)); do
        run_one_file "$batch" "$i" &

        while [ "$(jobs -rp | wc -l)" -ge "$MAX_CONCURRENT" ]; do
            sleep 1
        done
    done

    wait

    echo "[INFO] Finished batch $batch"
done

echo "[INFO] All batches finished."
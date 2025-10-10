#!/bin/bash
# ============================================================
#  N-Body Benchmark Script (Sequential vs Parallel)
# ============================================================

EXEC=./nbody_sim
OUTPUT=benchmark_results.csv

# Configurable test parameters
BODIES_LIST=(50 150 250 500 1000 2000)
STEPS_LIST=(10000 20000 50000 100000)
THREADS_LIST=(1 2 4 6 8 10 12)

# Create / reset output file
echo "bodies,threads,avg_force_time,avg_update_time,avg_total_step,total_time" > "$OUTPUT"

echo "Starting benchmark..."
echo "-------------------------------------"

for N in "${BODIES_LIST[@]}"; do
  for T in "${THREADS_LIST[@]}"; do
    for I in "${STEPS_LIST[@]}"; do
      echo "▶ Running: N=$N, Threads=$T"
      echo "-------------------------------------"

      # Run simulation and show live output
      RESULT=$($EXEC --bodies=$N --steps=$I --threads=$T --noview)

      # Print the full result to console
      echo "$RESULT"
      echo

      # Extract numeric values (remove trailing 's' and whitespace)
      AVG_FORCE=$(echo "$RESULT" | grep "Avg Force Time" | awk '{print $(NF-1)}')
      AVG_UPDATE=$(echo "$RESULT" | grep "Avg Update Time" | awk '{print $(NF-1)}')
      AVG_TOTAL=$(echo "$RESULT" | grep "Avg Total Step" | awk '{print $(NF-1)}')
      TOTAL=$(echo "$RESULT" | grep "Total Simulation Time" | awk '{print $(NF-1)}')

      # Append parsed results to CSV
      echo "$N,$T,$AVG_FORCE,$AVG_UPDATE,$AVG_TOTAL,$TOTAL" >> "$OUTPUT"
    done
  done
done

echo "-------------------------------------"
echo "Benchmark complete! Results saved to: $OUTPUT"

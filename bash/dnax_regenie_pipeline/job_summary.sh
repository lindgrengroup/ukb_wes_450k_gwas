#!/bin/bash

# User ID to query
USER_ID=${1:-"user-nbaya"}

# Number of recent jobs to analyze (adjust as needed)
NUM_JOBS=${2:-10000}

echo "Querying the last $NUM_JOBS jobs for user: $USER_ID..."

# 1. --no-subjobs: Prevents double counting costs (e.g. Workflow + Stage 1 + Stage 2...)
# 2. --json: Retrieves the metadata including totalPrice
dx find executions --user "$USER_ID" --no-subjobs --json -n "$NUM_JOBS" 2> /dev/null | \
jq -r '.[] | "\(.state) \(.totalPrice // 0)"' | \
awk '
BEGIN {
    print "---------------------------------------------"
    printf "%-15s %-10s %-15s\n", "STATE", "COUNT", "COST ($)"
    print "---------------------------------------------"
}
{
    count[$1]++
    cost[$1] += $2
    
    total_count++
    total_cost += $2
}
END {
    for (state in count) {
        printf "%-15s %-10d %-15.2f\n", state, count[state], cost[state]
    }
    print "---------------------------------------------"
    printf "%-15s %-10d %-15.2f\n", "TOTAL", total_count, total_cost
    print "---------------------------------------------"
}
'

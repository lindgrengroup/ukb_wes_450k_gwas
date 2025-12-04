#!/bin/bash

# --- 1. Check Input ---
if [[ -z "$1" ]]; then
  echo "Error: Missing bed file path." >&2
  echo "Usage: ./_get_INSTANCE_TYPE_for_bed_file.sh [BED_FILE_PATH]" >&2
  exit 1
fi

BED_FILE_PATH=$1

# --- 2. Extract Size Safely ---
# We capture the full line first to check if the file exists
line=$(dx ls -l "${BED_FILE_PATH}" 2> /dev/null)

if [[ -z "$line" ]]; then
  echo "Error: File not found or dx command failed: ${BED_FILE_PATH}" >&2
  exit 1
fi


INSTANCE_TYPE=$(echo "$line" | awk '{
    val = $4;
    unit = $5;
		
		# Validation: Ensure val is a valid number (integer or float)
    if (val !~ /^[0-9]+(\.[0-9]+)?$/) {
        # If not a number, exit awk without printing. 
        # The bash script checks for empty INSTANCE_TYPE below.
        exit 1
    }
    
    # Normalize everything to GB
    if (unit == "TB") {
        gb = val * 1024;
    } else if (unit == "GB") {
        gb = val;
    } else if (unit == "MB") {
        gb = val / 1024;
    } else if (unit == "KB") {
        gb = val / 1024 / 1024;
    } else {
        # Fallback: Assume Bytes if no unit matches
        gb = val / 1024 / 1024 / 1024;
    }

    # Select Instance Type based on GB size
    if (gb > 240) {
        print "mem1_ssd1_v2_x36";
    } else if (gb > 120) {
        print "mem1_ssd1_v2_x16";
    } else if (gb > 50) {
        print "mem1_ssd1_v2_x8";
    } else {
        print "mem1_ssd1_v2_x4";
    }
}')


echo $INSTANCE_TYPE

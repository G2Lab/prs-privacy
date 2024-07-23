#!/bin/bash

# Number of jobs to submit
total_chunks=501
chunk_size=5

# Maximum number of concurrent jobs
max_concurrent_jobs=10

# Job script name
job_script="solve.slurm"

# Function to count the number of running jobs
count_running_jobs() {
    squeue -u "$USER" --name "prs" -t R,PD | awk 'NR>1 {print $1}' | wc -l
}

# Function to check for job errors
check_job_errors() {
    completed_jobs=$(squeue -u "$USER" -t CD -o "%i" --noheader)
    for job_id in $completed_jobs; do
        status=$(sacct -j "$job_id" --format=State --noheader | grep -E "FAILED|CANCELLED")
        if [ -n "$status" ]; then
            echo "Job $job_id ended with error: $status"
        fi
    done
}

echo "Waiting to submit jobs..."
# Submit the jobs in batches
for ((i=9; i<total_chunks; i++)); do
    while [ $(count_running_jobs) -ge $max_concurrent_jobs ]; do
#        echo "Waiting for a job to finish..."
        sleep 60  # Wait for 60 seconds before checking again
        check_job_errors  # Check for job errors
    done

    echo "Submitting job: $job_script $i $chunk_size"
    sbatch $job_script $i $chunk_size
    sleep 5
done

# Final check for errors in case there are remaining jobs
while [ $(count_running_jobs) -gt 0 ]; do
    sleep 30
    check_job_errors
done

echo "All jobs submitted and monitored."

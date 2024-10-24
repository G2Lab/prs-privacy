import os
import subprocess
import re

# Directory containing log files
log_dir = "."  # Replace with your log directory

# Initialize variables
total_hours = 0
max_ram = 0
min_ram = float('inf')
total_ram = 0
num_jobs = 0


def convert_to_hours(elapsed):
    """Convert elapsed time (D-HH:MM:SS or HH:MM:SS) to hours."""
    if '-' in elapsed:
        days, time_part = elapsed.split('-')
        days = int(days)
        h, m, s = map(int, time_part.split(':'))
        return days * 24 + h + m / 60 + s / 3600
    else:
        h, m, s = map(int, elapsed.split(':'))
        return h + m / 60 + s / 3600


# Iterate over log files
for log_file in os.listdir(log_dir):
    if not log_file.startswith("prs_") or not log_file.endswith(".log"):
        continue

    # Read the log file
    with open(os.path.join(log_dir, log_file), 'r') as file:
        content = file.read()

    # Skip files containing the OOM killer message
    if "Some of your processes may have been killed by the cgroup out-of-memory handler" in content:
        continue

    # Extract job ID from file name
    job_id = re.findall(r'\d+', log_file)[0]

    # Retrieve elapsed time and max RAM usage with sacct
    result = subprocess.run(
        ["sacct", "-j", job_id, "--format=Elapsed,MaxRSS", "--parsable2", "--noheader", "--units=G"],
        capture_output=True, text=True
    )

    if result.returncode != 0 or not result.stdout.strip():
        continue

    # Parse the sacct output
    lines = result.stdout.strip().split('\n')
    elapsed, max_rss = lines[1].split('|')

    # Convert elapsed time to hours
    hours = convert_to_hours(elapsed)

    # Convert max RAM usage to GB
    if max_rss.endswith("G"):
        ram_gb = float(max_rss[:-1])
        total_ram += ram_gb
        num_jobs += 1
    else:
        ram_gb = 0

    # Add to total hours
    total_hours += hours

    # Determine max and min RAM usage
    max_ram = max(max_ram, ram_gb)
    if ram_gb > 0:
        min_ram = min(min_ram, ram_gb)

# Handle the case where no valid files are found
if min_ram == float('inf'):
    min_ram = 0

# Print results
print(f"Total runtime: {total_hours:.2f} hours")
print(f"Maximum RAM used: {max_ram:.2f} GB")
print(f"Minimum RAM used: {min_ram:.2f} GB")
print(f"Average RAM used: {total_ram / num_jobs:.2f} GB")

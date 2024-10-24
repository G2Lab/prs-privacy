import os

# Define the directory where the files are located
directory = 'results/sequential/'

# Define the expected range of numbers
expected_numbers = set(range(1, 500))

# Find the numbers of the files that actually exist
existing_numbers = set()

for filename in os.listdir(directory):
    if filename.startswith("guesses") and filename.endswith(".json"):
        # Extract the number part of the filename
        number = int(filename[7:-5])  # Extract the digits from the filename (from 'guesses' to '.json')
        existing_numbers.add(number)

# Find the missing numbers
missing_numbers = expected_numbers - existing_numbers

# Print the missing numbers
if missing_numbers:
    print("Missing chunks:")
    for num in sorted(missing_numbers):
        print(f"guesses{num:03}.json")
else:
    print("No files are missing.")

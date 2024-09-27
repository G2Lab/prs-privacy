import json

# Load the JSON data from files
def load_json(filename):
    with open(filename, 'r') as file:
        return json.load(file)

# Save the JSON data to a file
def save_json(data, filename):
    with open(filename, 'w') as file:
        json.dump(data, file, indent=4)

# Find missing results and update the second file
def merge_results(file1, file2, output_file):
    data1 = load_json(file1)
    data2 = load_json(file2)

    pgs_ids_1 = {result['PgsID'] for result in data1}

    # Find and copy missing results from the second file
    for result in data2:
        if result['PgsID'] not in pgs_ids_1:
            # Modify the values as required
            result['PredictedMedianAnonymitySet'] = -1
            result['PredictedPercentageUnique'] = -1
            data1.append(result)

    save_json(data1, output_file)

merge_results('results/uniqueness/scores_1000Genomes.json', 'results/uniqueness/scores_1000Genomes_full.json',
              'results/uniqueness/scores_1000Genomes_merged.json')
merge_results('results/uniqueness/scores_UKBB.json', 'results/uniqueness/scores_UKBiobank_full.json',
              'results/uniqueness/scores_UKBB_merged.json')
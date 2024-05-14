from collections import defaultdict
import csv
import fnmatch
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import statistics


def pairwise(filepath):
    correlations = []
    with open(filepath, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            row_values = [float(value) for value in row]
            correlations.append(row_values)

    xs = np.arange(1, len(correlations)+1, 1)
    ys = []
    for i in range(len(correlations)):
        ys.append([correlations[i][j] for j in range(i+1, len(correlations))])
        plt.scatter([xs[i] for _ in range(len(ys[i]))], ys[i], s=1)

    plt.xlabel("Distance between a pair of SNPs")
    plt.ylabel("Standard deviation of pairwise conditional probabilities")
    plt.tight_layout()
    # plt.show()
    plt.savefig('pairwise.png', dpi=300)


def score_distribution(filepath):
    scores = []
    with open(filepath, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            scores.append(float(row[1]))

    print("Median: ", np.median(scores))
    print("Percentiles: ", np.percentile(scores, range(5, 100, 5)))
    plt.hist(scores, bins=100)
    plt.xlabel("Score")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.show()
    # plt.savefig('score_distribution.png', dpi=300)


def solution_ranking(pgs_id):
    accuracies, positions, scores, likelihoods = [], [], [], []
    major_accuracies, major_scores, major_likelihoods = [], [], []
    directory = "results/"
    filepaths = [os.path.join(directory, filename) for filename in os.listdir(directory) if
                 fnmatch.fnmatch(filename, f"*{pgs_id}*")]
    print(filepaths)
    for filepath in filepaths:
        with open(filepath, 'r') as file:
            reader = csv.reader(file)
            try:
                next(reader)
            except StopIteration:
                # Empty file
                continue
            for row in reader:
                if len(row) < 11:
                    continue
                score, major_score, median_score = float(row[1]), float(row[2]), float(row[3])
                accuracy, major_accuracy = float(row[4]), float(row[5])
                position, total_solutions = float(row[6]), float(row[7])
                first_likelihood, likelihood, major_likelihood = float(row[8]), float(row[9]), float(row[10])
                accuracies.append(accuracy)
                major_accuracies.append(major_accuracy)
                positions.append(1-position/total_solutions)
                likelihoods.append(likelihood)
                major_likelihoods.append(major_likelihood*score/major_score)

    # plt.hist(accuracies, bins=50, histtype='step', label="Top-choice accuracy")
    # plt.hist(major_accuracies, bins=50, histtype='step', label="All-major accuracy")
    # plt.hist(positions, bins=50, histtype='step', label="True solution proximity to the top")
    plt.hist(likelihoods, bins=50, histtype='step', label="Likelihood ratio")
    plt.hist(major_likelihoods, bins=50, histtype='step', label="Likelihood ratio")
    plt.title(pgs_id)
    plt.xlabel("Ratio")
    plt.ylabel("Count")
    plt.legend()
    plt.tight_layout()
    plt.show()


def likelihood_distribution(pgs_id):
    directory = "results/scoreLikelihood/"
    filepaths = [os.path.join(directory, filename) for filename in os.listdir(directory) if
                 fnmatch.fnmatch(filename, f"*{pgs_id}*")]
    print(filepaths)
    labels = []
    scores = []
    likelihoods = []
    for filepath in filepaths:
        with open(filepath, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                label, score, *likelihood = row
                labels.append(label)
                scores.append(score)
                likelihoods.append(list(map(lambda x: float(x), likelihood)))

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    for i in range(0, int(len(labels)/2)):
        plt.hist(likelihoods[2*i], bins=8, histtype='step', linewidth=2, color=colors[i], label=f"{labels[i]}({scores[2*i]})")
        plt.hist(likelihoods[2*i+1], bins=8, histtype='step', linewidth=2, color=colors[i], label=f"{labels[i]}({scores[2*i+1]})")
        plt.plot(likelihoods[2*i][0], 1600, 'x', markersize=5, color=colors[i])
        plt.plot(likelihoods[2*i+1][0], 1600, 'x', markersize=5, color=colors[i])
    plt.title(pgs_id)
    plt.xlabel("Likelihood")
    plt.ylabel("Count")
    plt.legend()
    plt.tight_layout()
    plt.show()


def score_to_likelihood(pgs_id, num_snps):
    directory = "results/scoreLikelihood/"
    filepaths = [os.path.join(directory, filename) for filename in os.listdir(directory) if
                 fnmatch.fnmatch(filename, f"*{pgs_id}*")]
    print(filepaths)
    labels = []
    scores = []
    likelihoods = []
    distribution = []
    for filepath in filepaths:
        with open(filepath, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                if len(row) < 3:
                    continue
                label, score, likelihood, *distro = row
                labels.append(label)
                scores.append(float(score))
                likelihoods.append(float(likelihood))
                if len(distro) > 0:
                    distribution.append(list(map(lambda x: float(x) if x else 0.0, distro)))
                else:
                    distribution.append([])

    print("Number of samples: ", len(scores))
    transparency = 1
    marker_size = 0.5
    plt.scatter(scores, likelihoods, color='black', s=1, label="Correct solutions")
    plt.scatter([scores[0] for x in range(0, len(distribution[0]))], distribution[0], alpha=transparency, color='moccasin',
                s=marker_size, label="Incorrect solutions")
    for i in range(1, len(scores)):
        plt.scatter([scores[i] for x in range(0, len(distribution[i]))], distribution[i], alpha=transparency, color='moccasin',
                    s=marker_size)
    plt.title(pgs_id + f" ({num_snps} variants)")
    plt.xlabel("Score")
    plt.ylabel("Solution likelihood")
    plt.legend()
    plt.tight_layout()
    plt.show()


def score_to_likelihood_position(pgs_id, num_snps):
    directory = "results/scoreLikelihood/"
    filepaths = [os.path.join(directory, filename) for filename in os.listdir(directory) if
                 fnmatch.fnmatch(filename, f"*{pgs_id}*")]
    labels = []
    scores = []
    likelihoods = []
    distribution = []
    for filepath in filepaths:
        with open(filepath, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                if len(row) < 3:
                    continue
                label, score, likelihood, *distro = row
                labels.append(label)
                scores.append(float(score))
                likelihoods.append(float(likelihood))
                if len(distro) > 0:
                    distribution.append(list(map(lambda x: float(x) if x else 0.0, distro)))
                else:
                    distribution.append([])

    print("Number of samples: ", len(scores))
    proximity = []
    for i in range(0, len(scores)):
        if len(distribution[i]) > 0:
            position = bisect(sorted(distribution[i], reverse=True), likelihoods[i])
        else:
            position = 0
        proximity.append(float(position)*100/float(len(distribution[i])+1))
        # print(f"{position}/{len(distribution[i])+1}", end=" ")

    plt.scatter(scores, proximity, color='sandybrown', s=1, label="Correct solutions")
    plt.title(pgs_id + f" ({num_snps} variants)")
    plt.xlabel("Score")
    plt.ylabel("True solution in the top X% by likelihood")
    plt.tight_layout()
    # plt.show()
    plt.savefig(pgs_id + '-pos.pdf')


def bisect(array, value):
    for i in range(0, len(array)):
        if array[i] <= value:
            return i
    return len(array)


def accuracy_likelihood(pgs_id):
    directory = "results/accuracyLikelihood/"
    filepaths = [os.path.join(directory, filename) for filename in os.listdir(directory) if
                 fnmatch.fnmatch(filename, f"*{pgs_id}-*")]
    scores, accuracies, likelihoods = [], [], []
    chosen_acc, chosen_like = [], []
    true_like = []
    lowest_like, lowest_like_acc = [], []
    for filepath in filepaths:
        print(filepath)
        with open(filepath, 'r') as file:
            for row in file:
                # print(row)
                row = row.strip()
                data = json.loads(row)
                scores.append(float(data["Score"]))
                chosen_acc.append(float(data["Accuracies"][0]))
                chosen_like.append(float(data["Likelihoods"][0]))
                acc = list(map(lambda x: float(x), data["Accuracies"]))
                true_like.append(float(data["Likelihoods"][acc.index(1.0)]) if 1.0 in acc else 0.0)
                lowest_index = acc.index(min(acc))
                lowest_like.append(float(data["Likelihoods"][lowest_index]))
                lowest_like_acc.append(float(data["Accuracies"][lowest_index]))

    # plt.scatter(chosen_like, chosen_acc, color='black', s=1, label="Chosen solutions")
    # plt.scatter(true_like, [1.0 for x in range(0, len(true_like))], color='turquoise', alpha=0.2, s=0.5, label="True solutions")
    # plt.scatter(lowest_like, lowest_like_acc, color='salmon', alpha=0.2, s=0.5, label="Lowest-likelihood solutions")
    # plt.scatter(true_like, [1.0 for x in range(0, len(true_like))], color='turquoise', alpha=0.5, s=1, label="True solutions")
    # plt.scatter(lowest_like, lowest_like_acc, color='salmon', alpha=0.5, s=1, label="Lowest-likelihood solutions")
    plt.scatter(true_like, chosen_acc, color='purple', alpha=0.5, s=1)
    # for i in range(0, len(scores)):
    #     plt.scatter(likelihoods[i], accuracies[i], alpha=0.05, color='moccasin', s=0.1)
    plt.title(pgs_id)
    plt.xlabel("Likelihood")
    plt.ylabel("Accuracy")
    lgnd = plt.legend()
    for i in range(0, len(lgnd.legend_handles)):
        lgnd.legend_handles[i].set_sizes([10])
    plt.tight_layout()
    plt.show()
    # plt.savefig(pgs_id + '-eaf.pdf')


def true_position_cdf(pgs_ids):
    colors = ['dodgerblue', 'limegreen', 'tomato', 'blueviolet']
    directory = "results/accuracyLikelihood/"
    positions = []
    for pgs in pgs_ids:
        filepath = os.path.join(directory, pgs+".json")
        with open(filepath, 'r') as file:
            for row in file:
                row = row.strip()
                data = json.loads(row)
                accuracy = list(map(lambda x: float(x), data["Accuracies"]))
                positions.append(float(accuracy.index(1.0))/float(len(accuracy)))

        spos = [pos*100 for pos in np.sort(positions)]
        cdf = np.arange(1, len(spos) + 1) / len(spos)
        plt.plot(spos, cdf, label=pgs, linewidth=2, color=colors[pgs_ids.index(pgs)])

    # plt.title("CDF for the likelihood position of the true solution")
    plt.xlabel("Selected solution in the top X% by likelihood", fontsize=16)
    plt.ylabel("Fraction", fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend()
    plt.tight_layout()
    plt.grid(True)
    # plt.show()
    plt.savefig('cdf.pdf')


def accuracy_cdf(pgs_ids):
    colors = ['firebrick', 'sandybrown', 'forestgreen', 'mediumorchid']
    directory = "results/accuracyLikelihood/"
    accuracy = []
    for pgs in pgs_ids:
        filepath = os.path.join(directory, pgs+".json")
        with open(filepath, 'r') as file:
            for row in file:
                row = row.strip()
                data = json.loads(row)
                accuracy.append(float(data["Accuracies"][0]))

        spos = [pos for pos in np.sort(accuracy)]
        cdf = np.arange(1, len(spos) + 1) / len(spos)
        plt.plot(spos, cdf, label=pgs, linewidth=2, color=colors[pgs_ids.index(pgs)])

    # plt.title("CDF for the likelihood position of the true solution")
    plt.xlabel("Top-solution accuracy", fontsize=16)
    plt.ylabel("Fraction", fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.gca().invert_xaxis()
    plt.legend()
    plt.tight_layout()
    plt.grid(True)
    # plt.show()
    plt.savefig('accuracy.pdf')


def score_deviation():
    filepath = "results/distro.json"
    rows = []
    with open(filepath, 'r') as file:
        for row in file:
            row = row.strip()
            if row.startswith("[["):
                rows.append(list(map(lambda x: [float(x[0]), float(x[1])], json.loads(row))))
            else:
                rows.append(list(map(lambda x: float(x), json.loads(row))))
    maxTotal = rows[0]
    examples = rows[1]
    left = rows[2]
    right = rows[3]
    # diff = [left[i]-right[i] for i in range(0, len(left))]
    full = [left[i]+right[i] for i in range(0, len(left))]
    print(len(left), len(right))
    mf = statistics.median(full)
    leftLim, rightLim = mf - maxTotal[0], mf - maxTotal[1]
    plt.axvline(leftLim, color='purple', linewidth=1, label='Naive Min')
    plt.axvline(rightLim, color='purple', linewidth=1, label='Naive Max')

    # plt.hist(diff, bins=50, color='blue', alpha=0.7, label='Left-right diff')
    plt.hist(full, bins=50, color='orange', alpha=0.7)
    # plt.hist(left, bins=50, color='blue', alpha=0.7, label='Scores')
    # plt.hist(right, bins=50, color='orange', alpha=0.7, label='Scores')
    labels = ['Sampled Max', 'Sampled Min']
    colors = ['red', 'green', 'yellow', 'purple']
    for i in range(0, len(examples)):
        plt.axvline(examples[i], color=colors[i], linestyle='dashed', linewidth=1, label=labels[i])
    # for example in examples:
    #     plt.axvline(example[0], color=colors[0], linestyle='dashed', linewidth=1, label=labels[0])
    #     plt.axvline(example[1], color=colors[1], linestyle='dashed', linewidth=1, label=labels[1])
        # plt.axvline(example[0]+example[1], color=colors[1], linestyle='dashed', linewidth=1, label=labels[1])
    # for v in major:
    #     plt.axvline(v, color=colors[major.index(v)], linestyle='dashed', linewidth=1, label=labels[major.index(v)])
    # for v in minor:
    #     plt.axvline(v, color=colors[minor.index(v)+2], linestyle='dashed', linewidth=1, label=labels[minor.index(v)])
    plt.xlabel("Score")
    plt.ylabel("Frequency")
    plt.legend()
    plt.show()


def weight_to_af():
    # Load the JSON data
    # with open('results/weights/PGS000648.json', 'r') as f:
    with open('results/weights/PGS002302.json', 'r') as f:
        data = json.load(f)

    # Extract the Weights and AF data
    weights = data['Weights']
    af_data = data['AF']

    # Sort the weights in ascending order
    sorted_indices = sorted(range(len(weights)), key=lambda i: weights[i])
    sorted_weights = [weights[i] for i in sorted_indices]

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the sorted weights on the X-axis
    ax.plot(sorted_weights, [0] * len(sorted_weights), 'ko', markersize=2)
    ax.set_xlabel('Sorted Weights')

    # Plot the corresponding AF values for each key on the Y-axis
    colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    for i, (key, af_values) in enumerate(af_data.items()):
        sorted_af = [af_values[j] for j in sorted_indices]
        ax.plot(sorted_weights, sorted_af, f'{colors[i % len(colors)]}-', label=key)

    # Set the Y-axis label and legend
    ax.set_ylabel('AF Values')
    ax.legend()

    # Show the plot
    plt.show()


def likelihood_spectrum_select(pgs_id):
    num_variants = {
        "PGS000073": 10,
        "PGS000037": 15,
        "PGS000639": 20,
        "PGS002302": 28,
        "PGS000154": 30,
        "PGS000851": 37,
        "PGS003760": 49,
    }
    pops = ["AFR", "AMR", "EAS", "EUR", "SAS"]
    directory = "results/sorting/"
    filepath = os.path.join(directory, pgs_id+".json")
    with open(filepath, 'r') as f:
        data = json.load(f)
    individuals, ancestry, scores, true_lkl, true_spec, ref_acc = [], [], [], [], [], []
    lkl_acc, spec_acc, lkl_spec_acc = [], [], []
    for result in data:
        individuals.append(result["Individual"])
        ancestry.append(result["Ancestry"])
        scores.append(float(result["Score"]))
        true_lkl.append(float(result["TrueLikelihood"]))
        true_spec.append(float(result["TrueSpectrum"]))
        ref_acc.append(float(result["ReferenceAccuracy"]))
        lkl_acc.append(list(map(lambda x: float(x), result["LikelihoodAccuracies"])))
        spec_acc.append(list(map(lambda x: float(x), result["SpectrumAccuracies"])))
        lkl_spec_acc.append(list(map(lambda x: float(x), result["LikelihoodSpectrumAccuracies"])))

    # print(f"Median and mean likelihood accuracy: {np.median([l[0] for l in lkl_acc])}, {np.mean([l[0] for l in lkl_acc])}")
    # print(f"Median and mean spectrum accuracy: {np.median([l[0] for l in spec_acc])}, {np.mean([l[0] for l in spec_acc])}")
    # print(f"Median and mean likelihood+spectrum accuracy: {np.median([l[0] for l in lkl_spec_acc])}, {np.mean([l[0] for l in lkl_spec_acc])}")

    colors = ['firebrick', 'sandybrown', 'forestgreen', 'mediumorchid', 'dodgerblue', 'limegreen', 'tomato']
    markers = ['.', '*', 'v', "^", "+"]

    # Create a figure and axis
    fig1, ax1 = plt.subplots(figsize=(12, 4))
    for i, pop in enumerate(pops):
        score, acc = [], []
        for j, ind in enumerate(individuals):
            if pop == ancestry[j]:
                score.append(scores[j])
                acc.append(lkl_acc[j][0])
        ax1.scatter(score, acc, color=colors[i], s=10, label=pop)

    fig2, ax2 = plt.subplots(figsize=(12, 4))
    tlkl, acc = [], []
    for j, ind in enumerate(individuals):
        tlkl.append(true_lkl[j])
        acc.append(lkl_acc[j][0])

    ax2.scatter(tlkl, acc, color=colors[0], s=10)

    # Set axis labels and title
    ax1.set_xlabel("Score")
    ax1.set_ylabel("Accuracy")
    ax1.set_title(f"{pgs_id}({num_variants[pgs_id]} variants): The effect of ancestry and score")

    ax2.set_xlabel("True Likelihood (less is more common)")
    ax2.set_ylabel("Accuracy")
    ax2.set_title(f"{pgs_id}({num_variants[pgs_id]} variants): The effect of the true likelihood")

    # Adjust the legend
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    fig1.savefig(f'{pgs_id}_ancestry.png', dpi=300, bbox_inches='tight')
    fig2.savefig(f'{pgs_id}_likelihood.png', dpi=300, bbox_inches='tight')
    # Show the plot
    # plt.show()


def likelihood_spectrum_total():
    num_variants = {
        "PGS000073": 10,
        "PGS000037": 15,
        "PGS000639": 20,
        "PGS002302": 28,
        "PGS000154": 30,
        "PGS000851": 37,
        "PGS003760": 49,
    }
    directory = "results/sorting/"
    pgs_ids = ["PGS000073", "PGS000037", "PGS000639", "PGS002302", "PGS000154", "PGS000851", "PGS003760"]
    loci, reference_median, reference_mean = [], [], []
    likelihoods_median, spectra_median, likelihoods_spectra_median = [], [], []
    likelihoods_mean, spectra_mean, likelihoods_spectra_mean = [], [], []
    for pgs_id in pgs_ids:
        filepath = os.path.join(directory, pgs_id+".json")
        with open(filepath, 'r') as f:
            data = json.load(f)
        individuals, ancestry, scores, true_lkl, true_spec, ref_acc = [], [], [], [], [], []
        lkl_acc, spec_acc, lkl_spec_acc = [], [], []
        for result in data:
            individuals.append(result["Individual"])
            ancestry.append(result["Ancestry"])
            scores.append(float(result["Score"]))
            true_lkl.append(float(result["TrueLikelihood"]))
            true_spec.append(float(result["TrueSpectrum"]))
            ref_acc.append(float(result["ReferenceAccuracy"]))
            lkl_acc.append(list(map(lambda x: float(x), result["LikelihoodAccuracies"])))
            spec_acc.append(list(map(lambda x: float(x), result["SpectrumAccuracies"])))
            lkl_spec_acc.append(list(map(lambda x: float(x), result["LikelihoodSpectrumAccuracies"])))
        #
        loci.append(num_variants[pgs_id])
        reference_median.append(np.median(ref_acc))
        likelihoods_median.append(np.median([l[0] for l in lkl_acc if len(l) > 0]))
        spectra_median.append(np.median([l[0] for l in spec_acc if len(l) > 0]))
        likelihoods_spectra_median.append(np.median([l[0] for l in lkl_spec_acc if len(l) > 0]))
        reference_mean.append(np.mean(ref_acc))
        likelihoods_mean.append(np.mean([l[0] for l in lkl_acc if len(l) > 0]))
        spectra_mean.append(np.mean([l[0] for l in spec_acc if len(l) > 0]))
        likelihoods_spectra_mean.append(np.mean([l[0] for l in lkl_spec_acc if len(l) > 0]))
        if pgs_id == "PGS003760":
            print(f"PGS003760: {len(individuals)} samples, "
                  f"{len([l[0] for l in lkl_acc if len(l) > 0])} likelihoods, "
                  f"{len([l[0] for l in spec_acc if len(l) > 0])} spectra, "
                  f"{len([l[0] for l in lkl_spec_acc if len(l) > 0])} likelihood+spectra")

    colors = ['firebrick', 'sandybrown', 'forestgreen', 'mediumorchid', 'dodgerblue', 'limegreen', 'tomato']
    markers = ['.', '*', 'v', "^", "+"]

    # Create a figure and axis
    fig1, ax1 = plt.subplots(figsize=(6, 4))
    ax1.plot(loci, likelihoods_median, color=colors[0], marker=".", markersize=8, label="Likelihood")
    ax1.plot(loci, spectra_median, color=colors[1], marker=".", markersize=8,label="Spectrum")
    ax1.plot(loci, likelihoods_spectra_median, color=colors[2], marker=".", markersize=8, label="Likelihood+Spectrum")
    ax1.plot(loci, reference_median, color='k', linewidth=0.5, linestyle='--', marker=".", markersize=5, label="Reference")

    for i, tup in enumerate(zip(likelihoods_median, spectra_median, likelihoods_spectra_median)):
        print(f"{likelihoods_median[i]:.3f}/{likelihoods_mean[i]:.3f}\t{spectra_median[i]:.3f}/{spectra_mean[i]:.3f}\t"
              f"{likelihoods_spectra_median[i]:.3f}/{likelihoods_spectra_mean[i]:.3f}\t{reference_median[i]:.3f}/{reference_mean[i]:.3f}")

    # Set axis labels and title
    ax1.set_xlabel("Number of variants")
    ax1.set_ylabel("Accuracy")
    ax1.set_title("Median Accuracy")
    # ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.legend()

    # Plot for mean
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    ax2.plot(loci, likelihoods_mean, color=colors[0], marker=".", markersize=8, label="Likelihood")
    ax2.plot(loci, spectra_mean, color=colors[1], marker=".", markersize=8, label="Spectrum")
    ax2.plot(loci, likelihoods_spectra_mean, color=colors[2], marker=".", markersize=8, label="Likelihood+Spectrum")
    ax2.plot(loci, reference_mean, color='k', linewidth=0.5, linestyle='--', marker=".", markersize=5, label="Reference")

    # Set axis labels and title for mean plot
    ax2.set_xlabel("Number of variants")
    ax2.set_ylabel("Accuracy")
    ax2.set_title("Mean Accuracy")

    # Adjust the legend for mean plot
    # ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.legend()

    plt.tight_layout()
    fig1.savefig('median.png', dpi=300, bbox_inches='tight')
    fig2.savefig('mean.png', dpi=300, bbox_inches='tight')
    # Show the plot
    # plt.show()


def loci_coverage():
    with open('loci_pgp_coverage.json', 'r') as file:
        data = json.load(file)
    locus_counts = {}

    for location, studies in data.items():
        locus_counts[location] = len(studies)

    sorted_loci = sorted(locus_counts.keys(), key=lambda x: int(x.split(':')[0]))
    loci = []
    counts = []
    for locus in sorted_loci:
        loci.append(locus)
        counts.append(locus_counts[locus])

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(15, 5))

    # Plot the bar graph
    bar_width = 1
    ax.bar(range(len(counts)), counts, bar_width, align='center')

    current_chr = 1
    positions = []
    for locus in loci:
        chr_num = int(locus.split(':')[0])
        if chr_num != current_chr:
            ax.axvline(x=loci.index(locus), color='red', linestyle="dotted", linewidth=0.5)
            current_chr = chr_num
            positions.append(loci.index(locus))

    ax.set_xticks(positions)
    ax.set_xticklabels([l for l in range(2, 23)])

    ax.set_ylabel('Frequency')
    ax.set_xlabel('Chromosome')
    ax.set_title('Distribution of Loci per Chromosome')

    plt.tight_layout()
    fig.savefig('coverage.png', dpi=300, bbox_inches='tight')
    # plt.show()


def accuracy(pgs_ids):
    num_variants = {
            "PGS003181": 15,
            "PGS000778": 20,
            "PGS004249": 25,
            "PGS001868": 30,
            "PGS002270": 33,
            "PGS001835": 38,
        }
    world = "All"
    pops = ["AFR", "AMR", "EAS", "EUR", "SAS"]
    directory = "results/accuracy/"
    data = {}
    for pgs_id in pgs_ids+[world]:
        data[pgs_id] = {
                               'Ancestry': [],
                               'Accuracy': [],
                               'Score': []
                           }
    for pgs_id in pgs_ids:
        filepaths = [os.path.join(directory, filename) for filename in os.listdir(directory) if
                         fnmatch.fnmatch(filename, f"*{pgs_id}-*")]
        for filepath in filepaths:
            with open(filepath, 'r') as f:
                content = json.load(f)
            for result in content:
                if len(result["LikelihoodAccuracies"]) > 0:
                    data[pgs_id]['Ancestry'].append(result["Ancestry"])
                    data[pgs_id]['Accuracy'].append(float(result["LikelihoodAccuracies"][0]))
                    data[pgs_id]['Score'].append(float(result["Score"]))
                    data[world]['Ancestry'].append(result["Ancestry"])
                    data[world]['Accuracy'].append(float(result["LikelihoodAccuracies"][0]))
                    data[world]['Score'].append(float(result["Score"]))

    for pgs_id in pgs_ids + [world]:
        print(f"PGS: {pgs_id}, Mean accuracy: {np.mean(data[pgs_id]['Accuracy'])}, "
              f"Median accuracy: {np.median(data[pgs_id]['Accuracy'])}")
        df = pd.DataFrame(data[pgs_id])
        fig, ax = plt.subplots(figsize=(4, 3))
        sns.violinplot(x='Ancestry', y='Accuracy', data=df, ax=ax, palette='Pastel1', hue='Ancestry')
        if pgs_id != world:
            ax.set_title(f'Accuracy by Ancestry {pgs_id} ({num_variants[pgs_id]} variants)')
        else:
            ax.set_title(f'Accuracy by Ancestry {pgs_id}')
        fig.savefig(f'{pgs_id}-ancestry.png', dpi=300, bbox_inches='tight')
#         plt.tight_layout()
#         plt.show()

        decile_labels = ['10%', '20%', '30%', '40%', '50%',
                         '60%', '70%', '80%', '90%', '100%']

        df['Decile'] = pd.qcut(df['Score'], q=10, labels=decile_labels)
        fig, ax = plt.subplots(figsize=(4, 3))
#         sns.boxplot(x='Decile', y='Accuracy', data=df, ax=ax, width=0.5, fliersize=2, medianprops=dict(color='white'))
        sns.boxplot(x='Decile', y='Accuracy', data=df, ax=ax, width=0.5, fliersize=2)
        if pgs_id != world:
            ax.set_title(f'Accuracy by Score {pgs_id} ({num_variants[pgs_id]} variants)')
        else:
            ax.set_title(f'Accuracy by Score {pgs_id}')
        plt.xticks(rotation=45)
        fig.savefig(f'{pgs_id}-score.png', dpi=300, bbox_inches='tight')
#         plt.tight_layout()
#         plt.show()


def sequential():
    directory = "results/sequential/"
    files = os.listdir(directory)
    data = []
    for file_name in files:
        if "with_repair" not in file_name:
#         if "without_repair" not in file_name:
            continue
        filepath = os.path.join(directory, file_name)
        with open(filepath, 'r') as f:
            content = json.load(f)
        for result in content:
            individual = result['Individual']
            accuracies = result['Accuracies']
            for label, value in accuracies.items():
                data.append({'Individual': individual, 'Label': int(label), 'Value': value})

    df = pd.DataFrame(data)
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.boxplot(x='Label', y='Value', data=df, color='darkorange')
    plt.title('Accuracy distribution by SNP count')
    plt.xlabel('#SNPs targeted')
    plt.ylabel('Accuracy')
    plt.tight_layout()
    plt.ylim(0.87, 1.01)
    plt.show()
#     fig.savefig('sequential-unrepaired.png', dpi=300, bbox_inches='tight')
    fig.savefig('sequential-repaired.png', dpi=300, bbox_inches='tight')


def king_test():
	directory = "results/kinship/"
	filepath = os.path.join(directory, "truth.json")
	data = []
	with open(filepath, 'r') as f:
		content = json.load(f)
	for num_snps, results in content.items():
		print(np.median([int(x['Position']) for x in results]))
		for result in results:
			data.append({'SNPs': int(num_snps), 'TruePhi': float(result['TruePhi']), 'TopPhi': float(result['HighPhi']),
			 'Position': int(result['Position'])})

	df = pd.DataFrame(data)
	fig, ax = plt.subplots(figsize=(4, 3))
# 	sns.violinplot(x='SNPs', y='Position', ax=ax, data=df)
# 	plt.gca().invert_yaxis()
	sns.boxplot(x='SNPs', y='TruePhi', data=df)
	plt.axhline(y=0.25, linestyle='--', color='r')
	plt.axhline(y=0.125, linestyle='--', color='r')
	plt.xlabel('#SNPs')
	plt.title('Phi of a relative in the KING test')
	plt.ylabel('Phi')
# 	plt.title('Position of a relative by KING test')
# 	plt.ylabel('Position')
	plt.tight_layout()
# 	fig.savefig('king-position.png', dpi=300, bbox_inches='tight')
	fig.savefig('king-phi.png', dpi=300, bbox_inches='tight')
# 	plt.show()


def kinship_experiment():
	directory = "results/kinship/"
	num_snps = [2000, 2500]
# 	num_snps = [2000]
	data = []
	for snps in num_snps:
		filepath = os.path.join(directory, f"{snps}.json")
		with open(filepath, 'r') as f:
			content = json.load(f)
		print(f"{snps}: {np.median([int(x['Position']) for x in content])}")
		for result in content:
			data.append({'SNPs': snps, 'TruePhi': float(result['TruePhi']), 'TopPhi': float(result['HighPhi']),
			 'Position': int(result['Position']), 'Accuracy': float(result['Accuracy'])})

	df = pd.DataFrame(data)
	fig, ax = plt.subplots(figsize=(4, 3))
	sns.boxplot(x='SNPs', y='Accuracy', data=df)
# 	plt.gca().invert_yaxis()
# 	sns.boxplot(x='SNPs', y='Position', data=df)
# 	sns.boxplot(x='SNPs', y='TruePhi', data=df)
# 	plt.axhline(y=0.25, linestyle='--', color='r')
# 	plt.axhline(y=0.125, linestyle='--', color='r')
# 	plt.title('Phi of a relative by KING test')
	plt.title('Position of a relative by KING test')
	plt.xlabel('#SNPs')
# 	plt.ylabel('Phi')
	plt.ylabel('Position')
	plt.tight_layout()
# 	fig.savefig('kinship-position.png', dpi=300, bbox_inches='tight')
# 	fig.savefig('kinship-phi.png', dpi=300, bbox_inches='tight')
	plt.show()


def score_uniqueness():
	directory = "results/uniqueness/"
	filepath = os.path.join(directory, "scores.json")
	data = []
	with open(filepath, 'r') as f:
		content = json.load(f)
	for result in content:
		data.append({'Number of Variants': int(result['NumVariants']),
		'Mean Anonymity Set Size': np.mean([float(x) for x in result['AnonymitySets']])})

	df = pd.DataFrame(data)
	fig, ax = plt.subplots(figsize=(4, 3))
	sns.lineplot(x='Number of Variants', y='Mean Anonymity Set Size', data=df)
	plt.title('Score uniqueness in 1000 genomes')
	plt.xlabel('Number of Variants in PRS')
	plt.tight_layout()
	fig.savefig('uniqueness.png', dpi=300, bbox_inches='tight')
# 	plt.show()


def random_hist():
    data = np.random.beta(0.5, 0.5, 100)

    # Create a histogram using seaborn
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.histplot(data, bins=20, kde=False, color='darkorange')
#     plt.hist(data, bins=30, color='darkorange', linewidth=2)
    plt.xlabel('Allele Frequency')
    plt.ylabel('Count')

    # Show the plot
    fig.savefig('hist.png', dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    # pairwise("data/prior/PGS000040.pairwise")
    # score_distribution("PGS000040.scores")
    # solution_ranking("PGS000639")
    # solution_ranking("PGS002302")
    # solution_ranking("PGS000073")
    # solution_ranking("PGS000037")
    # score_to_likelihood("PGS000037", 15)
    # score_to_likelihood("PGS000073", 10)
    # score_to_likelihood("PGS000639", 20)
    # score_to_likelihood("PGS002302", 28)
    # score_to_likelihood("PGS001827", 33)
    # score_to_likelihood_position("PGS002302", 28)
    # score_to_likelihood_position("PGS000073", 10)
    # score_to_likelihood_position("PGS000037", 15)
    # score_to_likelihood_position("PGS000639", 20)
    # score_to_likelihood_position("PGS001827", 33)
    # accuracy_likelihood("PGS000040")
    # accuracy_likelihood("PGS000073")
    # accuracy_likelihood("PGS000037")
    # accuracy_likelihood("PGS002302")
    # true_position_cdf(["PGS000037", "PGS000639", "PGS000073", "PGS002302"])
    # accuracy_cdf(["PGS000037", "PGS000639", "PGS000073", "PGS002302"])
    # score_deviation()
    # weight_to_af()
    # likelihood_spectrum_select("PGS000154")
    # likelihood_spectrum_select("PGS002302")
    # likelihood_spectrum_select("PGS000851")
    # likelihood_spectrum_total()
#     loci_coverage()
#     accuracy(["PGS003181", "PGS000778", "PGS004249", "PGS001868", "PGS002270", "PGS001835"])
#     accuracy(["PGS003181", "PGS000778", "PGS004249", "PGS001868", "PGS002270", "PGS001835"])
# 	sequential()
#     king_test()
#     kinship_experiment()
# 	score_uniqueness()
    random_hist()

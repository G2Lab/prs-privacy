import csv
import fnmatch
import json
import random

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import statistics
from scipy.stats import norm


def pairwise(filepath):
    correlations = []
    with open(filepath, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            row_values = [float(value) for value in row]
            correlations.append(row_values)

    xs = np.arange(1, len(correlations) + 1, 1)
    ys = []
    for i in range(len(correlations)):
        ys.append([correlations[i][j] for j in range(i + 1, len(correlations))])
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
                positions.append(1 - position / total_solutions)
                likelihoods.append(likelihood)
                major_likelihoods.append(major_likelihood * score / major_score)

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
    for i in range(0, int(len(labels) / 2)):
        plt.hist(likelihoods[2 * i], bins=8, histtype='step', linewidth=2, color=colors[i],
                 label=f"{labels[i]}({scores[2 * i]})")
        plt.hist(likelihoods[2 * i + 1], bins=8, histtype='step', linewidth=2, color=colors[i],
                 label=f"{labels[i]}({scores[2 * i + 1]})")
        plt.plot(likelihoods[2 * i][0], 1600, 'x', markersize=5, color=colors[i])
        plt.plot(likelihoods[2 * i + 1][0], 1600, 'x', markersize=5, color=colors[i])
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
    plt.scatter([scores[0] for x in range(0, len(distribution[0]))], distribution[0], alpha=transparency,
                color='moccasin',
                s=marker_size, label="Incorrect solutions")
    for i in range(1, len(scores)):
        plt.scatter([scores[i] for x in range(0, len(distribution[i]))], distribution[i], alpha=transparency,
                    color='moccasin',
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
        proximity.append(float(position) * 100 / float(len(distribution[i]) + 1))
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
        filepath = os.path.join(directory, pgs + ".json")
        with open(filepath, 'r') as file:
            for row in file:
                row = row.strip()
                data = json.loads(row)
                accuracy = list(map(lambda x: float(x), data["Accuracies"]))
                positions.append(float(accuracy.index(1.0)) / float(len(accuracy)))

        spos = [pos * 100 for pos in np.sort(positions)]
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
        filepath = os.path.join(directory, pgs + ".json")
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
    full = [left[i] + right[i] for i in range(0, len(left))]
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
    filepath = os.path.join(directory, pgs_id + ".json")
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
        filepath = os.path.join(directory, pgs_id + ".json")
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
    ax1.plot(loci, spectra_median, color=colors[1], marker=".", markersize=8, label="Spectrum")
    ax1.plot(loci, likelihoods_spectra_median, color=colors[2], marker=".", markersize=8, label="Likelihood+Spectrum")
    ax1.plot(loci, reference_median, color='k', linewidth=0.5, linestyle='--', marker=".", markersize=5,
             label="Reference")

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
    ax2.plot(loci, reference_mean, color='k', linewidth=0.5, linestyle='--', marker=".", markersize=5,
             label="Reference")

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
    for pgs_id in pgs_ids + [world]:
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
        if "raccuracies" not in file_name:
            # if "with_repair" not in file_name:
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
#     fig.savefig('sequential-repaired.png', dpi=300, bbox_inches='tight')


def sequential_accuracy():
    directory = "results/guessAccuracy/"
    data = []
    filepath = os.path.join(directory, "accuracy.json")
    with open(filepath, 'r') as f:
        content = json.load(f)
    for result in content:
        individual = result['Individual']
        ancestry = result['Ancestry']
        snps = int(result['SNPs'])
        guess_acc = float(result['GuessAccuracy']) * 100
        reference_acc = float(result['ReferenceAccuracy']) * 100
        data.append({'Individual': individual, 'Ancestry': ancestry, 'SNPs': snps,
                     'GuessAccuracy': guess_acc, 'ReferenceAccuracy': reference_acc})

    df = pd.DataFrame(data)

    # Calculate and print the median number of SNPs per ancestry
    median_snps_per_ancestry = df.groupby('Ancestry')['SNPs'].median()
    print("Median number of SNPs per ancestry:")
    print(median_snps_per_ancestry)
    # Calculate and print the median number of SNPs total
    median_snps_total = df['SNPs'].median()
    print("Median number of SNPs total:", median_snps_total)

    fig, ax = plt.subplots(figsize=(4, 3))
    sns.boxplot(x='Ancestry', y='GuessAccuracy', data=df, palette='Pastel1', hue='Ancestry')

    median_ref_acc = df.groupby('Ancestry')['ReferenceAccuracy'].median()
    for i, (ancestry, median) in enumerate(median_ref_acc.items()):
        ax.plot([i - 0.35, i + 0.35], [median, median], linestyle='--', alpha=0.7, color='black')
        # Add "Baseline" label above the dashed line
        if i == 2:
            ax.text(i, median - 4, 'Baseline', ha='center', va='bottom', fontsize=10, color='black')

    plt.xlabel('Ancestry')
    plt.ylabel('Genotype recovery accuracy, %')
    plt.tight_layout()
    # plt.ylim(64, 100)
    # fig.savefig('recovery.pdf', dpi=300, bbox_inches='tight')
    plt.show()


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
    #     sns.violinplot(x='SNPs', y='Position', ax=ax, data=df)
    #     plt.gca().invert_yaxis()
    sns.boxplot(x='SNPs', y='TruePhi', data=df)
    plt.axhline(y=0.25, linestyle='--', color='r')
    plt.axhline(y=0.125, linestyle='--', color='r')
    plt.xlabel('#SNPs')
    plt.title('Phi of a relative in the KING test')
    plt.ylabel('Phi')
    #     plt.title('Position of a relative by KING test')
    #     plt.ylabel('Position')
    plt.tight_layout()
    #     fig.savefig('king-position.png', dpi=300, bbox_inches='tight')
    fig.savefig('king-phi.png', dpi=300, bbox_inches='tight')


#     plt.show()


def kinship_experiment():
    directory = "results/kinship/"
    num_snps = [2000, 2500]
    #     num_snps = [2000]
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
    #     plt.gca().invert_yaxis()
    #     sns.boxplot(x='SNPs', y='Position', data=df)
    #     sns.boxplot(x='SNPs', y='TruePhi', data=df)
    #     plt.axhline(y=0.25, linestyle='--', color='r')
    #     plt.axhline(y=0.125, linestyle='--', color='r')
    #     plt.title('Phi of a relative by KING test')
    plt.title('Position of a relative by KING test')
    plt.xlabel('#SNPs')
    #     plt.ylabel('Phi')
    plt.ylabel('Position')
    plt.tight_layout()
    #     fig.savefig('kinship-position.png', dpi=300, bbox_inches='tight')
    #     fig.savefig('kinship-phi.png', dpi=300, bbox_inches='tight')
    plt.show()


def score_uniqueness():
    directory = "results/uniqueness/"
    filepath = os.path.join(directory, "scores.json")
    data = []
    with open(filepath, 'r') as f:
        content = json.load(f)
    for result in content:
        data.append({'NumVariants': int(result['NumVariants']),
                     'Unique': 100*float(result['NumUniqueIndividuals'])/2504,
                     'AnonSize': np.median([float(x) for x in result['AnonymitySets']])})

    first_color = 'teal'
    second_color = 'coral'
    df = pd.DataFrame(data)
    fig, ax1 = plt.subplots(figsize=(4, 3))
    sns.lineplot(x='NumVariants', y='AnonSize', data=df, ax=ax1, color=first_color, linestyle='--')
    plt.xlabel('Number of variants in a PRS')
    plt.ylabel('Median anonymity set size', color=first_color)
    ax1.tick_params(axis='y', labelcolor=first_color)

    ax2 = ax1.twinx()
    sns.lineplot(x='NumVariants', y='Unique', data=df, ax=ax2, color=second_color)
    plt.ylabel('Uniquely scored individuals (%)', color=second_color)
    ax2.tick_params(axis='y', labelcolor=second_color)

    grouped = df.groupby('NumVariants')['AnonSize'].median().reset_index()
    for _, row in grouped.iterrows():
        num_variants = row['NumVariants']
        median_anon_size = row['AnonSize']
        print(f'Number of Variants: {num_variants}, Median Anonymity Set Size: {median_anon_size}')

    grouped = df.groupby('NumVariants')['Unique'].median().reset_index()
    for _, row in grouped.iterrows():
        print(f'Number of Variants: {row['NumVariants']}, Unique percentage: {row['Unique']}')

    plt.tight_layout()
    fig.savefig('uniqueness.pdf', dpi=300, bbox_inches='tight')
    # plt.show()


def random_hist():
    params = {
        'axes.labelsize': 36,
        'xtick.labelsize': 36,
        'ytick.labelsize': 36,
        'text.usetex': True,
        'font.family': 'serif'
              }
    mpl.rcParams.update(params)
    data = np.random.beta(0.5, 0.5, 100)

    # Create a histogram using seaborn
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.histplot(data, bins=10, kde=False, color='darkorange')
    plt.xlabel('Population AF', fontweight='bold')
    plt.ylabel('')
    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    plt.xlim(0, 1)
    plt.xticks([0, 0.5, 1])
    # Show the plot
    fig.savefig('allelefreq.pdf', dpi=300, bbox_inches='tight')
    # plt.tight_layout()
    # plt.show()


def af_hist():
    params = {
        'axes.labelsize': 36,
        'xtick.labelsize': 36,
        'ytick.labelsize': 36,
        'text.usetex': True,
        'font.family': 'serif'
    }
    mpl.rcParams.update(params)

    x = [0, 1, 2]
    y = [0.6, 0.3, 0.1]
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.barplot(x=x, y=y, ax=ax, color='darkorange')
    ax.set_xticks([0, 1, 2])
    ax.set_yticks([0, 0.5, 1])
    plt.ylabel('Frequency')
    plt.xlabel('Genotype')
    # plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    # Show the plot
    fig.savefig('allelefreq.pdf', dpi=300, bbox_inches='tight')
    # plt.tight_layout()
    # plt.show()


def guessed_mia():
    directory = "results/impute/"
    filepath = os.path.join(directory, "unimputed.json")
    data = []
    with open(filepath, 'r') as f:
        content = json.load(f)
    for idv, result in content.items():
        data.append({'SelfPos': int(result['SelfPos']), 'RelativePos': int(result['RelativePos']),
                     'SelfAccuracy': float(result['SelfAccuracy']), 'RelativeAccuracy': float(result['RelAccuracy'])})

    df = pd.DataFrame(data)
    fig, ax = plt.subplots(figsize=(4, 3))
    # melted_df = pd.melt(df, value_vars=['SelfPos', 'RelativePos'], var_name='PositionType', value_name='Value')
    # sns.boxplot(x='PositionType', y='Value', data=melted_df)
    # plt.title("Linking position")
    # plt.ylabel('Match position')
    # fig.savefig('linking-pos.png', dpi=300, bbox_inches='tight')
    # plt.gca().invert_yaxis()
    melted_df = pd.melt(df, value_vars=['SelfAccuracy', 'RelativeAccuracy'], var_name='Type', value_name='Value')
    sns.boxplot(x='Type', y='Value', data=melted_df)
    plt.title("SNP matching accuracy")
    plt.ylabel('Match comparison accuracy')
    fig.savefig('linking-acc.png', dpi=300, bbox_inches='tight')
    plt.xlabel('')
    plt.tight_layout()
    plt.show()


def full_imputed():
    directory = "results/impute/"
    filepath = os.path.join(directory, "full_imputed22.json")
    data = []
    with open(filepath, 'r') as f:
        content = json.load(f)
    for idv, result in content.items():
        data.append({'SelfPosCount': int(result['Self']['count']), 'SelfPosKing': int(result['Self']['king']),
                     'RelativePosCount': int(result['Relative']['count']),
                     'RelativePosKing': int(result['Relative']['king']),
                     'SelfAccuracy': float(result['CountAccuracy']["self"]),
                     'RelativeAccuracy': float(result['KingScore']["relative"])})

    df = pd.DataFrame(data)
    fig, ax = plt.subplots(figsize=(3, 3))
    # melted_df = pd.melt(df, value_vars=['SelfPosCount', 'SelfPosKing'], var_name='Method', value_name='Value')
    # melted_df = pd.melt(df, value_vars=['RelativePosCount', 'RelativePosKing'], var_name='Method', value_name='Value')
    sns.boxplot(x=["Relative" for x in df['RelativeAccuracy']], y='RelativeAccuracy', data=df)
    # plt.title("Self linking")
    # plt.title("Relative linking")
    plt.title("Accuracy")
    plt.ylabel('KING score of the relative')
    # plt.gca().invert_yaxis()
    # ax.set_xticklabels(['Self match ratio', 'Relative KING score'])
    plt.xlabel('')
    plt.tight_layout()
    fig.savefig('relative-accuracy.png', dpi=300, bbox_inches='tight')
    # fig.savefig('self-linking.png', dpi=300, bbox_inches='tight')
    # melted_df = pd.melt(df, value_vars=['SelfAccuracy', 'RelativeAccuracy'], var_name='Type', value_name='Value')
    # sns.boxplot(x='Type', y='Value', data=melted_df)
    # plt.title("SNP matching accuracy")
    # plt.ylabel('Match comparison accuracy')
    # fig.savefig('linking-acc.png', dpi=300, bbox_inches='tight')
    plt.show()


def linking_accuracy():
    directory = "results/linking/"
    filepaths = [os.path.join(directory, filename) for filename in os.listdir(directory) if
                 # fnmatch.fnmatch(filename, f"unimputed*.json")]
                 fnmatch.fnmatch(filename, f"imputed*.json")]
    data = {}
    for filepath in filepaths:
        with open(filepath, 'r') as f:
            content = json.load(f)
        for idv, relations in content.items():
            if idv not in data:
                data[idv] = {}
            for rel in relations:
                if rel["Target"] not in data[idv]:
                    data[idv][rel["Target"]] = {}
                    data[idv][rel["Target"]]["Count"] = [float(c) for c in rel["Count"]]
                    data[idv][rel["Target"]]["King"] = [float(c) for c in rel["King"]]
                    data[idv][rel["Target"]]["Information"] = [float(c) for c in rel["Information"] if float(c) < 1e20]
                else:
                    for i, c in enumerate(rel["Count"]):
                        data[idv][rel["Target"]]["Count"][i] += float(c)
                    for i, c in enumerate(rel["King"]):
                        data[idv][rel["Target"]]["King"][i] += float(c)
                    for i, c in enumerate(rel["Information"]):
                        data[idv][rel["Target"]]["Information"][i] += float(c) if float(c) < 1e20 else 0

    divided = {}
    for idv in data:
        divided[idv] = {"Count": [], "King": [], "Information": []}
        for target in data[idv]:
            divided[idv]["Count"].append((target, data[idv][target]["Count"][0] / data[idv][target]["Count"][1]))
            divided[idv]["King"].append((target, data[idv][target]["King"][0] / data[idv][target]["King"][1]))
            divided[idv]["Information"].append(
                (target, data[idv][target]["Information"][0] / data[idv][target]["Information"][1]))
        # sort by the second element of the tuple
        divided[idv]["Count"] = sorted(divided[idv]["Count"], key=lambda x: x[1], reverse=True)
        divided[idv]["King"] = sorted(divided[idv]["King"], key=lambda x: x[1], reverse=True)
        divided[idv]["Information"] = sorted(divided[idv]["Information"], key=lambda x: x[1], reverse=True)
        # print(f"IDV: {idv}, Count: {divided[idv]['Count'][0]}, King: {divided[idv]['King'][0]}")

    parsed = {"SelfCountPos": [], "SelfKingPos": [], "SelfInformationPos": [],
              "RelativeCountPos": [], "RelativeKingPos": [], "RelativeInformationPos": [],
              "ReferenceCountPos": [], "ReferenceKingPos": [], "ReferenceInformationPos": [],
              "SelfCountAcc": [], "SelfKingAcc": [], "SelfInformationAcc": [],
              "RelativeCountAcc": [], "RelativeKingAcc": [], "RelativeInformationAcc": [],
              "ReferenceCountAcc": [], "ReferenceKingAcc": [], "ReferenceInformationAcc": [],
              "EveryCountRatio": [], "EveryKingRatio": [], "EveryInformationRatio": []}
    related = read_related_individuals()
    for idv in divided:
        relative_found = False
        for i, tpl in enumerate(divided[idv]["Count"]):
            if tpl[0] == idv:
                parsed["SelfCountPos"].append(i)
                parsed["SelfCountAcc"].append(tpl[1])
            if tpl[0] == "reference":
                parsed["ReferenceCountPos"].append(i)
                parsed["ReferenceCountAcc"].append(tpl[1])
            if not relative_found and tpl[0] in related[idv]:
                parsed["RelativeCountPos"].append(i)
                parsed["RelativeCountAcc"].append(tpl[1])
                relative_found = True
            parsed["EveryCountRatio"].append(tpl[1])
    for idv in divided:
        relative_found = False
        for i, tpl in enumerate(divided[idv]["King"]):
            if tpl[0] == idv:
                parsed["SelfKingPos"].append(i)
                parsed["SelfKingAcc"].append(tpl[1])
            if tpl[0] == "reference":
                parsed["ReferenceKingPos"].append(i)
                parsed["ReferenceKingAcc"].append(tpl[1])
            if not relative_found and tpl[0] in related[idv]:
                parsed["RelativeKingPos"].append(i)
                parsed["RelativeKingAcc"].append(tpl[1])
                relative_found = True
            parsed["EveryKingRatio"].append(tpl[1])
    for idv in divided:
        relative_found = False
        for i, tpl in enumerate(divided[idv]["Information"]):
            if tpl[0] == idv:
                parsed["SelfInformationPos"].append(i)
                parsed["SelfInformationAcc"].append(tpl[1])
            if tpl[0] == "reference":
                parsed["ReferenceInformationPos"].append(i)
                parsed["ReferenceInformationAcc"].append(tpl[1])
            if not relative_found and tpl[0] in related[idv]:
                parsed["RelativeInformationPos"].append(i)
                parsed["RelativeInformationAcc"].append(tpl[1])
                relative_found = True
            parsed["EveryInformationRatio"].append(tpl[1])

    print(f"SelfCountAcc: {np.median(parsed['SelfCountAcc'])}, {np.mean(parsed['SelfCountAcc'])}")
    print(f"RelativeCountAcc: {np.median(parsed['RelativeCountAcc'])}, {np.mean(parsed['RelativeCountAcc'])}")
    print(f"ReferenceCountAcc: {np.median(parsed['ReferenceCountAcc'])}, {np.mean(parsed['ReferenceCountAcc'])}")
    print(f"SelfKingAcc: {np.median(parsed['SelfKingAcc'])}, {np.mean(parsed['SelfKingAcc'])}")
    print(f"RelativeKingAcc: {np.median(parsed['RelativeKingAcc'])}, {np.mean(parsed['RelativeKingAcc'])}")
    print(f"ReferenceKingAcc: {np.median(parsed['ReferenceKingAcc'])}, {np.mean(parsed['ReferenceKingAcc'])}")
    print(f"SelfInformationAcc: {np.median(parsed['SelfInformationAcc'])}, {np.mean(parsed['SelfInformationAcc'])}")
    print(
        f"RelativeInformationAcc: {np.median(parsed['RelativeInformationAcc'])}, {np.mean(parsed['RelativeInformationAcc'])}")
    print(
        f"ReferenceInformationAcc: {np.median(parsed['ReferenceInformationAcc'])}, {np.mean(parsed['ReferenceInformationAcc'])}")

    palette = {}
    for key in parsed:
        if key.startswith("Self"):
            palette[key] = sns.color_palette("pastel")[0]
        elif key.startswith("Relative"):
            palette[key] = sns.color_palette("pastel")[1]
        elif key.startswith("Reference"):
            palette[key] = sns.color_palette("pastel")[2]
        else:
            palette[key] = sns.color_palette("pastel")[3]

    keys = [["SelfCountPos", "RelativeCountPos", "ReferenceCountPos"],
            ["SelfKingPos", "RelativeKingPos", "ReferenceKingPos"],
            ["SelfInformationPos", "RelativeInformationPos", "ReferenceInformationPos"],
            ["SelfCountAcc", "RelativeCountAcc", "ReferenceCountAcc", "EveryCountRatio"],
            ["SelfKingAcc", "RelativeKingAcc", "ReferenceKingAcc", "EveryKingRatio"],
            ["SelfInformationAcc", "RelativeInformationAcc", "ReferenceInformationAcc", "EveryInformationRatio"]]
    df = []
    for key in keys:
        df.append(prepare_data(parsed, key))
    fig1, axes1 = plt.subplots(1, 3, figsize=(18, 5))
    for i, ax1 in enumerate(axes1):
        ax1.invert_yaxis()
        ax1.set_ylabel('Rank')
        ax1.set_xlabel('')
        ax1.set_xticklabels(["Self", "Relative", "Reference"])
        sns.boxplot(x='Versus', y='Value', data=df[i], hue='Versus', palette=palette, ax=axes1[i])
        medians = df[i].groupby(['Versus'])['Value'].median()
        print(f"Medians for df[{i}]:", medians)  # Debug print
        # for tick, label in zip(range(len(medians)), ax1.get_xticklabels()):
        #     ax1.annotate(f'{medians[label.get_text()]:.2f}',
        #                  xy=(tick, medians[label.get_text()]),
        #                  xytext=(0, 5),  # 5 points vertical offset
        #                  textcoords='offset points',
        #                  ha='center', va='center',
        #                  color='red', fontsize=10, fontweight='bold')
    axes1[0].set_title('Match Count Based')
    axes1[0].set_ylim(2550, -50)
    axes1[1].set_title('KING Based')
    axes1[1].set_ylim(2550, -50)
    axes1[2].set_title('Mutual Information Based')
    axes1[2].set_ylim(2550, -50)
    # fig1.savefig('unimputed-pos.png', dpi=300, bbox_inches='tight')
    fig1.savefig('imputed-pos.png', dpi=300, bbox_inches='tight')

    fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))
    for i, ax2 in enumerate(axes2):
        ax2.set_xlabel('')
        ax2.set_xticklabels(["Self", "Relative", "Reference", "All"])
        sns.boxplot(x='Versus', y='Value', data=df[i + 3], hue='Versus', palette=palette, ax=axes2[i])
        medians = df[i + 3].groupby(['Versus'])['Value'].median()
        print(f"Medians for df[{i}]:", medians)  # Debug print
        # for tick, label in zip(range(len(medians)), ax2.get_xticklabels()):
        #     ax2.annotate(f'{medians[label.get_text()]:.2f}',
        #                  xy=(tick, medians[label.get_text()]),
        #                  xytext=(0, 5),  # 5 points vertical offset
        #                  textcoords='offset points',
        #                  ha='center', va='center',
        #                  color='red', fontsize=10, fontweight='bold')
    axes2[0].set_title('Match Count Based')
    axes2[0].set_ylabel('Match ratio')
    axes2[1].set_ylabel('KING score')
    axes2[1].set_title('KING Based')
    axes2[2].set_ylabel('Mutual information')
    axes2[2].set_title('Mutual Information Based')
    # fig2.savefig('unimputed-acc.png', dpi=300, bbox_inches='tight')
    fig2.savefig('imputed-acc.png', dpi=300, bbox_inches='tight')

    plt.tight_layout()
    plt.show()


def king_accuracy():
    directory = "results/king/"
    filepath = os.path.join(directory, "relations.json")
    # filepaths = [os.path.join(directory, filename) for filename in os.listdir(directory) if
    #              fnmatch.fnmatch(filename, f"unimputed*.json")]

    data = {}
    ancestries = {}
    with open(filepath, 'r') as f:
        content = json.load(f)
    for idv, entry in content.items():
        if idv not in data:
            data[idv] = {}
        ancestries[idv] = entry["Ancestry"]
        for target, relation in entry["Relations"].items():
            if target not in data[idv]:
                data[idv][target] = {}
                data[idv][target]["Count"] = [float(c) for c in relation["Count"]]
                data[idv][target]["King"] = [float(c) for c in relation["King"]]
            else:
                for i, c in enumerate(relation["Count"]):
                    data[idv][target]["Count"][i] += float(c)
                for i, c in enumerate(relation["King"]):
                    data[idv][target]["King"][i] += float(c)

    divided = {}
    for idv in data:
        divided[idv] = {"Count": [], "King": []}
        for target in data[idv]:
            divided[idv]["Count"].append((target, data[idv][target]["Count"][0] / data[idv][target]["Count"][1]))
            divided[idv]["King"].append((target, data[idv][target]["King"][0] / data[idv][target]["King"][1]))
        # sort by the match count and king score
        divided[idv]["Count"] = sorted(divided[idv]["Count"], key=lambda x: x[1], reverse=True)
        divided[idv]["King"] = sorted(divided[idv]["King"], key=lambda x: x[1], reverse=True)
        # print(f"IDV: {idv}, Count: {divided[idv]['Count'][0]}, King: {divided[idv]['King'][0]}")

    parsed = {"SelfCountPos": [], "SelfKingPos": [],
              "RelativeCountPos": [], "RelativeKingPos": [],
              "ReferenceCountPos": [], "ReferenceKingPos": [],
              "SelfCountAcc": [], "SelfKingAcc": [],
              "RelativeCountAcc": [], "RelativeKingAcc": [],
              "ReferenceCountAcc": [], "ReferenceKingAcc": [],
              "EveryCountRatio": [], "EveryKingRatio": []}
    related = read_related_individuals()
    for idv in divided:
        relative_found = False
        for i, tpl in enumerate(divided[idv]["Count"]):
            if tpl[0] == idv:
                parsed["SelfCountPos"].append(i)
                parsed["SelfCountAcc"].append(tpl[1])
            if tpl[0] == "reference":
                parsed["ReferenceCountPos"].append(i)
                parsed["ReferenceCountAcc"].append(tpl[1])
            if not relative_found and idv in related and tpl[0] in related[idv]:
                parsed["RelativeCountPos"].append(i)
                parsed["RelativeCountAcc"].append(tpl[1])
                relative_found = True
            parsed["EveryCountRatio"].append(tpl[1])
    for idv in divided:
        relative_found = False
        for i, tpl in enumerate(divided[idv]["King"]):
            if tpl[0] == idv:
                parsed["SelfKingPos"].append(i)
                parsed["SelfKingAcc"].append(tpl[1])
            if tpl[0] == "reference":
                parsed["ReferenceKingPos"].append(i)
                parsed["ReferenceKingAcc"].append(tpl[1])
            if not relative_found and idv in related and tpl[0] in related[idv]:
                parsed["RelativeKingPos"].append(i)
                parsed["RelativeKingAcc"].append(tpl[1])
                relative_found = True
            parsed["EveryKingRatio"].append(tpl[1])

    print(f"SelfCountAcc: {np.median(parsed['SelfCountAcc'])}, {np.mean(parsed['SelfCountAcc'])}")
    print(f"RelativeCountAcc: {np.median(parsed['RelativeCountAcc'])}, {np.mean(parsed['RelativeCountAcc'])}")
    print(f"ReferenceCountAcc: {np.median(parsed['ReferenceCountAcc'])}, {np.mean(parsed['ReferenceCountAcc'])}")
    print(f"EveryCountRatio: {np.median(parsed['EveryCountRatio'])}, {np.mean(parsed['EveryCountRatio'])}")
    print(f"SelfKingAcc: {np.median(parsed['SelfKingAcc'])}, {np.mean(parsed['SelfKingAcc'])}")
    print(f"RelativeKingAcc: {np.median(parsed['RelativeKingAcc'])}, {np.mean(parsed['RelativeKingAcc'])}")
    print(f"ReferenceKingAcc: {np.median(parsed['ReferenceKingAcc'])}, {np.mean(parsed['ReferenceKingAcc'])}")
    print(f"EveryKingRatio: {np.median(parsed['EveryKingRatio'])}, {np.mean(parsed['EveryKingRatio'])}")

    palette = {}
    for key in parsed:
        if key.startswith("Self"):
            palette[key] = sns.color_palette("pastel")[0]
        elif key.startswith("Relative"):
            palette[key] = sns.color_palette("pastel")[1]
        elif key.startswith("Reference"):
            palette[key] = sns.color_palette("pastel")[2]
        else:
            palette[key] = sns.color_palette("pastel")[3]

    keys = [["SelfCountPos", "RelativeCountPos", "ReferenceCountPos"],
            ["SelfKingPos", "RelativeKingPos", "ReferenceKingPos"],
            ["SelfCountAcc", "RelativeCountAcc", "ReferenceCountAcc", "EveryCountRatio"],
            ["SelfKingAcc", "RelativeKingAcc", "ReferenceKingAcc", "EveryKingRatio"]]
    df = []
    for key in keys:
        df.append(prepare_data(parsed, key))
    fig1, axes1 = plt.subplots(1, 2, figsize=(12, 5))
    for i, ax1 in enumerate(axes1):
        ax1.invert_yaxis()
        ax1.set_ylabel('Rank')
        ax1.set_xlabel('')
        ax1.set_xticklabels(["Self", "Relative", "Reference"])
        sns.boxplot(x='Versus', y='Value', data=df[i], hue='Versus', palette=palette, ax=axes1[i])
        medians = df[i].groupby(['Versus'])['Value'].median()
        print(f"Medians for df[{i}]:", medians)  # Debug print
        # for tick, label in zip(range(len(medians)), ax1.get_xticklabels()):
        #     ax1.annotate(f'{medians[label.get_text()]:.2f}',
        #                  xy=(tick, medians[label.get_text()]),
        #                  xytext=(0, 5),  # 5 points vertical offset
        #                  textcoords='offset points',
        #                  ha='center', va='center',
        #                  color='red', fontsize=10, fontweight='bold')
    axes1[0].set_title('Match Count Based')
    axes1[0].set_ylim(2550, -50)
    axes1[1].set_title('KING Based')
    axes1[1].set_ylim(2550, -50)
    # fig1.savefig('unimputed-pos.png', dpi=300, bbox_inches='tight')
    # fig1.savefig('imputed-pos.png', dpi=300, bbox_inches='tight')

    fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
    for i, ax2 in enumerate(axes2):
        ax2.set_xlabel('')
        ax2.set_xticklabels(["Self", "Relative", "Reference", "All"])
        sns.boxplot(x='Versus', y='Value', data=df[i + int(len(df)/2)], hue='Versus', palette=palette, ax=axes2[i])
        medians = df[i + int(len(df)/2)].groupby(['Versus'])['Value'].median()
        print(f"Medians for df[{i}]:", medians)  # Debug print
        # for tick, label in zip(range(len(medians)), ax2.get_xticklabels()):
        #     ax2.annotate(f'{medians[label.get_text()]:.2f}',
        #                  xy=(tick, medians[label.get_text()]),
        #                  xytext=(0, 5),  # 5 points vertical offset
        #                  textcoords='offset points',
        #                  ha='center', va='center',
        #                  color='red', fontsize=10, fontweight='bold')
    axes2[0].set_title('Match Count Based')
    axes2[0].set_ylabel('Match ratio')
    axes2[1].set_ylabel('KING score')
    axes2[1].set_title('KING Based')
    # fig2.savefig('unimputed-acc.png', dpi=300, bbox_inches='tight')
    # fig2.savefig('imputed-acc.png', dpi=300, bbox_inches='tight')

    plt.tight_layout()
    plt.show()


def prepare_data(parsed, keys):
    data = []
    for key in keys:
        data.extend([(key, value) for value in parsed[key]])
    return pd.DataFrame(data, columns=['Versus', 'Value'])


def prepare_data_for_cdf(parsed, key):
    data = []
    for label in key:
        values = parsed[label]
        data.extend([(label, value) for value in values])
    return data


def read_related_individuals():
    file_path = "data/related_individuals.txt"
    related = {}

    try:
        with open(file_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            header = next(reader)

            sample_column = -1
            relatives_column = -1

            for i, field in enumerate(header):
                if field == "Sample":
                    sample_column = i
                if field == "Reason for exclusion":
                    relatives_column = i

            for record in reader:
                if len(record) > relatives_column:
                    relation, relatives = record[relatives_column].split(":")
                    if relation == "Parent" or relation == "Sibling" or relation == "trio":
                        relation = "1"
                    elif relation == "Second Order":
                        relation = "2"
                    related[record[sample_column]] = [(relative, relation) for relative in relatives.split(",")]
                    for relative in relatives.split(","):
                        related[relative] = [(record[sample_column], relation)]

    except Exception as e:
        print(f"Error: {e}")

    return related


def imputation_accuracy():
    directory = "results/impute/f1/"
    filepaths = [os.path.join(directory, filename) for filename in os.listdir(directory) if
                 fnmatch.fnmatch(filename, f"*.json")]
    data = []
    for filepath in filepaths:
        with open(filepath, 'r') as f:
            content = json.load(f)
        for distance, accuracies in content.items():
            if distance == "0":
                continue
            # if int(distance) > 10000:
            #     continue
            for acc in accuracies:
                data.append((int(distance), float(acc)))

    df = pd.DataFrame(data, columns=["Distance", "Accuracy"])
    range_size = 1000
    df['DistanceRange'] = (df['Distance'] // range_size) * range_size
    mean_accuracy_by_range = df.groupby('DistanceRange')['Accuracy'].mean().reset_index()
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(x='DistanceRange', y='Accuracy', data=mean_accuracy_by_range)
    plt.xlabel('Distance')
    plt.ylabel('F1 Accuracy score')
    plt.title('Imputation accuracy by distance from a known SNP')
    ax.xaxis.set_major_locator(plt.MaxNLocator(12))
    # xticks = mean_accuracy_by_range['DistanceRange']
    # ax.set_xticks(xticks)
    # ax.set_xticklabels([str(tick) if i % 10000 == 0 else '' for i, tick in enumerate(xticks)])
    plt.show()
    # plt.savefig('imputation.png', dpi=300, bbox_inches='tight')


def percentile_prediction():
    filepath = "results/predict/prediction.json"
    data = []
    with open(filepath, 'r') as f:
        content = json.load(f)
    for result in content.values():
        for i, predicted in enumerate(result['Predicted']):
            data.append({'Known': f'{result["Known"]}/{result["Total"]}',
                         'Percentile Difference': 100 * abs(float(result['Predicted'][i]) - float(result['Real'][i]))})

    sorted(data, key=lambda x: float(x['Known'].split("/")[0]) / float(x['Known'].split("/")[1]))
    df = pd.DataFrame(data)
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.invert_yaxis()
    sns.boxplot(x='Known', y='Percentile Difference', data=df)
    plt.title('Percentile difference between predicted and real PRS')
    plt.xlabel('Known/Total SNPs')
    plt.tight_layout()
    fig.savefig('prediction.png', dpi=300, bbox_inches='tight')
    # plt.show()


def random_bars():
    # Data
    labels = ['TP53', 'BRCA1', 'CFTR']
    values1 = [10, 15, 20]
    values2 = [x + random.randint(1, 5) for x in values1]  # Slightly higher values

    # Parameters
    x = np.arange(len(labels))  # the label locations
    width = 0.22  # narrower width of the bars

    # Set seaborn color palette
    palette = sns.color_palette("colorblind")

    fig, ax = plt.subplots(figsize=(5, 4))

    # Plotting the bars with seaborn colors
    bars1 = ax.bar(x - width / 2, values1, width, label='Old technique', color=palette[2])
    bars2 = ax.bar(x + width / 2, values2, width, label='My awesome technique', color=palette[1])

    # Adding labels, title, and legend
    ax.set_xlabel('Genes')
    ax.set_ylabel('Predictive power')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    fig.savefig('randbars.pdf', dpi=300, bbox_inches='tight')
    # plt.show()


def plot_normal_distribution_with_fill(mu=0, sigma=1):
    params = {
        'axes.labelsize': 36,
        'xtick.labelsize': 36,
        'ytick.labelsize': 36,
        'text.usetex': True,
        'font.family': 'serif'
    }
    mpl.rcParams.update(params)

    x = np.linspace(mu - 4*sigma, mu + 4*sigma, 1000)
    y = norm.pdf(x, mu, sigma)
    plt.figure(figsize=(4, 3))
    plt.plot(x, y, 'k')  # 'k' for black color
    # Calculate the cutoff value for the top 10%
    cutoff = norm.ppf(0.85, mu, sigma)

    x_fill = np.linspace(cutoff, mu + 4*sigma, 1000)
    y_fill = norm.pdf(x_fill, mu, sigma)
    plt.fill_between(x_fill, y_fill, color='red', alpha=0.5)
    plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                   labelbottom=False, labelleft=False)

    # Add labels and legend
    plt.xlabel('Risk score')
    plt.ylabel('Frequency')
    plt.savefig('risk_distro.pdf', dpi=300, bbox_inches='tight')
    # Show the plot
    # plt.tight_layout()
    # plt.show()


def ibd_accuracy():
    indf = pd.read_csv("ibd/refined.ibd", sep='\t', index_col=0)
    pruned_df = indf.loc[indf.index.str.startswith('$'), ~indf.columns.str.startswith('$')]
    results = []
    for idv in pruned_df.index:
        target = idv[1:]
        if target not in pruned_df.columns:
            continue
        kinship = pruned_df.at[idv, target]
        all_kinship_values = pruned_df[target]
        rank = all_kinship_values.rank(ascending=False).loc[idv]
        results.append({'ID': target, 'Kinship': kinship, 'Rank': rank})
        print(f"ID: {target}, Kinship: {kinship}, Rank: {rank}, Top: {all_kinship_values.idxmax()}")

    df = pd.DataFrame(results)
    df['Self'] = 'Self'
    plt.figure(figsize=(4, 3))
    print(f'Median rank: {df['Rank'].median()}')
    # sns.boxplot(x='Self', y='Rank', data=df)
    sns.boxplot(x='Self', y='Kinship', data=df)
    # plt.xticks(rotation=90)  # Rotate x labels for better readability
    plt.tight_layout()
    plt.show()


def plink_accuracy():
    # Load IDs from a separate file
    # with open("ibd/plink.rel.id", "r") as f:
    with open("ibd/plink.king.id", "r") as f:
        ids = [line.strip() for line in f]
    ids = ids[1:]  # Remove the first line

    # Load the kinship matrix
    # indf = pd.read_csv("ibd/plink.rel", sep='\t', header=None)
    indf = pd.read_csv("ibd/plink.king", sep='\t', header=None)
    print(len(ids), indf.shape)

    # Use the IDs to rename columns and index of the DataFrame
    indf.columns = ids
    indf.index = ids

    # Prune the DataFrame to include only rows starting with '$' and columns not starting with '$'
    pruned_df = indf.loc[indf.index.str.startswith('$'), ~indf.columns.str.startswith('$')]

    results = []
    for idv in pruned_df.index:
        target = idv[1:]  # Remove the '$' symbol
        if target not in pruned_df.columns:
            continue
        kinship = pruned_df.at[idv, target]
        all_kinship_values = pruned_df[target]
        rank = all_kinship_values.rank(ascending=False).loc[idv]
        results.append({'ID': target, 'Kinship': kinship, 'Rank': rank})
        print(f"ID: {target}, Kinship: {kinship}, Rank: {rank}, Top: {all_kinship_values.idxmax()}")

    df = pd.DataFrame(results)
    df['Self'] = 'Self'
    plt.figure(figsize=(4, 3))
    print(f'Median rank: {df["Rank"].median()}, Median kinship: {df["Kinship"].median()}')
    print(f'Max rank: {df["Rank"].max()}, Max kinship: {df["Kinship"].max()}')
    print(f'Min rank: {df["Rank"].min()}, Min kinship: {df["Kinship"].min()}')
    # plt.ylim(0, df['Kinship'].max() + 0.1)
    sns.boxplot(x='Self', y='Rank', data=df)
    # sns.boxplot(x='Self', y='Kinship', data=df)
    # plt.xticks(rotation=90)  # Rotate x labels for better readability
    plt.tight_layout()
    plt.show()


def load_results(id_file, data_file, method):
    if id_file is not None:
        with open(id_file, "r") as f:
            ids = [line.strip() for line in f]
        ids = ids[1:]
        indf = pd.read_csv(data_file, sep='\t', header=None)
        indf.columns = ids
        indf.index = ids
    else:
        indf = pd.read_csv("ibd/refined.ibd", sep='\t', index_col=0)

    # Prune the DataFrame to include only rows starting with '$' and columns not starting with '$'
    pruned_df = indf.loc[indf.index.str.startswith('$'), ~indf.columns.str.startswith('$')]
    related = read_related_individuals()

    results = []
    for idv in pruned_df.index:
        target = idv[1:]  # Remove the '$' symbol
        if target not in pruned_df.columns:
            continue
        kinship = pruned_df.at[idv, target]
        self_rank = pruned_df.loc[idv].rank(ascending=False)[target]
        results.append({'Category': 'Self', 'Kinship': kinship, 'Rank': self_rank, 'Method': method})
        relatives = []
        if target in related:
            relatives = []
            for relative in related[target]:
                if relative[0] not in pruned_df.columns:
                    continue
                relative_kinship = pruned_df.at[idv, relative[0]]
                relative_rank = pruned_df.loc[idv].rank(ascending=False)[relative[0]]
                category = ''
                if relative[1] == '1':
                    category = '1st degree'
                elif relative[1] == '2':
                    category = '2nd degree'
                results.append({'Category': category, 'Kinship': relative_kinship, 'Rank': relative_rank, 'Method': method})
                relatives.append(relative[0])
        other_columns = pruned_df.columns[~pruned_df.columns.isin([target] + relatives)]
        kinships = pruned_df.loc[idv, other_columns]
        ranks = pruned_df.loc[idv].rank(ascending=False)[other_columns]
        # high_kinship_positions = kinships[kinships > 0.22]
        # for col in high_kinship_positions.index:
        #     print(f"{method}: {target}->{col} {high_kinship_positions[col]}, rank: {ranks[col]}")
        temp_df = pd.DataFrame({
            'Category': 'Everyone else',
            'Kinship': kinships.values,
            'Rank': ranks.values,
            'Method': method
        })
        results.extend(temp_df.to_dict('records'))

    return pd.DataFrame(results)


def plot_rank_cdf(df, category, ax):
    colors = sns.color_palette("deep")
    # Filter the dataframe for the given category
    category_df = df[df['Category'] == category]
    for i, method in enumerate(category_df['Method'].unique()):
        method_df = category_df[category_df['Method'] == method]
        # Sort the ranks for the cumulative distribution
        sorted_ranks = np.sort(method_df['Rank'])
        # Calculate the cumulative distribution
        cdf = np.arange(1, len(sorted_ranks) + 1) / len(sorted_ranks) * 100
        ax.plot(sorted_ranks, cdf, label=f'{method}', color=colors[i], linewidth=2)
    ax.set_title(f'{category}')
    ax.set_xlim(-5, 300)


def deanonymization_accuracy():
    king_df = load_results("ibd/plink.king.id", "ibd/plink.king", 'KING')
    rel_df = load_results("ibd/plink.rel.id", "ibd/plink.rel", 'Plink IBD')
    rel_df['Kinship'] = rel_df['Kinship'] / 2 # Relatedness -> Kinship
    rel_df['Kinship'] = [0.5 if k > 0.5 else k for k in rel_df['Kinship']]
    refined_df = load_results(None, "ibd/refined.ibd", 'Refined IBD')
    print(f'Median KING rank: {king_df["Rank"].median()}, Median KING kinship: {king_df["Kinship"].median()}')
    print(f'Median Plink IBD rank: {rel_df["Rank"].median()}, Median Plink IBD kinship: {rel_df["Kinship"].median()}')
    print(f'Median Refined IBD rank: {refined_df["Rank"].median()}, Median Refined IBD kinship: {refined_df["Kinship"].median()}')
    df = pd.concat([king_df, rel_df, refined_df])
    # print(df[df['Category'] == '2nd degree'])

    params = {
        'font.size': 13,
        'axes.labelsize': 13,
        'xtick.labelsize': 13,
        'ytick.labelsize': 13,
        'legend.fontsize': 12,
        'legend.title_fontsize': 12,
        'text.usetex': True,
        'font.family': 'serif'
    }
    mpl.rcParams.update(params)
    self_cutoff = 0.3535
    first_degree_cutoff = 0.1768
    second_degree_cutoff = 0.089
    fig1, ax1 = plt.subplots(figsize=(7, 3))
    sns.boxplot(x='Category', y='Kinship', data=df, hue='Method', palette='pastel', ax=ax1,
                order=['Self', '1st degree', '2nd degree', 'Everyone else'])
    ax1.text(x=ax1.get_xlim()[1] + 0.03, y=0.5, s='Criteria', va='center', ha='left', color='black')
    ax1.axhline(y=self_cutoff, color='lightgray', linestyle='--', linewidth=1)
    ax1.text(x=ax1.get_xlim()[1] + 0.03, y=self_cutoff, s='Self', va='center', ha='left', color='black')
    ax1.axhline(y=first_degree_cutoff, color='lightgray', linestyle='--', linewidth=1)
    ax1.text(x=ax1.get_xlim()[1] + 0.03, y=first_degree_cutoff, s='1st', va='center', ha='left', color='black')
    ax1.axhline(y=second_degree_cutoff, color='lightgray', linestyle='--', linewidth=1)
    ax1.text(x=ax1.get_xlim()[1] + 0.03, y=second_degree_cutoff, s='2nd', va='center', ha='left', color='black')
    ax1.set_ylabel('Kinship')
    ax1.set_xlabel('')
    ax1.legend(title="")

    fig2, ax2 = plt.subplots(figsize=(5, 3))
    categories = ['Self', '1st degree', '2nd degree']
    sns.boxplot(x='Category', y='Rank', data=df[df['Category'].isin(categories)], hue='Method', palette='pastel',
                ax=ax2, order=categories)
    ax2.invert_yaxis()
    ax2.set_ylabel('Rank (/2535)')
    ax2.set_xlabel('')
    ax2.legend(title="")

    fig3, axes3 = plt.subplots(1, 3, figsize=(6, 3), sharey=True)
    # Plot CDFs for each category
    categories = ['Self', '1st degree', '2nd degree']
    for i, category in enumerate(categories):
        plot_rank_cdf(df, category, axes3[i])
    axes3[0].set_ylabel('Percentage')
    axes3[1].set_xlabel('Rank (/2535)')
    handles, labels = axes3[0].get_legend_handles_labels()
    fig3.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.05), title='', ncol=len(labels), frameon=False)
    # plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.tight_layout()
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
    # loci_coverage()
    # accuracy(["PGS003181", "PGS000778", "PGS004249", "PGS001868", "PGS002270", "PGS001835"])
    # accuracy(["PGS003181", "PGS000778", "PGS004249", "PGS001868", "PGS002270", "PGS001835"])
    # king_test()
    # kinship_experiment()
    # random_hist()
    # guessed_mia()
    # full_imputed()
    # linking_accuracy()
    # imputation_accuracy()
    # percentile_prediction()
    # random_bars()
    # score_uniqueness()
    # sequential()
    # sequential_accuracy()
    # af_hist()
    # plot_normal_distribution_with_fill()
    # king_accuracy()
    # ibd_accuracy()
    # plink_accuracy()
    # load_results("ibd/plink.rel.id", "ibd/plink.rel")
    deanonymization_accuracy()
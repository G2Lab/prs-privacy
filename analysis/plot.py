import csv
import fnmatch
import json
import matplotlib.pyplot as plt
import numpy as np
import os
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
                 fnmatch.fnmatch(filename, f"*{pgs_id}-s-eaf*")]
    scores, accuracies, likelihoods = [], [], []
    chosen_acc, chosen_like = [], []
    true_like = []
    lowest_index = 0
    lowest_like, lowest_like_acc = [], []
    for filepath in filepaths:
        with open(filepath, 'r') as file:
            for row in file:
                row = row.strip()
                data = json.loads(row)
                scores.append(float(data["Score"]))
                chosen_acc.append(float(data["Accuracies"][0]))
                chosen_like.append(float(data["Likelihoods"][0]))
                acc = list(map(lambda x: float(x), data["Accuracies"]))
                true_like.append(float(data["Likelihoods"][acc.index(1.0)]))
                lowest_index = acc.index(min(acc))
                lowest_like.append(float(data["Likelihoods"][lowest_index]))
                lowest_like_acc.append(float(data["Accuracies"][lowest_index]))
                # if len(data["Accuracies"]) > 1:
                #     accuracies.append(list(map(lambda x: float(x), data["Accuracies"][1:])))
                #     likelihoods.append(list(map(lambda x: float(x), data["Likelihoods"][1:])))
                # else:
                #     accuracies.append([])
                #     likelihoods.append([])

    plt.scatter(chosen_like, chosen_acc, color='black', s=1, label="Chosen solutions")
    # plt.scatter(true_like, [1.0 for x in range(0, len(true_like))], color='turquoise', alpha=0.2, s=0.5, label="True solutions")
    # plt.scatter(lowest_like, lowest_like_acc, color='salmon', alpha=0.2, s=0.5, label="Lowest-likelihood solutions")
    plt.scatter(true_like, [1.0 for x in range(0, len(true_like))], color='turquoise', alpha=0.5, s=1, label="True solutions")
    plt.scatter(lowest_like, lowest_like_acc, color='salmon', alpha=0.5, s=1, label="Lowest-likelihood solutions")
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
    # accuracy_likelihood("PGS000639")
    # accuracy_likelihood("PGS000073")
    # accuracy_likelihood("PGS000037")
    # accuracy_likelihood("PGS002302")
    # true_position_cdf(["PGS000037", "PGS000639", "PGS000073", "PGS002302"])
    # accuracy_cdf(["PGS000037", "PGS000639", "PGS000073", "PGS002302"])
    score_deviation()

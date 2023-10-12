import csv
import matplotlib.pyplot as plt
import numpy as np


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


def solution_ranking(filepath):
    accuracies, positions = [], []
    likelihoods = []
    with open(filepath, 'r') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            accuracies.append(1-float(row[1]))
            positions.append(float(row[2])/float(row[3]))
            likelihoods.append(float(row[6])/float(row[5]))

    plt.hist(accuracies, bins=50, histtype='step', label="Inverse accuracy")
    plt.hist(positions, bins=50, histtype='step', label="Relative position")
    plt.hist(likelihoods, bins=50, histtype='step', label="Likelihood ratio")
    plt.xlabel("Ratio")
    plt.ylabel("Count")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # pairwise("data/prior/PGS000040.pairwise")
    # score_distribution("PGS000040.scores")
    solution_ranking("PGS000639.csv")
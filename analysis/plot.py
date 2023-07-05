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


if __name__ == "__main__":
    pairwise("data/prior/PGS000040.pairwise")
import csv
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import norm

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
    fig.savefig('coverage.pdf', dpi=300, bbox_inches='tight')
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
    fig.savefig('allelefreq.pdf', dpi=300, bbox_inches='tight')


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
import csv
import fnmatch
import json
import random
from cProfile import label

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import pandas as pd
import seaborn as sns
import statistics
from scipy.stats import norm
from scipy.stats import rankdata


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


def load_idv_accuracy(filename):
    directory = "results/guessAccuracy/"
    data = []
    filepath = os.path.join(directory, filename)
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
    return pd.DataFrame(data)


def sequential_idv_accuracy():
    df = load_idv_accuracy("individuals.json")
    df_gaf = load_idv_accuracy("individuals_gaf.json")

    # Calculate and print the median number of SNPs per ancestry
    median_snps_per_ancestry = df.groupby('Ancestry')['SNPs'].median()
    print("Median number of SNPs per ancestry:")
    print(median_snps_per_ancestry)
    # Calculate and print the median number of SNPs total
    median_snps_total = df['SNPs'].median()
    print("Median number of SNPs total:", median_snps_total)

    fig, ax = plt.subplots(figsize=(4, 3))
    sns.boxplot(x='Ancestry', y='GuessAccuracy', data=df, palette='Pastel1', hue='Ancestry')

    # Calculate median GuessAccuracy for each ancestry in df_gaf
    median_gaf_acc = df_gaf.groupby('Ancestry')['GuessAccuracy'].median()
    for i, (ancestry, median) in enumerate(median_gaf_acc.items()):
        ax.plot([i - 0.35, i + 0.35], [median, median], linestyle=':', color='red', linewidth=1)

    median_ref_acc = df.groupby('Ancestry')['ReferenceAccuracy'].median()
    for i, (ancestry, median) in enumerate(median_ref_acc.items()):
        ax.plot([i - 0.35, i + 0.35], [median, median], linestyle='--', alpha=0.7, color='black')
        # Add "Baseline" label above the dashed line
        # if i == 2:
        #     ax.text(i, median - 4, 'Baseline', ha='center', va='bottom', fontsize=10, color='black')

    overall_median_accuracy = df['GuessAccuracy'].median()
    ax.axhline(overall_median_accuracy, color='gray', linestyle='-', linewidth=1, label='Overall Median')
    ax.text(-0.7, overall_median_accuracy, f'{overall_median_accuracy:.1f}', ha='right', va='center',
            fontsize=TICK_SIZE, color='black')

    print("Overall median accuracy:", overall_median_accuracy)
    print("Median accuracy per ancestry:", df.groupby('Ancestry')['GuessAccuracy'].median())
    print("Total samples:", df['Ancestry'].value_counts())

    red_line = Line2D([0], [0], color='red', linestyle=':', label='One-AF loss')
    gray_line = Line2D([0], [0], alpha=0.7, color='black', linestyle='--', label='Baseline')
    # Add the custom legend entry to the plot
    ax.legend(handles=[gray_line, red_line], loc='lower center', ncols=2)

    plt.xlabel('Ancestry')
    plt.ylabel(r'Genotype recovery accuracy, %')
    plt.tight_layout()
    plt.ylim(57, 100)
    fig.savefig('recovery.pdf', dpi=300, bbox_inches='tight')
    # plt.show()


def sequential_loci_accuracy():
    directory = "results/guessAccuracy/"
    data = []
    filepath = os.path.join(directory, "loci.json")
    with open(filepath, 'r') as f:
        content = json.load(f)
    for locus, result in content.items():
        for ancestry, eaf in result['EAF'].items():
            data.append({'EAF': float(eaf), 'Ancestry': ancestry, 'PGS': int(result['SmallestPGS']),
                         'Density': float(result['Density']), 'Weight': float(result['EffectWeight']), 'Accuracy':
                             100*float(result['CorrectGuesses'][ancestry])/float(result['TotalGuesses'][ancestry]),
                         'GWAS': float(result['GWASEAF'])})

    df = pd.DataFrame(data)
    df['EAF_bin'] = pd.cut(df['EAF'], bins=35)
    eaf_df = df.groupby(['EAF_bin', 'Ancestry'], observed=False)['Accuracy'].mean().reset_index()
    eaf_df['EAF_mid'] = eaf_df['EAF_bin'].apply(lambda x: x.mid)
    fig1, ax1 = plt.subplots(figsize=(4, 3))
    sns.lineplot(data=eaf_df, x='EAF_mid', y='Accuracy', ax=ax1, color='teal')
    ax1.set_ylim(80, 100)
    plt.xlabel('Effect Allele Frequency')
    plt.ylabel('Guess accuracy, %')
    plt.tight_layout()
    fig1.savefig('accuracy_eaf.pdf', dpi=300, bbox_inches='tight')

    fig2, ax2 = plt.subplots(figsize=(4, 3))
    # df['PGS_bin'] = pd.cut(df['PGS'], bins=20)
    # pgs_df = df.groupby(['PGS_bin', 'Ancestry'], observed=False)['Accuracy'].mean().reset_index()
    # pgs_df['PGS_mid'] = pgs_df['PGS_bin'].apply(lambda x: x.mid)
    sns.lineplot(x='PGS', y='Accuracy', data=df, ax=ax2, color='coral')
    plt.xlabel('Size of the smallest PGS with the locus')
    plt.ylabel('Guess accuracy, %')
    plt.tight_layout()
    fig2.savefig('accuracy_pgs.pdf', dpi=300, bbox_inches='tight')

    fig3, ax3 = plt.subplots(figsize=(4, 3))
    df['Density_bin'] = pd.cut(df['Density'], bins=30)
    density_df = df.groupby(['Density_bin', 'Ancestry'], observed=False)['Accuracy'].mean().reset_index()
    density_df['Density_mid'] = density_df['Density_bin'].apply(lambda x: x.mid)
    sns.lineplot(x='Density_mid', y='Accuracy', data=density_df, ax=ax3, color='#79ff50')
    plt.xlabel('Density of the PGS with the locus')
    plt.ylabel('Guess accuracy, %')
    plt.tight_layout()
    fig3.savefig('accuracy_density.pdf', dpi=300, bbox_inches='tight')

    fig4, ax4 = plt.subplots(figsize=(4, 3))
    gwas_df = df[df['GWAS'] != 0]
    sns.kdeplot(data=df, x='EAF', hue='Ancestry', ax=ax4, common_norm=False)
    sns.kdeplot(data=gwas_df, x='GWAS', ax=ax4, linestyle='--', color='black', common_norm=False, label='GWAS')
    plt.xlabel('EAF')
    plt.ylabel('Density')
    plt.xlim(0, 1)
    ancestries = df['Ancestry'].unique()
    legend_labels = [*ancestries, "GWAS"]
    legend_colors = sns.color_palette()[:len(ancestries)] + ['black']
    handles = [Line2D([0, 1], [0, 1], color=color) for color in legend_colors[:-1]]
    gwas_handle = Line2D([0, 1], [0, 1], color='black', linestyle='--')
    handles.append(gwas_handle)
    ax4.legend(handles, legend_labels, ncol=2, loc="lower left", columnspacing=0.5)
    # ax4.legend(ncol=2, loc="lower left")
    plt.tight_layout()
    fig4.savefig('accuracy_ancestry.pdf', dpi=300, bbox_inches='tight')
    # plt.show()

    fig5, ax5 = plt.subplots(figsize=(4, 3))
    df['Weight_bin'] = pd.cut(df['Weight'], bins=40)
    weight_df = df.groupby(['Weight_bin', 'Ancestry'], observed=False)['Accuracy'].mean().reset_index()
    weight_df['Weight_mid'] = weight_df['Weight_bin'].apply(lambda x: x.mid)
    sns.lineplot(data=weight_df, x='Weight_mid', y='Accuracy', ax=ax5)
    ax5.set_ylim(70, 100)
    plt.xlabel('Effect Weight')
    plt.ylabel('Accuracy')
    plt.tight_layout()
    # fig5.savefig('accuracy_weight.pdf', dpi=300, bbox_inches='tight')
    # plt.show()


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
    filepath = os.path.join(directory, "scores_1000G.json")
    data = []
    with open(filepath, 'r') as f:
        content = json.load(f)
    for result in content:
        data.append({'NumVariants': int(result['NumVariants']),
                     'TotalPresentScores': int(result['TotalPresentScores']),
                        'TotalPossibleScores': int(result['TotalPossibleScores']),
                        'RealPercentageUnique': int(result['RealPercentageUnique']),
                        'PredictedPercentageUnique': int(result['PredictedPercentageUnique']),
                     'AnonSize': np.median([float(x) for x in result['AnonymitySets']])})

    first_color = 'teal'
    second_color = 'coral'
    df = pd.DataFrame(data)
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.lineplot(x='TotalPossibleScores', y='RealPercentageUnique', data=df, ax=ax, color=first_color, linestyle='-',
                 label='Dataset')
    sns.lineplot(x='TotalPossibleScores', y='PredictedPercentageUnique', data=df, ax=ax, color=second_color,
                 linestyle='--', label='Predicted')
    plt.xlabel('Number of possible scores')
    plt.ylabel('Percentage of unique scores')
    ax.set_xscale('log')

    plt.legend()
    plt.tight_layout()
    # fig.savefig('uniqueness.pdf', dpi=300, bbox_inches='tight')
    plt.show()


def score_uniqueness_old():
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

    data = []
    related = read_related_individuals()
    with open(filepath, 'r') as f:
        content = json.load(f)
    for idv, entry in content.items():
        # if idv not in data:
        #     data[idv] = {}
        # ancestries[idv] = entry["Ancestry"]
        for target, relation in entry["Relations"].items():
            if target == idv:
                data.append({'Category': 'Self', 'Ancestry': entry["Ancestry"], 'Method': 'King',
                            'Kinship': float(relation["King"][0])/float(relation["King"][1])})
                data.append({'Category': 'Self', 'Ancestry': entry["Ancestry"], 'Method': 'Count',
                            'Kinship': float(relation["Count"][0])/float(relation["Count"][1])})
                continue
            if idv in related:
                for relative in related[idv]:
                    if target == relative[0]:
                        category = ''
                        if relative[1] == '1':
                            category = '1st degree'
                        elif relative[1] == '2':
                            category = '2nd degree'
                        data.append({'Category': category, 'Ancestry': entry["Ancestry"], 'Method': 'King',
                                     'Kinship': float(relation["King"][0])/float(relation["King"][1])})
                        data.append({'Category': category, 'Ancestry': entry["Ancestry"], 'Method': 'Count',
                                    'Kinship': float(relation["Count"][0])/float(relation["Count"][1])})
                continue
            data.append({'Category': 'Everyone else', 'Ancestry': entry["Ancestry"], 'Method': 'King',
                        'Kinship': float(relation["King"][0])/float(relation["King"][1])})
            data.append({'Category': 'Everyone else', 'Ancestry': entry["Ancestry"], 'Method': 'Count',
                        'Kinship': float(relation["Count"][0])/float(relation["Count"][1])})

    df = pd.DataFrame(data)
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


    # divided = {}
    # for idv in data:
    #     divided[idv] = {"Count": [], "King": []}
    #     for target in data[idv]:
    #         divided[idv]["Count"].append((target, data[idv][target]["Count"][0] / data[idv][target]["Count"][1]))
    #         divided[idv]["King"].append((target, data[idv][target]["King"][0] / data[idv][target]["King"][1]))
    #     # sort by the match count and king score
    #     divided[idv]["Count"] = sorted(divided[idv]["Count"], key=lambda x: x[1], reverse=True)
    #     divided[idv]["King"] = sorted(divided[idv]["King"], key=lambda x: x[1], reverse=True)
    #     # print(f"IDV: {idv}, Count: {divided[idv]['Count'][0]}, King: {divided[idv]['King'][0]}")
    #
    # parsed = {"SelfCountPos": [], "SelfKingPos": [],
    #           "RelativeCountPos": [], "RelativeKingPos": [],
    #           "ReferenceCountPos": [], "ReferenceKingPos": [],
    #           "SelfCountAcc": [], "SelfKingAcc": [],
    #           "RelativeCountAcc": [], "RelativeKingAcc": [],
    #           "ReferenceCountAcc": [], "ReferenceKingAcc": [],
    #           "EveryCountRatio": [], "EveryKingRatio": []}
    #
    # for idv in divided:
    #     relative_found = False
    #     for i, tpl in enumerate(divided[idv]["Count"]):
    #         if tpl[0] == idv:
    #             parsed["SelfCountPos"].append(i)
    #             parsed["SelfCountAcc"].append(tpl[1])
    #         if tpl[0] == "reference":
    #             parsed["ReferenceCountPos"].append(i)
    #             parsed["ReferenceCountAcc"].append(tpl[1])
    #         if not relative_found and idv in related and tpl[0] in related[idv]:
    #             parsed["RelativeCountPos"].append(i)
    #             parsed["RelativeCountAcc"].append(tpl[1])
    #             relative_found = True
    #         parsed["EveryCountRatio"].append(tpl[1])
    # for idv in divided:
    #     relative_found = False
    #     for i, tpl in enumerate(divided[idv]["King"]):
    #         if tpl[0] == idv:
    #             parsed["SelfKingPos"].append(i)
    #             parsed["SelfKingAcc"].append(tpl[1])
    #         if tpl[0] == "reference":
    #             parsed["ReferenceKingPos"].append(i)
    #             parsed["ReferenceKingAcc"].append(tpl[1])
    #         if not relative_found and idv in related and tpl[0] in related[idv]:
    #             parsed["RelativeKingPos"].append(i)
    #             parsed["RelativeKingAcc"].append(tpl[1])
    #             relative_found = True
    #         parsed["EveryKingRatio"].append(tpl[1])
    #
    # print(f"SelfCountAcc: {np.median(parsed['SelfCountAcc'])}, {np.mean(parsed['SelfCountAcc'])}")
    # print(f"RelativeCountAcc: {np.median(parsed['RelativeCountAcc'])}, {np.mean(parsed['RelativeCountAcc'])}")
    # print(f"ReferenceCountAcc: {np.median(parsed['ReferenceCountAcc'])}, {np.mean(parsed['ReferenceCountAcc'])}")
    # print(f"EveryCountRatio: {np.median(parsed['EveryCountRatio'])}, {np.mean(parsed['EveryCountRatio'])}")
    # print(f"SelfKingAcc: {np.median(parsed['SelfKingAcc'])}, {np.mean(parsed['SelfKingAcc'])}")
    # print(f"RelativeKingAcc: {np.median(parsed['RelativeKingAcc'])}, {np.mean(parsed['RelativeKingAcc'])}")
    # print(f"ReferenceKingAcc: {np.median(parsed['ReferenceKingAcc'])}, {np.mean(parsed['ReferenceKingAcc'])}")
    # print(f"EveryKingRatio: {np.median(parsed['EveryKingRatio'])}, {np.mean(parsed['EveryKingRatio'])}")
    #
    # palette = {}
    # for key in parsed:
    #     if key.startswith("Self"):
    #         palette[key] = sns.color_palette("pastel")[0]
    #     elif key.startswith("Relative"):
    #         palette[key] = sns.color_palette("pastel")[1]
    #     elif key.startswith("Reference"):
    #         palette[key] = sns.color_palette("pastel")[2]
    #     else:
    #         palette[key] = sns.color_palette("pastel")[3]
    #
    # keys = [["SelfCountPos", "RelativeCountPos", "ReferenceCountPos"],
    #         ["SelfKingPos", "RelativeKingPos", "ReferenceKingPos"],
    #         ["SelfCountAcc", "RelativeCountAcc", "ReferenceCountAcc", "EveryCountRatio"],
    #         ["SelfKingAcc", "RelativeKingAcc", "ReferenceKingAcc", "EveryKingRatio"]]
    # df = []
    # for key in keys:
    #     df.append(prepare_data(parsed, key))
    # fig1, axes1 = plt.subplots(1, 2, figsize=(12, 5))
    # for i, ax1 in enumerate(axes1):
    #     ax1.invert_yaxis()
    #     ax1.set_ylabel('Rank')
    #     ax1.set_xlabel('')
    #     ax1.set_xticklabels(["Self", "Relative", "Reference"])
    #     sns.boxplot(x='Versus', y='Value', data=df[i], hue='Versus', palette=palette, ax=axes1[i])
    #     medians = df[i].groupby(['Versus'])['Value'].median()
    #     print(f"Medians for df[{i}]:", medians)  # Debug print
    #     # for tick, label in zip(range(len(medians)), ax1.get_xticklabels()):
    #     #     ax1.annotate(f'{medians[label.get_text()]:.2f}',
    #     #                  xy=(tick, medians[label.get_text()]),
    #     #                  xytext=(0, 5),  # 5 points vertical offset
    #     #                  textcoords='offset points',
    #     #                  ha='center', va='center',
    #     #                  color='red', fontsize=10, fontweight='bold')
    # axes1[0].set_title('Match Count Based')
    # axes1[0].set_ylim(2550, -50)
    # axes1[1].set_title('KING Based')
    # axes1[1].set_ylim(2550, -50)
    # # fig1.savefig('unimputed-pos.png', dpi=300, bbox_inches='tight')
    # # fig1.savefig('imputed-pos.png', dpi=300, bbox_inches='tight')
    #
    # fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
    # for i, ax2 in enumerate(axes2):
    #     ax2.set_xlabel('')
    #     ax2.set_xticklabels(["Self", "Relative", "Reference", "All"])
    #     sns.boxplot(x='Versus', y='Value', data=df[i + int(len(df)/2)], hue='Versus', palette=palette, ax=axes2[i])
    #     medians = df[i + int(len(df)/2)].groupby(['Versus'])['Value'].median()
    #     print(f"Medians for df[{i}]:", medians)  # Debug print
    #     # for tick, label in zip(range(len(medians)), ax2.get_xticklabels()):
    #     #     ax2.annotate(f'{medians[label.get_text()]:.2f}',
    #     #                  xy=(tick, medians[label.get_text()]),
    #     #                  xytext=(0, 5),  # 5 points vertical offset
    #     #                  textcoords='offset points',
    #     #                  ha='center', va='center',
    #     #                  color='red', fontsize=10, fontweight='bold')
    # axes2[0].set_title('Match Count Based')
    # axes2[0].set_ylabel('Match ratio')
    # axes2[1].set_ylabel('KING score')
    # axes2[1].set_title('KING Based')
    # # fig2.savefig('unimputed-acc.png', dpi=300, bbox_inches='tight')
    # # fig2.savefig('imputed-acc.png', dpi=300, bbox_inches='tight')

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
        indf = pd.read_csv(data_file, sep='\t', index_col=0)

    # Prune the DataFrame to include only rows starting with '$' and columns not starting with '$'
    pruned_df = indf.loc[indf.index.str.startswith('$'), ~indf.columns.str.startswith('$')]
    truth_df = indf.loc[~indf.index.str.startswith('$'), ~indf.columns.str.startswith('$')]
    related = read_related_individuals()

    results = []
    for idv in pruned_df.index:
        target = idv[1:]  # Remove the '$' symbol
        if target not in pruned_df.columns:
            continue
        kinship = pruned_df.at[idv, target]
        self_rank = pruned_df.loc[idv].rank(ascending=False)[target]
        results.append({'Category': 'Self', 'Kinship': kinship, 'Rank': self_rank, 'Method': method})
        results.append({'Category': 'Self', 'Kinship': truth_df.at[target, target],
                        'Rank': truth_df.loc[target].rank(ascending=False)[target], 'Method': 'True'+method})
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
                results.append({'Category': category, 'Kinship': truth_df.at[target, relative[0]],
                                'Rank': truth_df.loc[target].rank(ascending=False)[relative[0]], 'Method': 'True'+method})
                relatives.append(relative[0])
        other_columns = pruned_df.columns[~pruned_df.columns.isin([target] + relatives)]
        kinships = pruned_df.loc[idv, other_columns]
        ranks = pruned_df.loc[idv].rank(ascending=False)[other_columns]
        temp_df = pd.DataFrame({
            'Category': 'Everyone else',
            'Kinship': kinships.values,
            'Rank': ranks.values,
            'Method': method
        })
        results.extend(temp_df.to_dict('records'))
        #
        other_columns = truth_df.columns[~truth_df.columns.isin([target] + relatives)]
        kinships = truth_df.loc[target, other_columns]
        ranks = truth_df.loc[target].rank(ascending=False)[other_columns]
        temp_df = pd.DataFrame({
            'Category': 'Everyone else',
            'Kinship': kinships.values,
            'Rank': ranks.values,
            'Method': 'True'+method
        })
    results.extend(temp_df.to_dict('records'))

    return pd.DataFrame(results)


def plot_rank_cdf(df, category, ax):
    colors = sns.color_palette("deep")
    # Filter the dataframe for the given category
    category_df = df[df['Category'] == category]
    methods = ['KING', 'GCTA']
    for i, method in enumerate(category_df['Method'].unique()):
        if method not in methods:
            continue
        method_df = category_df[category_df['Method'] == method]
        # Sort the ranks for the cumulative distribution
        sorted_ranks = np.sort(method_df['Rank'])
        # Calculate the cumulative distribution
        cdf = np.arange(1, len(sorted_ranks) + 1) / len(sorted_ranks) * 100
        ax.plot(sorted_ranks, cdf, label=f'{method}', color=colors[i], linewidth=2)
    ax.set_title(f'{category}')
    ax.set_xlim(-5, 100)


def deanonymization_accuracy():
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

    king_df = load_results("ibd/nophase/plink.king.id", "ibd/nophase/plink.king", 'KING')
    rel_df = load_results("ibd/nophase/plink.rel.id", "ibd/nophase/plink.rel", 'GCTA')
    rel_df['Kinship'] = rel_df['Kinship'] / 2 # Relatedness -> Kinship
    rel_df['Kinship'] = [0.5 if k > 0.5 else k for k in rel_df['Kinship']]
    # refined_df = load_results(None, "ibd/refined.ibd", 'Refined IBD')
    # df = pd.concat([king_df, rel_df, refined_df])
    df = pd.concat([king_df, rel_df])

    palette = sns.color_palette("pastel")
    for category in df['Category'].unique():
        print(f'==== {category} ====')
        for method in df['Method'].unique():
            category_df = df[(df['Category'] == category) & (df['Method'] == method)]
            print(f'{method} kinship: {category_df["Kinship"].median()}, '
                  f'rank: {category_df["Rank"].median()}')
    self_cutoff = 0.3535
    first_degree_cutoff = 0.1768
    second_degree_cutoff = 0.089
    fig1, ax1 = plt.subplots(figsize=(7, 3))
    guess_methods = ['KING', 'GCTA']
    ground_truth_methods = ['TrueKING', 'TrueGCTA']
    method_df = df[df['Method'].isin(guess_methods)]
    sns.boxplot(x='Category', y='Kinship', data=method_df, hue='Method', palette='pastel', ax=ax1,
                order=['Self', '1st degree', '2nd degree', 'Everyone else'])

    box_width = 0.8 / len(ground_truth_methods)
    for method in ground_truth_methods:
        for category in ['Self', '1st degree', '2nd degree', 'Everyone else']:
            category_data = df[(df['Method'] == method) & (df['Category'] == category)]
            median_value = category_data['Kinship'].median()
            # Find the position of the category on the X axis
            category_pos = ['Self', '1st degree', '2nd degree', 'Everyone else'].index(category)
            method_position_offset = ground_truth_methods.index(method) * box_width
            x_position = category_pos - 0.4 + method_position_offset + box_width / 2
            print(f'{method} median: {median_value}, position: {category_pos}')
            ax1.plot([x_position - box_width / 2, x_position + box_width / 2], [median_value, median_value],
                     color=palette[3], linestyle='-')

    # Adding the legend for the red lines
    ax1.plot([], [], color=palette[3], linestyle='-', label='Max accuracy')

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
    fig1.savefig('kinship_score.pdf', dpi=300, bbox_inches='tight')

    fig2, ax2 = plt.subplots(figsize=(5, 3))
    categories = ['Self', '1st degree', '2nd degree']
    sns.boxplot(x='Category', y='Rank', data=method_df[method_df['Category'].isin(categories)], hue='Method', palette='pastel',
                ax=ax2, order=categories)
    ax2.invert_yaxis()
    ax2.set_ylabel('Rank (/2535)')
    ax2.set_xlabel('')
    ax2.legend(title="")
    fig2.savefig('kinship_rank.pdf', dpi=300, bbox_inches='tight')

    fig3, axes3 = plt.subplots(1, 3, figsize=(6, 3), sharey=True)
    # Plot CDFs for each category
    categories = ['Self', '1st degree', '2nd degree']
    for i, category in enumerate(categories):
        plot_rank_cdf(method_df, category, axes3[i])
    axes3[0].set_ylabel('Percentage')
    axes3[1].set_xlabel('Rank (/2535)')
    handles, labels = axes3[0].get_legend_handles_labels()
    fig3.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.05), title='', ncol=len(labels), frameon=False)
    # plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.tight_layout()
    fig3.savefig('kinship_cdf.pdf', dpi=300, bbox_inches='tight')
    # plt.show()

    results = pd.DataFrame()
    pd.set_option('display.max_columns', None)  # Show all columns
    total_count = df.groupby(['Category', 'Method']).size().reset_index(name='Total_Count')
    for rank_value in range(1, 9):
        rank_df = df[df['Rank'] <= rank_value].groupby(['Category', 'Method']).size().reset_index(name=f'Rank_<={rank_value}_Count')
        merged_df = pd.merge(rank_df, total_count, on=['Category', 'Method'], how='right')
        merged_df[f'Rank_<={rank_value}_Count'].fillna(0, inplace=True)
        merged_df[f'Percentage_Rank_<={rank_value}'] = (merged_df[f'Rank_<={rank_value}_Count'] / merged_df['Total_Count']) * 100
        if results.empty:
            results = merged_df[['Category', 'Method', f'Percentage_Rank_<={rank_value}']]
        else:
            results = pd.merge(results, merged_df[['Category', 'Method', f'Percentage_Rank_<={rank_value}']],
                               on=['Category', 'Method'], how='outer')
    print(results)

    fig4, ax4 = plt.subplots(figsize=(4, 3))
    df_king = df[df['Method'] == 'KING']
    categories = ['Self', '1st degree', '2nd degree', 'Everyone else']
    colors = {'Self': 'black', '1st degree': 'coral', '2nd degree': 'green', 'Everyone else': 'purple'}
    patterns = {'Self': '-', '1st degree': 'dotted', '2nd degree': '--', 'Everyone else': '-.'}
    for category in categories:
        sns.kdeplot(df_king[df_king['Category'] == category]['Kinship'], label=category, common_norm=False,
                    ax=ax4, color=colors[category], linestyle=patterns[category])
    ax4.set_xlabel('Kinship coefficient')
    ax4.set_ylabel('Density')
    ax4.set_xlim(-0.4, 0.5)
    ax4.legend(title='', loc="upper left")
    plt.tight_layout()
    fig4.savefig('kinship_original.pdf', dpi=300, bbox_inches='tight')
    # plt.show()


def prs_prediction():
    filepath = "results/predict/prediction.json"
    data = []
    with open(filepath, 'r') as f:
        content = json.load(f)
    for prs_id, results in content.items():
        for i, prd in enumerate(results['Predicted']):
            data.append({'PRS': prs_id, 'Known': float(results["Known"]), 'Total': float(results["Total"]),
                         'Max': float(results["Max"]), 'Min': float(results["Min"]),
                         'Mean': float(results["Means"][results["Ancestries"][i]]),
                         'Std': float(results['Stds'][results["Ancestries"][i]]),
                         'Ancestry': results["Ancestries"][i],
                         'Predicted': float(results['Predicted'][i]), 'Real': float(results['Real'][i])})

    df = pd.DataFrame(data)
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    bin_size = 0.05
    df['Ratio'] = (df['Known'] / df['Total']).round(3)
    df['Binned_Ratio'] = (bin_size * np.round(df['Ratio'] / bin_size)).round(3)
    df['Threshold'] = df['Mean'] + df['Std']
    df['Above_Threshold'] = df['Real'] > df['Threshold']
    df['Predicted_Above_Threshold'] = df['Predicted'] > df['Threshold']

    df['Real_Percentile'] = norm.cdf(df['Real'], loc=df['Mean'], scale=df['Std'])
    df['Predicted_Percentile'] = norm.cdf(df['Predicted'], loc=df['Mean'], scale=df['Std'])
    df['Percentile_Diff'] = abs(df['Real_Percentile'] - df['Predicted_Percentile']) * 100

    def calculate_metrics(group):
        tp = np.sum((group['Above_Threshold'] == True) & (group['Predicted_Above_Threshold'] == True))
        fp = np.sum((group['Above_Threshold'] == False) & (group['Predicted_Above_Threshold'] == True))
        fn = np.sum((group['Above_Threshold'] == True) & (group['Predicted_Above_Threshold'] == False))

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        std_dev = group['Predicted'].std()
        percentile_diff = group['Percentile_Diff'].median()
        percentile_diff_std = group['Percentile_Diff'].std()
        high_risk_percentile_diff = group[group['Real'] > group['Threshold']]['Percentile_Diff'].median()
        high_risk_percentile_diff_std = group[group['Real'] > group['Threshold']]['Percentile_Diff'].std()

        return pd.Series({'Precision': precision, 'Recall': recall, 'Std_Dev': std_dev,
                          'Mean_Percentile_Diff': percentile_diff, 'High_Risk_Percentile_Diff': high_risk_percentile_diff,
                          'Percentile_std': percentile_diff_std, 'High_Risk_Percentile_Std': high_risk_percentile_diff_std})

    grouped = df.groupby('Binned_Ratio').apply(calculate_metrics).reset_index()

    plt.figure(figsize=(4, 3))
    sns.lineplot(data=grouped, x='Binned_Ratio', y='Precision', label='Precision', color='coral')
    sns.lineplot(data=grouped, x='Binned_Ratio', y='Recall', label='Recall', color='teal')
    plt.fill_between(grouped['Binned_Ratio'], grouped['Precision'] - grouped['Std_Dev'],
                     grouped['Precision'] + grouped['Std_Dev'], color='teal', alpha=0.1)
    plt.fill_between(grouped['Binned_Ratio'], grouped['Recall'] - grouped['Std_Dev'],
                     grouped['Recall'] + grouped['Std_Dev'], color='coral', alpha=0.1)
    plt.xlabel('Guessed SNPs ratio')
    plt.ylabel('')
    plt.gca().invert_xaxis()
    plt.ylim(0, 1)
    plt.legend(loc='lower left')
    # plt.savefig('prediction_highrisk_by_population.pdf', dpi=300, bbox_inches='tight')
    # plt.show()

    plt.figure(figsize=(4, 3))
    sns.lineplot(data=grouped, x='Binned_Ratio', y='Mean_Percentile_Diff', label='All individuals', linestyle='-',
                 color='teal')
    # plt.fill_between(grouped['Binned_Ratio'], grouped['Mean_Percentile_Diff'] - grouped['Percentile_std'],
    #                  grouped['Mean_Percentile_Diff'] + grouped['Percentile_std'], color='teal', alpha=0.1)
    sns.lineplot(data=grouped, x='Binned_Ratio', y='High_Risk_Percentile_Diff', label=r'High risk ($>\mu$+$\sigma$)',
                 linestyle='--', color='coral')
    # plt.fill_between(grouped['Binned_Ratio'], grouped['High_Risk_Percentile_Diff'] - grouped['High_Risk_Percentile_Std'],
    #                  grouped['High_Risk_Percentile_Diff'] + grouped['High_Risk_Percentile_Std'], color='coral', alpha=0.1)
    plt.xlabel('Guessed SNPs ratio')
    plt.ylabel('Median percentile error')
    plt.legend(loc='upper left')
    plt.ylim(0, 50)
    plt.xticks(np.arange(0.1, 0.9, 0.1))
    plt.gca().invert_xaxis()
    # plt.show()
    plt.savefig('prediction_median_percentile.pdf', dpi=300, bbox_inches='tight')


LABEL_SIZE = 13
TICK_SIZE = 13
LEGEND_SIZE = 11
def prepare_for_latex():
    params = {
        'axes.labelsize': LABEL_SIZE,
        'xtick.labelsize': TICK_SIZE,
        'ytick.labelsize': TICK_SIZE,
        'legend.fontsize': LEGEND_SIZE,
        'text.usetex': True,
        'font.family': 'serif',
        'font.serif': 'Computer Modern',
        'savefig.bbox': 'tight',
        'savefig.format': 'pdf',
    }
    mpl.rcParams.update(params)


if __name__ == "__main__":
    prepare_for_latex()
    # pairwise("data/prior/PGS000040.pairwise")
    # score_distribution("PGS000040.scores")
    # solution_ranking("PGS000639")
    # solution_ranking("PGS002302")
    # solution_ranking("PGS000073")
    # solution_ranking("PGS000037")
    # true_position_cdf(["PGS000037", "PGS000639", "PGS000073", "PGS002302"])
    # accuracy_cdf(["PGS000037", "PGS000639", "PGS000073", "PGS002302"])
    # score_deviation()
    # loci_coverage()
    # king_test()
    # kinship_experiment()
    # random_hist()
    # guessed_mia()
    # full_imputed()
    # linking_accuracy()
    # imputation_accuracy()
    # percentile_prediction()
    # random_bars()
    # score_uniqueness_old()
    score_uniqueness()
    # sequential()
    # sequential_idv_accuracy()
    # sequential_loci_accuracy()
    # af_hist()
    # plot_normal_distribution_with_fill()
    # king_accuracy()
    # ibd_accuracy()
    # plink_accuracy()
    # deanonymization_accuracy()
    # prs_prediction()
import csv
import fnmatch
import json
import random
from ast import Index
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
    fig.savefig('coverage.png', dpi=300, bbox_inches='tight')
    # plt.show()


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

    # # Calculate median GuessAccuracy for each ancestry in df_gaf
    # median_gaf_acc = df_gaf.groupby('Ancestry')['GuessAccuracy'].median()
    # print("Global AF:", median_gaf_acc)
    # print("Global AF median:", df_gaf['GuessAccuracy'].median())
    # for i, (ancestry, median) in enumerate(median_gaf_acc.items()):
    #     ax.plot([i - 0.35, i + 0.35], [median, median], linestyle=':', color='red', linewidth=1)

    median_ref_acc = df.groupby('Ancestry')['ReferenceAccuracy'].median()
    for i, (ancestry, median) in enumerate(median_ref_acc.items()):
        ax.plot([i - 0.35, i + 0.35], [median, median], linestyle='--', alpha=0.7, color='black')
        # Add "Baseline" label above the dashed line
        if i == 2:
            ax.text(i, median - 4, 'Baseline', ha='center', va='bottom', fontsize=LABEL_SIZE, color='black')

    overall_median_accuracy = df['GuessAccuracy'].median()
    ax.axhline(overall_median_accuracy, color='gray', linestyle='-', linewidth=1, label='Overall Median')
    ax.text(-0.7, overall_median_accuracy, f'{overall_median_accuracy:.1f}', ha='right', va='center',
            fontsize=TICK_SIZE, color='black')

    print("Overall median accuracy:", overall_median_accuracy)
    print("Median accuracy per ancestry:", df.groupby('Ancestry')['GuessAccuracy'].median())
    print("Total samples:", df['Ancestry'].value_counts())

    # red_line = Line2D([0], [0], color='red', linestyle=':', label='One-AF loss')
    # gray_line = Line2D([0], [0], alpha=0.7, color='black', linestyle='--', label='Baseline')
    ## Add the custom legend entry to the plot
    # ax.legend(handles=[gray_line, red_line], loc='lower center', ncols=2)

    plt.xlabel('Ancestry')
    # plt.ylabel(r'Genotype recovery accuracy, \%')
    plt.ylabel(r'Recovery accuracy, \%')
    plt.tight_layout()
    plt.ylim(64, 100)
    # fig.savefig('recovery.pdf', dpi=300, bbox_inches='tight')
    plt.show()


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
    ax1.set_ylim(64, 100)
    ax1.set_xlim(0, 1)
    plt.xticks([0, 0.25, 0.5, 0.75, 1], ['0', '0.25', '0.5', '0.75', '1'])
    plt.xlabel('Effect Allele Frequency')
    plt.ylabel(r'Recovery accuracy, \%')
    plt.tight_layout()
    fig1.savefig('accuracy_eaf.pdf', dpi=300, bbox_inches='tight')

    fig2, ax2 = plt.subplots(figsize=(4, 3))
    # df['PGS_bin'] = pd.cut(df['PGS'], bins=20)
    # pgs_df = df.groupby(['PGS_bin', 'Ancestry'], observed=False)['Accuracy'].mean().reset_index()
    # pgs_df['PGS_mid'] = pgs_df['PGS_bin'].apply(lambda x: x.mid)
    sns.lineplot(x='PGS', y='Accuracy', data=df, ax=ax2, color='coral')
    plt.xlabel('Size of the smallest PGS with the locus')
    plt.ylabel(r'Recovery accuracy, \%')
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
    plt.xlabel('Effect Allele Frequency')
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
    plt.xticks([0, 0.25, 0.5, 0.75, 1], ['0', '0.25', '0.5', '0.75', '1'])
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


def load_uniqueness_data(dataset):
    directory = "results/uniqueness/"
    filepath = os.path.join(directory, f"scores_{dataset}.json")
    data = []
    with open(filepath, 'r') as f:
        content = json.load(f)
    for result in content:
        data.append({'Dataset': dataset, 'NumVariants': int(result['NumVariants']),
                     'TotalPresentScores': int(result['TotalPresentScores']),
                     'TotalPossibleScores': int(result['TotalPossibleScores']),
                     'RealPercentageUnique': int(result['RealPercentageUnique']),
                     'PredictedPercentageUnique': int(result['PredictedPercentageUnique']),
                     'RealAnonSize': int(result['RealMedianAnonymitySet']),
                     'PredictedAnonSize': int(result['PredictedMedianAnonymitySet']) if
                     int(result['PredictedMedianAnonymitySet']) > 1 or int(result['PredictedMedianAnonymitySet']) == -1
                     else 1})
    return pd.DataFrame(data)


def score_uniqueness():
    ggdf = load_uniqueness_data('1000Genomes')
    ukbdf = load_uniqueness_data('UKBiobank')
    df = pd.concat([ggdf, ukbdf])

    # Create bins for powers of 10
    # bins = [0] + [10**i for i in range(1, int(np.log10(df['TotalPossibleScores'].max())) + 3)]
    bins = []
    for i in range(1, int(np.log10(df['TotalPossibleScores'].max())) + 2):
        bins.append(10**(i-1))
        bins.append(10**(i-1) * 10**0.5)
    bins.append(10**i)
    df['TotalPossibleScoresBin'] = pd.cut(df['TotalPossibleScores'], bins=bins, right=False)
    # Group by bins and calculate mean and std for real values
    grouped_real = df.groupby(['TotalPossibleScoresBin', 'Dataset'], observed=True).agg(
        # RealPercentageUnique_mean=('RealPercentageUnique', 'mean'),
        RealPercentageUnique_median=('RealPercentageUnique', 'median'),
        RealPercentageUnique_std=('RealPercentageUnique', 'std'),
        # RealAnonSize_mean=('RealAnonSize', 'mean'),
        RealAnonSize_median=('RealAnonSize', 'median'),
        RealAnonSize_std=('RealAnonSize', 'std'),
        TotalPossibleScores_mean=('TotalPossibleScores', 'mean'),
        NumVariants_median=('NumVariants', 'median')
    ).reset_index()

    df['PredictedPercentageUnique_filtered'] = df['PredictedPercentageUnique'].where(df['PredictedPercentageUnique'] > 0)
    grouped_predicted = df.groupby([
        'TotalPossibleScoresBin', 'Dataset'], observed=True).agg(
        # PredictedPercentageUnique_mean=('PredictedPercentageUnique_filtered', 'mean'),
        PredictedPercentageUnique_median=('PredictedPercentageUnique_filtered', 'median'),
        PredictedPercentageUnique_std=('PredictedPercentageUnique_filtered', 'std'),
        # PredictedAnonSize_mean=('PredictedAnonSize', 'mean'),
        PredictedAnonSize_median=('PredictedAnonSize', 'median'),
        PredictedAnonSize_std=('PredictedAnonSize', 'std'),
        NumVariants_median=('NumVariants', 'median')
    ).reset_index()

    num_unique = {}
    anon_group = {}
    for dataset in df['Dataset'].unique():
        dataset_group = grouped_real[grouped_real['Dataset'] == dataset]
        try:
            unique_95_bin = dataset_group[dataset_group['RealPercentageUnique_median'] >= 95].iloc[0]
            num_unique[dataset] = unique_95_bin['TotalPossibleScores_mean']
            print(f"{dataset} - Bin: {unique_95_bin['TotalPossibleScoresBin']}, Median NumVariants: {unique_95_bin['NumVariants_median']} (RealPercentageUnique >= 95)")
        except IndexError:
            pass
        try:
            anon_2_bin = dataset_group[dataset_group['RealAnonSize_median'] <= 2].iloc[0]
            anon_group[dataset] = anon_2_bin['TotalPossibleScores_mean']
            print(f"{dataset} - Bin: {anon_2_bin['TotalPossibleScoresBin']}, Median NumVariants: {anon_2_bin['NumVariants_median']} (RealAnonSize <= 2)")
        except IndexError:
            pass

        # Print the results for each dataset
        print(f"{dataset} - Mean TotalPossibleScores at 95% unique: {num_unique[dataset]}")
        print(f"{dataset} - Mean TotalPossibleScores when RealAnonSize reaches 2: {anon_group[dataset]}")

    first_color, second_color = 'teal', 'coral'
    first_marker, second_marker = '.', '*'
    observed, estimated = 'observed', 'estimated'

    fig1, ax1 = plt.subplots(figsize=(5, 4))
    for dataset, color, label_prefix, marker in [('1000Genomes', first_color, '1000G', first_marker),
                                                 ('UKBiobank', second_color, 'UKBB', second_marker)]:
        # Plot real values for percentage uniqueness
        subset_real = grouped_real[grouped_real['Dataset'] == dataset]
        x_real = subset_real['TotalPossibleScoresBin'].apply(lambda x: x.mid)
        y_real = subset_real['RealPercentageUnique_median']
        y_real_std = np.clip(subset_real['RealPercentageUnique_std'], 0, 100 - y_real)

        ax1.plot(x_real, y_real, color=color, label=f'{label_prefix} {observed}', marker=marker)
        ax1.fill_between(x_real, y_real - y_real_std, y_real + y_real_std, color=color, alpha=0.1)

        # Plot predicted values
        subset_pred = grouped_predicted[grouped_predicted['Dataset'] == dataset]
        if not subset_pred.empty:
            x_pred = subset_pred['TotalPossibleScoresBin'].apply(lambda x: x.mid)
            y_pred = subset_pred['PredictedPercentageUnique_median']
            ax1.plot(x_pred, y_pred, color=color, linestyle=':', label=f'{label_prefix} {estimated}', marker=marker)

    ax1.text(x=0.2, y=93, s='95', va='center', ha='right', color='gray', fontsize=TICK_SIZE)
    ax1.axhline(y=95, color='lightgray', linestyle='--', linewidth=1.5)
    # for dataset in df['Dataset'].unique():
    #     ax1.vlines(x=num_unique[dataset], ymin=0.8, ymax=95, color='lightgray', linestyle=':', linewidth=1.5)
    ax1.set_ylim(0, 105)
    ax1.set_xscale('log')
    ax1.set_xlabel('Number of possible scores')
    ax1.set_ylabel(r'Identifiable individuals, \%')
    ax1.legend()
    plt.tight_layout()

    # Second graph for anonymity set size
    fig2, ax2 = plt.subplots(figsize=(5, 4))
    for dataset, color, label_prefix, marker in [('1000Genomes', first_color, '1000G', first_marker),
                                                 ('UKBiobank', second_color, 'UKBB', second_marker)]:
        subset_anon_real = grouped_real[grouped_real['Dataset'] == dataset]
        x_anon_real = subset_anon_real['TotalPossibleScoresBin'].apply(lambda x: x.mid)
        y_real_anon = subset_anon_real['RealAnonSize_median']
        y_real_anon_std = np.clip(subset_anon_real['RealAnonSize_std'], 0, y_real_anon - 1)

        ax2.plot(x_anon_real, y_real_anon, color=color, label=f'{label_prefix} {observed}', marker=marker)
        ax2.fill_between(x_anon_real, np.maximum(y_real_anon - y_real_anon_std, 1), y_real_anon + y_real_anon_std,
                         color=color, alpha=0.1)

        # Plot predicted anonymity set size
        subset_anon_pred = grouped_predicted[grouped_predicted['Dataset'] == dataset]
        if not subset_anon_pred.empty:
            x_anon_pred = subset_anon_pred['TotalPossibleScoresBin'].apply(lambda x: x.mid)
            y_pred_anon = subset_anon_pred['PredictedAnonSize_median']
            ax2.plot(x_anon_pred, y_pred_anon, color=color, linestyle=':', label=f'{label_prefix} {estimated}',
                     marker=marker)

    ax2.text(x=0.2, y=2, s='2', va='center', ha='right', color='gray', fontsize=TICK_SIZE)
    ax2.axhline(y=2, color='lightgray', linestyle='--', linewidth=1.5)
    # for dataset in df['Dataset'].unique():
    #     ax2.vlines(x=anon_group[dataset], ymin=0.8, ymax=2, color='lightgray', linestyle=':', linewidth=1.5)
    ax2.set_xlabel('Number of possible scores')
    ax2.set_ylabel('Median anonymity set size')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim(0.8, 30000)
    ax2.legend()
    plt.tight_layout()
    # plt.show()
    fig1.savefig('uniqueness.pdf', dpi=300, bbox_inches='tight')
    fig2.savefig('anonymity.pdf', dpi=300, bbox_inches='tight')


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
    plt.tight_layout()
    plt.show()


def prepare_data(parsed, keys):
    data = []
    for key in keys:
        data.extend([(key, value) for value in parsed[key]])
    return pd.DataFrame(data, columns=['Versus', 'Value'])


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


def load_ancestry_baselines(id_file, data_file, method):
    if id_file is not None:
        with open(id_file, "r") as f:
            ids = [line.strip() for line in f]
        ids = ids[1:]
        indf = pd.read_csv(data_file, sep='\t', header=None)
        indf.columns = ids
        indf.index = ids
    else:
        indf = pd.read_csv(data_file, sep='\t', index_col=0)

    pruned_df = indf.loc[indf.index.str.startswith('$'), ~indf.columns.str.startswith('$')]
    related = read_related_individuals()
    individual_to_ancestry = read_ancestry()

    ancestries = ['EUR', 'AFR', 'EAS', 'SAS', 'AMR']
    results = []
    for idv in pruned_df.index:
        target = idv[1:]  # Remove the '$' symbol
        if target in ancestries:
            continue
        anc = "$"+individual_to_ancestry[target]
        kinship = pruned_df.at[anc, target]
        self_rank = pruned_df.loc[anc].rank(ascending=False)[target]
        results.append({'Category': 'Self', 'Kinship': kinship, 'Rank': self_rank, 'Method': method})
        relatives = []
        category = ''
        if target in related:
            relatives = []
            for relative in related[target]:
                if relative[0] not in pruned_df.columns:
                    continue
                relative_kinship = pruned_df.at[anc, relative[0]]
                relative_rank = pruned_df.loc[anc].rank(ascending=False)[relative[0]]
                category = ''
                if relative[1] == '1':
                    category = '1st degree'
                elif relative[1] == '2':
                    category = '2nd degree'
                results.append({'Category': category, 'Kinship': relative_kinship, 'Rank': relative_rank, 'Method': method})
                relatives.append(relative[0])
        other_columns = pruned_df.columns[~pruned_df.columns.isin([target] + relatives)]
        kinships = pruned_df.loc[anc, other_columns]
        ranks = pruned_df.loc[anc].rank(ascending=False)[other_columns]
        temp_df = pd.DataFrame({
            'Category': 'Unrelated',
            'Kinship': kinships.values,
            'Rank': ranks.values,
            'Method': method
        })
        results.extend(temp_df.to_dict('records'))
        if len(relatives) > 0:
            temp_df = pd.DataFrame({
                'Category': f'{category} unrelated',
                'Kinship': kinships.values,
                'Rank': ranks.values,
                'Method': method
            })
            results.extend(temp_df.to_dict('records'))

    return pd.DataFrame(results)


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

    ancestries = ['EUR', 'AFR', 'EAS', 'SAS', 'AMR']
    first_degree_pairs = set()
    second_degree_pairs = set()
    results = []
    for idv in pruned_df.index:
        target = idv[1:]  # Remove the '$' symbol
        if target in ancestries:
            continue
        if target not in pruned_df.columns:
            continue
        kinship = pruned_df.at[idv, target]
        self_rank = pruned_df.loc[idv].rank(ascending=False)[target]
        results.append({'Category': 'Self', 'Kinship': kinship, 'Rank': self_rank, 'Method': method})
        results.append({'Category': 'Self', 'Kinship': truth_df.at[target, target],
                        'Rank': truth_df.loc[target].rank(ascending=False)[target], 'Method': 'True'+method})
        relatives = []
        category = ''
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
                    first_degree_pairs.add(tuple(sorted((target, relative[0]))))
                elif relative[1] == '2':
                    category = '2nd degree'
                    second_degree_pairs.add(tuple(sorted((target, relative[0]))))
                results.append({'Category': category, 'Kinship': relative_kinship, 'Rank': relative_rank, 'Method': method})
                results.append({'Category': category, 'Kinship': truth_df.at[target, relative[0]],
                                'Rank': truth_df.loc[target].rank(ascending=False)[relative[0]], 'Method': 'True'+method})
                relatives.append(relative[0])
        other_columns = pruned_df.columns[~pruned_df.columns.isin([target] + relatives)]
        kinships = pruned_df.loc[idv, other_columns]
        ranks = pruned_df.loc[idv].rank(ascending=False)[other_columns]
        temp_df = pd.DataFrame({
            'Category': 'Unrelated',
            'Kinship': kinships.values,
            'Rank': ranks.values,
            'Method': method
        })
        results.extend(temp_df.to_dict('records'))
        if len(relatives) > 0:
            temp_df = pd.DataFrame({
                'Category': f'{category} unrelated',
                'Kinship': kinships.values,
                'Rank': ranks.values,
                'Method': method
            })
            results.extend(temp_df.to_dict('records'))
        # Add the true values for unrelated individuals
        other_columns = truth_df.columns[~truth_df.columns.isin([target] + relatives)]
        kinships = truth_df.loc[target, other_columns]
        ranks = truth_df.loc[target].rank(ascending=False)[other_columns]
        temp_df = pd.DataFrame({
            'Category': 'Unrelated',
            'Kinship': kinships.values,
            'Rank': ranks.values,
            'Method': 'True'+method
        })
        results.extend(temp_df.to_dict('records'))

    print(f"Number of unique 1st degree pairs: {len(first_degree_pairs)}")
    print(first_degree_pairs)
    print(f"Number of unique 2nd degree pairs: {len(second_degree_pairs)}")
    print(second_degree_pairs)

    return pd.DataFrame(results)


def deanonymization_accuracy():
    params = {
        # 'font.size': 13,
        # 'axes.labelsize': 13,
        # 'xtick.labelsize': 13,
        # 'ytick.labelsize': 13,
        'legend.fontsize': LABEL_SIZE,
    }
    mpl.rcParams.update(params)

    king_df = load_results("ibd/plink.king.id", "ibd/plink.king", 'KING')
    gcta_df = load_results("ibd/plink.rel.id", "ibd/plink.rel", 'GCTA')
    gcta_df['Kinship'] = gcta_df['Kinship'] / 2 # Relatedness -> Kinship
    gcta_df['Kinship'] = [0.5 if k > 0.5 else k for k in gcta_df['Kinship']]
    # refined_df = load_results(None, "ibd/refined.ibd", 'Refined IBD')
    # df = pd.concat([king_df, gcta_df, refined_df])

    # Set cutoffs for accuracy levels
    cutoffs = {
        'Self': 0.3535,
        '1st degree': 0.1768,
        '2nd degree': 0.089,
    }
    palette = sns.color_palette("pastel")
    plot_kinship('KING', king_df, cutoffs)
    plot_kinship('GCTA', gcta_df, cutoffs)

    df = pd.concat([king_df, gcta_df])
    fig1, ax1 = plt.subplots(figsize=(7, 3))
    guess_methods = ['KING', 'GCTA']
    ground_truth_methods = ['TrueKING', 'TrueGCTA']
    method_df = df[df['Method'].isin(guess_methods)]
    method_df = method_df[method_df['Category'].isin(['Self', '1st degree', '2nd degree', 'Unrelated'])]
    sns.boxplot(x='Category', y='Kinship', data=method_df, hue='Method', palette='pastel', ax=ax1,
                order=['Self', '1st degree', '2nd degree', 'Unrelated'])

    box_width = 0.8 / len(ground_truth_methods)
    for method in ground_truth_methods:
        for category in ['Self', '1st degree', '2nd degree', 'Unrelated']:
            category_data = df[(df['Method'] == method) & (df['Category'] == category)]
            median_value = category_data['Kinship'].median()
            # Find the position of the category on the X axis
            category_pos = ['Self', '1st degree', '2nd degree', 'Unrelated'].index(category)
            method_position_offset = ground_truth_methods.index(method) * box_width
            x_position = category_pos - 0.4 + method_position_offset + box_width / 2
            print(f'{method} median: {median_value}, position: {category_pos}')
            ax1.plot([x_position - box_width / 2, x_position + box_width / 2], [median_value, median_value],
                     color=palette[3], linestyle='-')

    # Adding the legend for the red lines
    ax1.plot([], [], color=palette[3], linestyle='-', label='True-genotypes baseline')

    ax1.axhline(y=cutoffs['Self'], color='lightgray', linestyle='--', linewidth=1)
    ax1.text(x=ax1.get_xlim()[1] + 0.03, y=cutoffs['Self'], s='Self', va='center', ha='left', color='black', fontsize=LABEL_SIZE)
    ax1.axhline(y=cutoffs['1st degree'], color='lightgray', linestyle='--', linewidth=1)
    ax1.text(x=ax1.get_xlim()[1] + 0.03, y=cutoffs['1st degree'], s='1st', va='center', ha='left', color='black', fontsize=LABEL_SIZE)
    ax1.axhline(y=cutoffs['2nd degree'], color='lightgray', linestyle='--', linewidth=1)
    ax1.text(x=ax1.get_xlim()[1] + 0.03, y=cutoffs['2nd degree'], s='2nd', va='center', ha='left', color='black', fontsize=LABEL_SIZE)
    ax1.set_yticks([-0.5, -0.25, 0, 0.25, 0.5])
    ax1.set_yticklabels(['-0.5', '0.25', '0', '0.25', '0.5'])
    ax1.set_ylabel('Kinship')
    ax1.set_xlabel('')
    ax1.legend(title="")
    fig1.savefig('kinship_score.pdf', dpi=300, bbox_inches='tight')

    categories = ['Self', '1st degree', '2nd degree', 'Unrelated']
    plot_precision_recall(df, 'KING', categories, cutoffs)
    plot_precision_recall(df, 'GCTA', categories, cutoffs)

    df_king = df[df['Method'] == 'KING']
    fig4, ax4 = plt.subplots(figsize=(4, 3))
    colors = {'Self': 'black', '1st degree': 'coral', '2nd degree': 'green', 'Unrelated': 'purple'}
    patterns = {'Self': '-', '1st degree': 'dotted', '2nd degree': '--', 'Unrelated': '-.'}
    for category in categories:
        sns.kdeplot(df_king[df_king['Category'] == category]['Kinship'], label=category, common_norm=False,
                    ax=ax4, color=colors[category], linestyle=patterns[category])
    ax4.set_xlabel('Kinship coefficient')
    ax4.set_ylabel('Density')
    ax4.set_xlim(-0.4, 0.5)
    ax4.set_xticks([-0.25, 0, 0.25, 0.5])
    ax4.set_xticklabels(['-0.25', '0', '0.25', '0.5'])
    ax4.legend(title='', loc="upper left")
    plt.tight_layout()
    fig4.savefig('kinship_original.pdf', dpi=300, bbox_inches='tight')
    # plt.show()


def plot_kinship(method, df, cutoffs):
    categories = ['Self', '1st degree', '2nd degree', 'Unrelated']
    fig, ax = plt.subplots(figsize=(6, 3))
    sns.boxplot(x='Category', y='Kinship', data=df[(df['Method']==method) & (df['Category'].isin(categories))], ax=ax,
                legend=False, palette='pastel', order=categories, hue='Category')

    for category in categories:
        category_data = df[(df['Method'] == 'True'+method) & (df['Category'] == category)]
        if len(category_data) > 0:
            median_value = category_data['Kinship'].median()
            category_pos = ['Self', '1st degree', '2nd degree', 'Unrelated'].index(category)
            x_position = category_pos
            ax.plot([x_position - 0.4, x_position + 0.4], [median_value, median_value], color='red',
                    linestyle='-')
    ax.plot([], [], color='red', linestyle='-', label='True-genotypes baseline')

    ax.axhline(y=cutoffs['Self'], color='lightgray', linestyle='--', linewidth=1)
    ax.text(x=ax.get_xlim()[1] + 0.03, y=cutoffs['Self'], s='Self', va='center', ha='left', color='black', fontsize=LABEL_SIZE)
    ax.axhline(y=cutoffs['1st degree'], color='lightgray', linestyle='--', linewidth=1)
    ax.text(x=ax.get_xlim()[1] + 0.03, y=cutoffs['1st degree'], s='1st', va='center', ha='left', color='black', fontsize=LABEL_SIZE)
    ax.axhline(y=cutoffs['2nd degree'], color='lightgray', linestyle='--', linewidth=1)
    ax.text(x=ax.get_xlim()[1] + 0.03, y=cutoffs['2nd degree'], s='2nd', va='center', ha='left', color='black', fontsize=LABEL_SIZE)
    ax.set_yticks([-0.5, -0.25, 0, 0.25, 0.5])
    ax.set_yticklabels(['-0.5', '-0.25', '0', '0.25', '0.5'])
    ax.set_ylabel('Kinship')
    ax.set_xlabel('')
    ax.legend(title="", loc='lower left')
    fig.savefig(f'kinship_{method.lower()}.pdf', dpi=300, bbox_inches='tight')


def plot_precision_recall(df_orig, method, categories, cutoffs):
    df = df_orig[df_orig['Method'] == method]
    precision, recall = calculate_precision_and_recall(df, categories, cutoffs)
    if method == 'KING':
        baselines = load_ancestry_baselines("ibd/plink.king.id", "ibd/plink.king", method)
    elif method == 'GCTA':
        baselines = load_ancestry_baselines("ibd/plink.rel.id", "ibd/plink.rel", method)
    else:
        return
    baseline_precision, baseline_recall = calculate_precision_and_recall(baselines, categories, cutoffs)

    print(f"=== {method}: Precision and recall ===")
    print(precision)
    print(recall)
    print(f"=== {method}: Baseline precision and recall ===")
    print(baseline_precision)
    print(baseline_recall)

    fig2, ax2 = plt.subplots(figsize=(7, 4))
    palette = sns.color_palette("pastel")
    bar_width = 0.3
    index = np.arange(len(categories))

    precision_values = [precision[cat] for cat in categories]
    recall_values = [recall[cat] for cat in categories]
    ax2.bar(index, precision_values, bar_width, label='Precision', color=palette[0])
    ax2.bar(index + bar_width, recall_values, bar_width, label='Recall', color=palette[1])
    for i, category in enumerate(categories):
        plw, rlw = 5, 5
        if baseline_precision[category] > 0.01:
            plw = 3
        if baseline_recall[category] > 0.01:
            rlw = 3
        ax2.hlines(y=baseline_precision[category], xmin=index[i] - bar_width/2, xmax=index[i] + bar_width/2, color='red',
                   linestyle='-', linewidth=plw)
        ax2.hlines(y=baseline_recall[category], xmin=index[i] + bar_width - bar_width/2,
                   xmax=index[i] + bar_width + bar_width/2, color='red', linestyle='-', linewidth=rlw)
    ax2.plot([], [], color='red', linestyle='-', label='Predicting major genotypes')

    ax2.set_xlabel('')
    ax2.set_ylabel('')
    ax2.set_xticks(index + bar_width / 2)
    ax2.set_xticklabels(categories)
    ax2.legend(title="", loc='lower right')
    plt.tight_layout()
    fig2.savefig(f'kinship_precision_{method.lower()}.pdf', dpi=300, bbox_inches='tight')


def calculate_precision_and_recall(df, categories, cutoffs):
    precision = {}
    recall = {}
    for category in categories:
        if category == 'Unrelated':
            true_positive = len(df[(df['Category'] == category) & (df['Kinship'] < cutoffs['2nd degree'])])
            false_positive = len(df[(df['Category'].isin(['Self', '1st degree', '2nd degree'])) & (df['Kinship'] < cutoffs['2nd degree'])])
            false_negative = len(df[(df['Category'] == category) & (df['Kinship'] >= cutoffs['2nd degree'])])
        else:
            true_positive = len(df[(df['Category'] == category) & (df['Kinship'] >= cutoffs[category])])
            false_negative = len(df[(df['Category'] == category) & (df['Kinship'] < cutoffs[category])])
            if category == 'Self':
                false_positive = len(df[(df['Category'].isin(['1st degree', '2nd degree', 'Unrelated'])) & (df['Kinship'] >= cutoffs[category])])
            if category == '1st degree':
                false_positive = len(df[(df['Category'] == '1st degree unrelated') & (df['Kinship'] >= cutoffs[category])])
            if category == '2nd degree':
                false_positive = len(df[(df['Category'] == '2nd degree unrelated') & (df['Kinship'] >= cutoffs[category])])
        precision[category] = true_positive / (true_positive + false_positive) if (true_positive + false_positive) > 0 else 0
        recall[category] = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else 0
    return precision, recall


def read_ancestry():
    with open('data/superpopulations.json', 'r') as file:
        data = json.load(file)
    for idv, anc in data.items():
        if "," in anc:
            data[idv] = anc.split(',')[0]
    return data


def plot_defense(pgs):
    params = {
        # 'axes.labelsize': 12,
        # 'xtick.labelsize': 12,
        # 'ytick.labelsize': 12,
        'legend.fontsize': 13,
        'legend.title_fontsize': 13,
    }
    mpl.rcParams.update(params)
    directory = "results/defense/"
    filepath = os.path.join(directory, pgs+".json")
    with open(filepath, 'r') as f:
        content = json.load(f)
    data = []
    for result in content:
        data.append({'Precision': int(result['Precision']), 'Density': float(result['Density']),
                     'TotalPossibleScores': int(result['TotalPossibleScores']),
                     'PercentageUnique': float(result['PercentageUnique']),
                     'AnonSize': int(result['MedianAnonymitySetSize'])})
    df = pd.DataFrame(data)

    fig1, ax1 = plt.subplots(figsize=(5, 3))
    sns.lineplot(x='Precision', y='PercentageUnique', data=df, color='coral', marker='.', ax=ax1)
    ax1.set_ylim(0, 105)
    plt.gca().invert_xaxis()
    ax1.set_xticks([p for p in range(17, 0, -2)])
    ax1.set_xlabel('Weight precision, digits')
    ax1.set_ylabel(r'Identifiable individuals, \%')
    ax1.spines['right'].set_visible(False)
    ax2 = ax1.twinx()
    sns.lineplot(x='Precision', y='AnonSize', data=df, color='coral', marker='v', linestyle="--", ax=ax2)
    ax2.set_ylabel('Median Anonymity Size')
    ax2.spines['right'].set_linestyle((0, (5, 5)))
    ax2.set_yscale('log')
    ax2.yaxis.set_minor_locator(plt.NullLocator())
    identifiability_legend = Line2D([0], [0], color='coral', marker='.', label='Identifiability')
    anonymity_legend = Line2D([0], [0], color='coral', marker='v', linestyle="--", label='Anonymity')
    ax1.legend(handles=[identifiability_legend, anonymity_legend], title='UKBB', loc='center left')
    plt.tight_layout()
    fig1.savefig('precision_uniqueness.pdf', dpi=300, bbox_inches='tight')

    fig2, ax3 = plt.subplots(figsize=(4, 3))
    sns.lineplot(x='Precision', y='Density', data=df, color='coral', marker='.', ax=ax3, label='UKBB')
    ax3.set_xticks([p for p in range(17, 0, -2)])
    ax3.set_yticks([p for p in range(0, 30, 5)])
    ax3.set_xlabel('Weight precision, digits')
    ax3.set_ylabel('Density')
    ax3.legend()
    plt.gca().invert_xaxis()
    plt.tight_layout()
    fig2.savefig('precision_density.pdf', dpi=300, bbox_inches='tight')


# def plot_defended_distribution():
#     filepath = "results/defense/distro.json"
#     with open(filepath, 'r') as f:
#         content = json.load(f)
#     data = []
#     for precision, scores in content.items():
#         for score, count in scores.items():
#             data.append({'Precision': int(precision), 'Score': float(score), 'Count': int(count)})
#     df = pd.DataFrame(data)
#     filtered_df = df[df['Precision'].isin([4, 17])]
#
#     # Plot the distribution of the scores
#     plt.figure(figsize=(4, 3))
#     sns.histplot(data=filtered_df, x='Score', weights='Count', hue='Precision', multiple='dodge', bins=100)
#     plt.xlabel('PRS')
#     plt.ylabel('Count')
#     plt.tick_params(axis='y', which='both', labelleft=False, left=False, right=False)
#     plt.xlim(0.10, 0.17)
#     plt.tight_layout()
#     # plt.savefig('precision_distribution.pdf', dpi=300, bbox_inches='tight')
#     plt.show()

def plot_defended_distribution():
    params = {
        'legend.fontsize': 13,
        'legend.title_fontsize': 13,
    }
    mpl.rcParams.update(params)
    filepath = "results/defense/scores.json"
    with open(filepath, 'r') as f:
        content = json.load(f)
    data = {}
    for precision, scores in content.items():
        data[int(precision)] = [float(x) for x in scores]

    # Plot the scatter plot
    plt.figure(figsize=(4, 3))
    linestyles = ['-', '--', '-.', ':']
    for i, precision in enumerate([1, 2, 5, 17]):
        sns.kdeplot(data=data[precision], label=f'{precision}', linestyle=linestyles[i])
    # sns.scatterplot(x=data[17], y=data[1], color='coral', alpha=0.5)
    plt.xlabel('PRS')
    plt.ylabel('Count')
    plt.legend(title='Precision')
    plt.tight_layout()
    plt.savefig('precision_distribution.pdf', dpi=300, bbox_inches='tight')
    # plt.show()


# def prs_prediction():
#     filepath = "results/predict/prediction2.json"
#     data = []
#     with open(filepath, 'r') as f:
#         content = json.load(f)
#     for prs_id, results in content.items():
#         for i, prd in enumerate(results['Predicted']):
#             data.append({'PRS': prs_id, 'Known': float(results["Known"]), 'Total': float(results["Total"]),
#                          'Max': float(results["Max"]), 'Min': float(results["Min"]),
#                          'Mean': float(results["Means"][results["Ancestries"][i]]),
#                          'Std': float(results['Stds'][results["Ancestries"][i]]),
#                          'Ancestry': results["Ancestries"][i],
#                          'Predicted': float(results['Predicted'][i]), 'Real': float(results['Real'][i])})
#
#     df = pd.DataFrame(data)
#     print(df['PRS'].nunique())
#     pd.set_option('display.max_rows', None)
#     pd.set_option('display.max_columns', None)
#     bin_size = 0.05
#     df['Ratio'] = (df['Known'] / df['Total']).round(3)
#     df['Binned_Ratio'] = (bin_size * np.round(df['Ratio'] / bin_size)).round(3)
#     df['Threshold'] = df['Mean'] + df['Std']
#     df['Above_Threshold'] = df['Real'] > df['Threshold']
#     df['Predicted_Above_Threshold'] = df['Predicted'] > df['Threshold']
#
#     df['Real_Percentile'] = norm.cdf(df['Real'], loc=df['Mean'], scale=df['Std'])
#     df['Predicted_Percentile'] = norm.cdf(df['Predicted'], loc=df['Mean'], scale=df['Std'])
#     df['Percentile_Diff'] = abs(df['Real_Percentile'] - df['Predicted_Percentile']) * 100
#
#     def calculate_metrics(group):
#         tp = np.sum((group['Above_Threshold'] == True) & (group['Predicted_Above_Threshold'] == True))
#         fp = np.sum((group['Above_Threshold'] == False) & (group['Predicted_Above_Threshold'] == True))
#         fn = np.sum((group['Above_Threshold'] == True) & (group['Predicted_Above_Threshold'] == False))
#
#         precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
#         recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
#         std_dev = group['Predicted'].std()
#         percentile_diff = group['Percentile_Diff'].median()
#         percentile_diff_std = group['Percentile_Diff'].std()
#         high_risk_percentile_diff = group[group['Real'] > group['Threshold']]['Percentile_Diff'].median()
#         high_risk_percentile_diff_std = group[group['Real'] > group['Threshold']]['Percentile_Diff'].std()
#
#         return pd.Series({'Precision': precision, 'Recall': recall, 'Std_Dev': std_dev,
#                           'Mean_Percentile_Diff': percentile_diff, 'High_Risk_Percentile_Diff': high_risk_percentile_diff,
#                           'Percentile_std': percentile_diff_std, 'High_Risk_Percentile_Std': high_risk_percentile_diff_std})
#
#     grouped = df.groupby('Binned_Ratio').apply(calculate_metrics).reset_index()
#
#     plt.figure(figsize=(4, 3))
#     sns.lineplot(data=grouped, x='Binned_Ratio', y='Precision', label='Precision', color='coral')
#     sns.lineplot(data=grouped, x='Binned_Ratio', y='Recall', label='Recall', color='teal')
#     plt.fill_between(grouped['Binned_Ratio'], grouped['Precision'] - grouped['Std_Dev'],
#                      grouped['Precision'] + grouped['Std_Dev'], color='teal', alpha=0.1)
#     plt.fill_between(grouped['Binned_Ratio'], grouped['Recall'] - grouped['Std_Dev'],
#                      grouped['Recall'] + grouped['Std_Dev'], color='coral', alpha=0.1)
#     plt.xlabel('Guessed SNPs ratio')
#     plt.ylabel('')
#     plt.gca().invert_xaxis()
#     plt.ylim(0, 1)
#     plt.legend(loc='lower left')
#     # plt.savefig('prediction_highrisk_by_population.pdf', dpi=300, bbox_inches='tight')
#     # plt.show()
#
#     plt.figure(figsize=(4, 3))
#     sns.lineplot(data=grouped, x='Binned_Ratio', y='Mean_Percentile_Diff', label='All individuals', linestyle='-',
#                  color='teal')
#     # plt.fill_between(grouped['Binned_Ratio'], grouped['Mean_Percentile_Diff'] - grouped['Percentile_std'],
#     #                  grouped['Mean_Percentile_Diff'] + grouped['Percentile_std'], color='teal', alpha=0.1)
#     sns.lineplot(data=grouped, x='Binned_Ratio', y='High_Risk_Percentile_Diff', label=r'High risk ($>\mu$+$\sigma$)',
#                  linestyle='--', color='coral')
#     # plt.fill_between(grouped['Binned_Ratio'], grouped['High_Risk_Percentile_Diff'] - grouped['High_Risk_Percentile_Std'],
#     #                  grouped['High_Risk_Percentile_Diff'] + grouped['High_Risk_Percentile_Std'], color='coral', alpha=0.1)
#     plt.xlabel('Guessed SNPs ratio')
#     plt.ylabel('Median percentile error')
#     plt.legend(loc='upper left')
#     plt.ylim(0, 50)
#     plt.xticks(np.arange(0.1, 0.9, 0.1))
#     plt.gca().invert_xaxis()
#     # plt.show()
#     plt.savefig('prediction_median_percentile.pdf', dpi=300, bbox_inches='tight')


def prs_prediction():
    filepath = "results/predict/prediction2.json"
    data = []
    with open(filepath, 'r') as f:
        content = json.load(f)
    for prs_id, results in content.items():
        for i, prd in enumerate(results['Predicted']):
            data.append({'PRS': prs_id, 'Known': int(results["Known"]), 'Total': int(results["Total"]),
                         'Known/Total': f'{results["Known"]}/{results["Total"]}',
                         'Trait': results["Trait"], 'Max': float(results["Max"]), 'Min': float(results["Min"]),
                         'Mean': float(results["Means"][results["Ancestries"][i]]),
                         'Std': float(results['Stds'][results["Ancestries"][i]]),
                         'Ancestry': results["Ancestries"][i], 'Correct': float(results['CorrectSnps'][i]),
                         'Predicted': float(results['Predicted'][i]), 'Real': float(results['Real'][i])})

    df = pd.DataFrame(data)
    # pd.set_option('display.max_rows', None)
    # pd.set_option('display.max_columns', None)
    # results = []
    # for prs_id, group in df.groupby('PRS'):
    #     # Use the Mean and Std from the data
    #     group['HighRiskThreshold'] = group['Mean'] + group['Std']
    #
    #     high_risk_indices = group[group['Real'] >= group['HighRiskThreshold']].index
    #
    #     true_positives = np.sum(group.loc[high_risk_indices, 'Predicted'] >= group.loc[high_risk_indices, 'HighRiskThreshold'])
    #     false_positives = np.sum(group['Predicted'] >= group['HighRiskThreshold']) - true_positives
    #     false_negatives = len(high_risk_indices) - true_positives
    #
    #     precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    #     recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    #
    #     results.append({
    #         'PRS': prs_id,
    #         'Trait': group['Trait'].iloc[0],
    #         'Known/Total': f"{group['Known'].iloc[0]}/{group['Total'].iloc[0]}",
    #         'Precision': precision,
    #         'Recall': recall
    #     })
    #
    # df = pd.DataFrame(results)
    df['Ratio'] = df['Known/Total'].apply(lambda x: float(x.split('/')[0])*100 / float(x.split('/')[1]))
    # Calculate percentiles and percentile differences
    df['Real_Percentile'] = norm.cdf(df['Real'], loc=df['Mean'], scale=df['Std'])
    df['Predicted_Percentile'] = norm.cdf(df['Predicted'], loc=df['Mean'], scale=df['Std'])
    df['Percentile_Diff'] = abs(df['Real_Percentile'] - df['Predicted_Percentile']) * 100
    # results_df_melted = df.melt(id_vars=['PRS', 'Trait', 'Ratio'], value_vars=['Precision', 'Recall'], var_name='Metric', value_name='Value')

    df['HighRisk'] = df['Real'] > (df['Mean'] + df['Std'])
    df['RiskCategory'] = np.where(df['HighRisk'], 'High Risk', 'All')

    plt.figure(figsize=(9, 5))
    # sns.boxplot(x='Trait', y='Value', hue='Metric', data=results_df_melted)
    sns.boxplot(x='Trait', y='Correct', data=df)
    # sns.boxplot(x='Trait', y='Percentile_Diff', data=df, hue='RiskCategory')
    xtick_labels = df.groupby('Trait').apply(lambda x: '\n'.join([f"{prs} ({kt})" for prs, kt in
                                                                  zip(x['PRS'].unique(), x['Known/Total'].unique())])).reset_index()
    xtick_labels = xtick_labels.set_index('Trait')[0].to_dict()
    plt.xticks(ticks=plt.xticks()[0], labels=[f"{trait}\n{xtick_labels[trait]}" for trait in xtick_labels])
    plt.legend(title="")
    # plt.ylabel('Percentile difference')
    plt.ylabel('Ratio of correctly predicted SNPs')
    plt.xlabel('')
    # plt.ylim(0, 1)
    plt.tight_layout()
    # plt.savefig('prediction_percentile_box.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('prediction_snps.pdf', dpi=300, bbox_inches='tight')
    # plt.show()


LABEL_SIZE = 15
TICK_SIZE = 15
LEGEND_SIZE = 15
def prepare_for_latex():
    params = {
        'axes.labelsize': LABEL_SIZE,
        'xtick.labelsize': TICK_SIZE,
        'ytick.labelsize': TICK_SIZE,
        'legend.fontsize': LEGEND_SIZE,
        'legend.title_fontsize': LEGEND_SIZE,
        'text.usetex': True,
        # 'font.family': 'serif',
        # 'font.serif': 'Computer Modern',
        'font.family': 'sans-serif',
        'font.serif': 'Helvetica',
        'savefig.bbox': 'tight',
        'savefig.format': 'pdf',
    }
    mpl.rcParams.update(params)


if __name__ == "__main__":
    prepare_for_latex()
    # score_uniqueness()
    # plot_defense("PGS000869")
    # plot_defended_distribution()
    # sequential_idv_accuracy()
    # sequential_loci_accuracy()
    # af_hist()
    # plot_normal_distribution_with_fill()
    deanonymization_accuracy()
    # prs_prediction()
import argparse
import csv
import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns


FIGURE_FOLDER = 'analysis/figures/'

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
        'font.family': 'sans-serif',
        'font.serif': 'Helvetica',
        'savefig.bbox': 'tight',
        'savefig.format': 'pdf',
    }
    mpl.rcParams.update(params)


def load_idv_accuracy(filename):
    directory = "results/recoveryAccuracy/"
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


def solving_idv_accuracy():
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
    # fig.savefig(FIGURE_FOLDER + 'recovery.pdf', dpi=300, bbox_inches='tight')
    plt.show()


def solving_loci_accuracy():
    directory = "results/recoveryAccuracy/"
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
    fig1.savefig(FIGURE_FOLDER + 'accuracy_eaf.pdf', dpi=300, bbox_inches='tight')

    fig2, ax2 = plt.subplots(figsize=(4, 3))
    gwas_df = df[df['GWAS'] != 0]
    sns.kdeplot(data=df, x='EAF', hue='Ancestry', ax=ax2, common_norm=False)
    sns.kdeplot(data=gwas_df, x='GWAS', ax=ax2, linestyle='--', color='black', common_norm=False, label='GWAS')
    plt.xlabel('Effect Allele Frequency')
    plt.ylabel('Density')
    plt.xlim(0, 1)
    ancestries = df['Ancestry'].unique()
    legend_labels = [*ancestries, "GWAS"]
    legend_colors = sns.color_palette()[:len(ancestries)] + ['black']
    handles = [Line2D([0, 1], [0, 1], color=color) for color in legend_colors[:-1]]
    gwas_handle = Line2D([0, 1], [0, 1], color='black', linestyle='--')
    handles.append(gwas_handle)
    ax2.legend(handles, legend_labels, ncol=2, loc="lower left", columnspacing=0.5)
    # ax2.legend(ncol=2, loc="lower left")
    plt.xticks([0, 0.25, 0.5, 0.75, 1], ['0', '0.25', '0.5', '0.75', '1'])
    plt.tight_layout()
    fig2.savefig(FIGURE_FOLDER + 'accuracy_ancestry.pdf', dpi=300, bbox_inches='tight')
    # plt.show()


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

    bins = []
    for i in range(1, int(np.log10(df['TotalPossibleScores'].max())) + 2):
        bins.append(10**(i-1))
        bins.append(10**(i-1) * 10**0.5)
    bins.append(10**i)
    df['TotalPossibleScoresBin'] = pd.cut(df['TotalPossibleScores'], bins=bins, right=False)
    grouped_real = df.groupby(['TotalPossibleScoresBin', 'Dataset'], observed=True).agg(
        RealPercentageUnique_median=('RealPercentageUnique', 'median'),
        RealPercentageUnique_std=('RealPercentageUnique', 'std'),
        RealAnonSize_median=('RealAnonSize', 'median'),
        RealAnonSize_std=('RealAnonSize', 'std'),
        TotalPossibleScores_mean=('TotalPossibleScores', 'mean'),
        NumVariants_median=('NumVariants', 'median')
    ).reset_index()

    df['PredictedPercentageUnique_filtered'] = df['PredictedPercentageUnique'].where(df['PredictedPercentageUnique'] > 0)
    grouped_predicted = df.groupby([
        'TotalPossibleScoresBin', 'Dataset'], observed=True).agg(
        PredictedPercentageUnique_median=('PredictedPercentageUnique_filtered', 'median'),
        PredictedPercentageUnique_std=('PredictedPercentageUnique_filtered', 'std'),
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

        print(f"{dataset} - Mean TotalPossibleScores at 95% unique: {num_unique[dataset]}")
        print(f"{dataset} - Mean TotalPossibleScores when RealAnonSize reaches 2: {anon_group[dataset]}")

    first_color, second_color = 'teal', 'coral'
    first_marker, second_marker = '.', '*'
    observed, estimated = 'observed', 'estimated'

    fig1, ax1 = plt.subplots(figsize=(5, 4))
    for dataset, color, label_prefix, marker in [('1000Genomes', first_color, '1000G', first_marker),
                                                 ('UKBiobank', second_color, 'UKBB', second_marker)]:
        subset_real = grouped_real[grouped_real['Dataset'] == dataset]
        x_real = subset_real['TotalPossibleScoresBin'].apply(lambda x: x.mid)
        y_real = subset_real['RealPercentageUnique_median']
        y_real_std = np.clip(subset_real['RealPercentageUnique_std'], 0, 100 - y_real)

        ax1.plot(x_real, y_real, color=color, label=f'{label_prefix} {observed}', marker=marker)
        ax1.fill_between(x_real, y_real - y_real_std, y_real + y_real_std, color=color, alpha=0.1)

        subset_pred = grouped_predicted[grouped_predicted['Dataset'] == dataset]
        if not subset_pred.empty:
            x_pred = subset_pred['TotalPossibleScoresBin'].apply(lambda x: x.mid)
            y_pred = subset_pred['PredictedPercentageUnique_median']
            ax1.plot(x_pred, y_pred, color=color, linestyle=':', label=f'{label_prefix} {estimated}', marker=marker)

    ax1.text(x=0.2, y=93, s='95', va='center', ha='right', color='gray', fontsize=TICK_SIZE)
    ax1.axhline(y=95, color='lightgray', linestyle='--', linewidth=1.5)
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

        subset_anon_pred = grouped_predicted[grouped_predicted['Dataset'] == dataset]
        if not subset_anon_pred.empty:
            x_anon_pred = subset_anon_pred['TotalPossibleScoresBin'].apply(lambda x: x.mid)
            y_pred_anon = subset_anon_pred['PredictedAnonSize_median']
            ax2.plot(x_anon_pred, y_pred_anon, color=color, linestyle=':', label=f'{label_prefix} {estimated}',
                     marker=marker)

    ax2.text(x=0.2, y=2, s='2', va='center', ha='right', color='gray', fontsize=TICK_SIZE)
    ax2.axhline(y=2, color='lightgray', linestyle='--', linewidth=1.5)
    ax2.set_xlabel('Number of possible scores')
    ax2.set_ylabel('Median anonymity set size')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim(0.8, 30000)
    ax2.legend()
    plt.tight_layout()
    # plt.show()
    fig1.savefig(FIGURE_FOLDER + 'uniqueness.pdf', dpi=300, bbox_inches='tight')
    fig2.savefig(FIGURE_FOLDER + 'anonymity.pdf', dpi=300, bbox_inches='tight')


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


def linking():
    params = {
        'legend.fontsize': LABEL_SIZE,
    }
    mpl.rcParams.update(params)

    king_df = load_results("results/linking/plink.king.id", "results/linking/plink.king", 'KING')
    gcta_df = load_results("results/linking/plink.rel.id", "results/linking/plink.rel", 'GCTA')
    gcta_df['Kinship'] = gcta_df['Kinship'] / 2 # Relatedness -> Kinship
    gcta_df['Kinship'] = [0.5 if k > 0.5 else k for k in gcta_df['Kinship']]
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
    fig1.savefig(FIGURE_FOLDER + 'kinship_joint.pdf', dpi=300, bbox_inches='tight')

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
    fig4.savefig(FIGURE_FOLDER + 'king_original.pdf', dpi=300, bbox_inches='tight')
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
    fig.savefig(FIGURE_FOLDER + f'kinship_{method.lower()}.pdf', dpi=300, bbox_inches='tight')


def plot_precision_recall(df_orig, method, categories, cutoffs):
    df = df_orig[df_orig['Method'] == method]
    precision, recall = calculate_precision_and_recall(df, categories, cutoffs)
    if method == 'KING':
        baselines = load_ancestry_baselines("results/linking/plink.king.id", "results/linking/plink.king", method)
    elif method == 'GCTA':
        baselines = load_ancestry_baselines("results/linking/plink.rel.id", "results/linking/plink.rel", method)
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
    fig2.savefig(FIGURE_FOLDER + f'kinship_precision_{method.lower()}.pdf', dpi=300, bbox_inches='tight')


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


def plot_rounding(pgs):
    params = {
        'legend.fontsize': 13,
        'legend.title_fontsize': 13,
    }
    mpl.rcParams.update(params)
    directory = "results/rounding/"
    filepath = os.path.join(directory, pgs+"_stats.json")
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
    fig1.savefig(FIGURE_FOLDER + 'precision_uniqueness.pdf', dpi=300, bbox_inches='tight')
    #
    fig3, ax3 = plt.subplots(figsize=(4, 3))
    sns.lineplot(x='Precision', y='Density', data=df, color='coral', marker='.', ax=ax3, label='UKBB')
    ax3.set_xticks([p for p in range(17, 0, -2)])
    ax3.set_yticks([p for p in range(0, 30, 5)])
    ax3.set_xlabel('Weight precision, digits')
    ax3.set_ylabel('Density')
    ax3.legend()
    plt.gca().invert_xaxis()
    plt.tight_layout()
    fig3.savefig(FIGURE_FOLDER + 'precision_density.pdf', dpi=300, bbox_inches='tight')
    #
    filepath = os.path.join(directory, pgs+"_scores.json")
    with open(filepath, 'r') as f:
        content = json.load(f)
    data = {}
    for precision, scores in content.items():
        data[int(precision)] = [float(x) for x in scores]
    fig4, ax4 = plt.subplots(figsize=(4, 3))
    linestyles = ['-', '--', '-.', ':']
    for i, precision in enumerate([1, 2, 5, 17]):
        sns.kdeplot(data=data[precision], label=f'{precision}', linestyle=linestyles[i])
    ax4.set_xlabel('PRS')
    ax4.set_ylabel('Count')
    ax4.legend(title='Precision')
    plt.tight_layout()
    fig4.savefig(FIGURE_FOLDER + 'precision_distribution.pdf', dpi=300, bbox_inches='tight')


def plot_snps_to_scores():
    filepath = "results/uniqueness/num_snps.json"
    with open(filepath, 'r') as f:
        content = json.load(f)
    data = []
    for row in content:
        data.append({'NumVariants': int(row['NumVariants']), 'Scores': int(row['TotalPossibleScores']),
                     'Precision': int(row['Precision'])})
    df = pd.DataFrame(data)
    grouped_df = df.groupby('NumVariants', observed=False).agg(
        Scores_median=('Scores', 'median'),   # Median of Scores for each NumVariants
        Scores_std=('Scores', 'std')          # Standard deviation of Scores for each NumVariants
    ).reset_index()
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.lineplot(x='NumVariants', y='Scores_median', data=grouped_df, ax=ax, color='#FFD166')
    # ax.fill_between(grouped_df['NumVariants'],
    #                 grouped_df['Scores_median'] - grouped_df['Scores_std'],
    #                 grouped_df['Scores_median'] + grouped_df['Scores_std'],
    #                 color='#FFD166', alpha=0.2)

    plt.ylabel('Total possible scores')
    plt.xlabel('Number of SNPs')
    ax.set_yscale('log')
    ax.set_yticks([10**i for i in [3, 7, 11, 15, 19]])
    ax.set_xticks([i for i in range(5, 50, 10)])
    plt.tight_layout()
    fig.savefig(FIGURE_FOLDER + 'snps_to_scores.pdf', dpi=300, bbox_inches='tight')
    # plt.show()


def plot_pgs_distribution():
    filepath = "results/validated_pgs.json"
    with open(filepath, 'r') as f:
        content = json.load(f)
    data = []
    for pgs, num_snps in content.items():
        data.append({'SNPs': int(num_snps)})
    df = pd.DataFrame(data)
    df['SNPs'] = pd.to_numeric(df['SNPs'])

    fig, ax = plt.subplots(figsize=(4, 3))
    sns.histplot(df['SNPs'], bins=47, kde=False, color='#FFD166', ax=ax)
    plt.ylabel('Number of PGS')
    plt.xlabel('Number of SNPs')
    ax.set_xlim(0, 50)
    ax.set_xticks([2, 10, 20, 30, 40, 50])
    ax.set_xticklabels([2, 10, 20, 30, 40, 50])
    plt.tight_layout()
    fig.savefig(FIGURE_FOLDER + 'pgs_distribution.pdf', dpi=300, bbox_inches='tight')
    # plt.show()


if __name__ == "__main__":
    # -----------Argument Parser-------------
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-e",
        "--expr",
        type=str,
        help="experiment to plot: recovery, linking, uniqueness, rounding, all",
        required=True)
    args = parser.parse_args()
    EXPR = args.expr
    prepare_for_latex()
    if EXPR == "recovery":
        solving_idv_accuracy()
        solving_loci_accuracy()
    elif EXPR == "linking":
        linking()
    elif EXPR == "uniqueness":
        score_uniqueness()
        plot_snps_to_scores()
    elif EXPR == "rounding":
        plot_rounding("PGS000869")
    elif EXPR == "all":
        solving_idv_accuracy()
        solving_loci_accuracy()
        linking()
        score_uniqueness()
        plot_snps_to_scores()
        plot_rounding("PGS000869")
    else:
        print("Unknown experiment: choose between the available options")
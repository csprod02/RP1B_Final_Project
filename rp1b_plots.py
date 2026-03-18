import pandas as pd
import numpy as np

import seaborn as sns
from matplotlib import pyplot as plt

from sklearn.metrics.pairwise import cosine_similarity


def read_data():
    
    typh_single_base = pd.read_csv('snp_count_df_typh.csv')
    typh_context = pd.read_csv('context_count_df_typh.csv')

    aga_single_base = pd.read_csv('snp_count_df_agalactiae.csv')
    aga_context = pd.read_csv('context_count_df_agalactiae.csv')

    epiderm_single_base = pd.read_csv('snp_count_df_epidermidis.csv')
    epiderm_context = pd.read_csv('context_count_df_epidermidis.csv')

    return typh_single_base, typh_context, aga_single_base, aga_context, epiderm_single_base, epiderm_context


def plot_single_base(typh_df, aga_df, epi_df):
    
    palette = {
        'transversion': '#FF8F00',
        'transition': '#26A69A'
        }

    species_dfs = [typh_df, aga_df, epi_df]
    species_names = ['S. typhimurium', 'S. agalactiae', 'S. epidermidis']

    metrics = ['count', 'norm', 'frac']

    label_map = {
        'count': 'Raw count',
        'norm': 'Normalised rate',
        'frac': 'Standardised rate'
    }

    fig, axes = plt.subplots(3, 3, figsize=(20, 14))
    axes = axes.reshape(3, 3)

    
    for row, metric in enumerate(metrics):

        row_max = max(df[metric].max() for df in species_dfs)
        
        for col, (df, species) in enumerate(zip(species_dfs, species_names)):
            
            ax = axes[row, col]

            sns.barplot(
                data=df,
                x='collapsed_group',
                y=metrics[row],
                hue='mutation_type',
                palette=palette,
                ax=ax
            )

            for container in ax.containers:
                ax.bar_label(container, fmt="%.2f", fontsize=8)

       
            ax.set_ylim(0, row_max * 1.05)   # 5% headroom


            # Titles for species
            if row == 0:
                ax.set_title(species)

            # Row labels
            if col == 0:
                ax.set_ylabel(label_map[metric])
            else:
                ax.set_ylabel("")

            # X labels only on bottom row
            if row == 2:
                ax.tick_params(axis="x", rotation=45)
            else:
                ax.set_xlabel("")

            if ax.get_legend():
                ax.get_legend().remove()


    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels,
               title='Mutation Type',
               loc="right",
               bbox_to_anchor=(1.02, 1)
               )

    plt.tight_layout()

    
    plt.savefig('mutation_grid.png', dpi=300, bbox_inches='tight')

    return


def plot_context(typh_df, aga_df, epi_df):

    palette = [
        "#C62828",  
        "#FF8F00",  
        "#FBC02D",  
        "#4CAF50",  
        "#26A69A",
        "#3F51B5"
    ]

    all_mutations = pd.concat([typh_df, aga_df, epi_df])['sb_mutation'].unique()

    color_dict = dict(zip(all_mutations, palette))

    species_list = [typh_df, aga_df, epi_df]
    species_names = ['S. typhimurium', 'S. agalactiae', 'S. epidermidis']

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(25, 17), sharex=True, sharey=True)

    linewidth_dict = {False: 1, True: 3}  # thicker edges for transversions
    alpha_dict = {False: 1.0, True: 0.8}  # slightly transparent for transversions

    
    for ax, df_sp, name in zip(axes, species_list, species_names):
        # Sort so bars are in the same order
        df_sp['sort_key'] = df_sp['sb_mutation'] + "_" + df_sp['context_collapsed_group']
        df_sp.sort_values(['sb_mutation', 'context_collapsed_group'], inplace=True)

        
        bars = sns.barplot(
            x='context_collapsed_group', 
            y='frac', 
            data=df_sp,
            hue='sb_mutation',
            palette=color_dict,
            dodge=False,  # one bar per mutation
            ax=ax,
            order=df_sp['context_collapsed_group']
        )

        for bar, ttv in zip(bars.patches, df_sp['flanking_equal']):
            bar.set_edgecolor('black')
            bar.set_linewidth(linewidth_dict[ttv])
#            bar.set_alpha(alpha_dict[ttv])

        
        ax.set_title(name, fontsize=14)
        ax.set_ylabel('Standardised Mutation Rate')
        ax.set_xlabel('')  # leave blank for all except bottom
        ax.legend_.remove()  # remove legend for individual axes

    # Single legend below all plots
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, title='Central Mutation', ncol=len(all_mutations),
               bbox_to_anchor=(0.5, -0.02), loc='upper center')

    plt.xticks(rotation=90)
    plt.tight_layout()

    plt.savefig('context_plot.png', dpi=300, bbox_inches='tight')





def plot_trinucleotide_grid(df, species_name):

    palette = [
        "#C62828",  
        "#FF8F00",  
        "#FBC02D",  
        "#4CAF50",  
        "#26A69A",
        "#3F51B5"
    ]


    bases = ['A','C','G','T']  # for rows and columns
    
    fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(20,16), sharey=True)
    
    # Get color palette for sb_mutation
    mutations = df['sb_mutation'].unique()
    color_dict = dict(zip(mutations, palette))
    
    for i, three in enumerate(bases):
        for j, five in enumerate(bases):
            ax = axes[i, j]
            
            # subset dataframe for this trinucleotide context
            subset = df[(df['fivePrime']==five) & (df['threePrime']==three)]
            
            if len(subset) == 0:
                ax.axis('off')  # hide empty plots
                continue
            
            sns.barplot(
                x='context_collapsed_group',
                y='frac',
                data=subset,
                hue='sb_mutation',
                palette=color_dict,
                dodge=False,
                ax=ax
            )
            

            if i == 0:
                ax.set_title(f"5′ = {five}", fontsize=12)
            
            # row labels only on first column
            if j == 0:
                ax.set_ylabel(f"3′ = {three}", fontsize=12)
            else:
                ax.set_ylabel('')
            
            # only show x tick labels on bottom row
            if i == len(bases) - 1:
                ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
            else:
                ax.set_xticklabels([])
            
            ax.set_xlabel('')
            ax.legend_.remove()
    
    # Single legend for the whole figure
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, title='Central Mutation', ncol=len(mutations),
               bbox_to_anchor=(0.5, -0.02), loc='upper center')
    
    plt.tight_layout()
    plt.suptitle(species_name, fontsize=16)
    plt.subplots_adjust(top=0.92, bottom=0.08)

    plt.savefig(f'context_{species_name}.png', dpi=300, bbox_inches='tight')

    return



def create_df_for_cosine(typh_df, aga_df, epi_df):

    dfs = {'S. typhimurium': typh_df, 'S. agalactiae': aga_df, 'S. epidermidis': epi_df}
    
    species_vectors = {}

    for species, df in dfs.items():
        # set index as context_collapsed_group and take frac
        v = df.set_index('context_collapsed_group')['frac']
        species_vectors[species] = v

    # Combine into a single dataframe
    species_matrix = pd.DataFrame(species_vectors).T.fillna(0)


    return species_matrix


def plot_cosine_similarity(species_matrix):

    cos_sim = cosine_similarity(species_matrix)
    cos_sim_df = pd.DataFrame(cos_sim, index=species_matrix.index, columns=species_matrix.index)


    plt.figure()

    sns.heatmap(cos_sim_df, annot=True, cmap='viridis')
    plt.title("Cosine similarity of mutation spectra across species")

    plt.savefig('cosine.png', dpi=300, bbox_inches='tight')


def main():

    typh_single_base, typh_context, aga_single_base, aga_context, epiderm_single_base, epiderm_context = read_data()

    plot_single_base(typh_single_base, aga_single_base, epiderm_single_base)

    plot_context(typh_context, aga_context, epiderm_context)

    species_matrix = create_df_for_cosine(typh_context, aga_context, epiderm_context)

    plot_trinucleotide_grid(typh_context, 'S. typhimurium')
    plot_trinucleotide_grid(aga_context, 'S. agalactiae')
    plot_trinucleotide_grid(epiderm_context, 'S. epidermidis')

    plot_cosine_similarity(species_matrix)


    

    


if __name__ == '__main__':
    main()

    

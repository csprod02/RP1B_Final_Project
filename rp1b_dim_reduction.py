import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import RobustScaler
from sklearn.manifold import TSNE
import umap

from sklearn.metrics.pairwise import cosine_similarity


def create_dataframe_for_PCA():

    species = ['S. typhimurium', 'S. agalactiae', 'S. epidermidis']

    typh_df = pd.read_csv('mutation_freqs_per_sample_typh.csv')
    agalac_df = pd.read_csv('mutation_freqs_per_sample_agalactiae.csv')
    epiderm_df = pd.read_csv('mutation_freqs_per_sample_epidermidis.csv')

    typh_df['species'] = 'S. typhimurium'
    agalac_df['species'] = 'S. agalactiae'
    epiderm_df['species'] = 'S. epidermidis'

    print(typh_df.head())

    combined_df = pd.concat([typh_df, agalac_df, epiderm_df]).drop(columns=['sample'])

    return combined_df


def run_pca(df):

    X = df.drop(columns='species')
    labels = df['species']

    X_logged = np.log1p(df.drop(columns="species"))

    X_scaled = StandardScaler().fit_transform(X_logged)

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    pca_df = pd.DataFrame(
        X_pca,
        columns=["PC1", "PC2"]
        )

    pca_df["species"] = labels.values

    return pca_df


def plot_pca(pca_df):

    plt.figure()

    for species in pca_df["species"].unique():
        subset = pca_df[pca_df["species"] == species]
        plt.scatter(subset["PC1"], subset["PC2"], label=species)

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.title("Mutation spectrum PCA")

    plt.savefig('pca.png', dpi=300, bbox_inches='tight')

    plt.close()


def tsne(df):

    plt.figure()

    X = df.drop(columns='species')
    labels = df['species']

    X_logged = np.log1p(df.drop(columns="species"))
    X_scaled = StandardScaler().fit_transform(X_logged)
    

    X_tsne = TSNE(n_components=2, random_state=1).fit_transform(X_scaled)


    for species in labels.unique():
        idx = labels == species
        plt.scatter(
            X_tsne[idx, 0],
            X_tsne[idx, 1],
            label=species,
            s=10
        )

    plt.legend()
        
    plt.savefig('tsne.png', dpi=300, bbox_inches='tight')

    plt.close()



def run_umap(df):

    plt.figure()

    X = df.drop(columns='species')
    labels = df['species']

    X_logged = np.log1p(X)
    X_scaled = StandardScaler().fit_transform(X_logged)

    reducer = umap.UMAP(random_state=1)
    X_umap = reducer.fit_transform(X_scaled)

    for species in labels.unique():
        idx = labels == species
        plt.scatter(
            X_umap[idx, 0],
            X_umap[idx, 1],
            label=species,
            s=10
        )

    plt.legend(title="Species", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")

    plt.savefig('umap.png', dpi=300, bbox_inches='tight')

    plt.close()
    


def main():
    df = create_dataframe_for_PCA()
    pca = run_pca(df)
    plot_pca(pca)
    tsne(df)

    run_umap(df)





if __name__ == '__main__':
    main()

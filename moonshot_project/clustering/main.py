from pathlib import Path
from plotnine import ggplot, geom_point, aes, geom_smooth
from plotnine_prism import *
import pandas as pd
import math
from sklearn.metrics import silhouette_score


from moonshot_project import DatabaseManager
from moonshot_project import DistanceCalculator
from moonshot_project import Cluster


def minus_log_10(assay_value: float):
    return -math.log10(assay_value * (10 ** (-6)))


def scatter_plot(df: pd.DataFrame, n_values: int, n_clusters: int):
    # Create a scatter plot using the relevant assays to see if structure
    # corresponds to assay
    plot = (ggplot(df, aes(x="r_avg_pIC50", y="f_avg_pIC50"))
            + geom_point(aes(color="factor(Cluster)"))
            + geom_smooth(method="lm", se=False)
            + theme_prism()
            + scale_colour_prism(palette="stained_glass")
            )
    plot.draw()
    plot.save(f"cluster_plots/clusters_{n_values}_{n_clusters}.png")


def create_data_frame(data, ids):
    df = pd.DataFrame(columns=["SMILES", "r_avg_IC50", "f_avg_IC50"])
    for key in data.keys():
        df.loc[key] = data[key]

    df["Cluster"] = ids
    df["r_avg_pIC50"] = df["r_avg_IC50"].apply(minus_log_10)
    df["f_avg_pIC50"] = df["f_avg_IC50"].apply(minus_log_10)
    return df


if __name__ == '__main__':
    all_data_file = Path('covid_submissions_all_info.csv')

    manager = DatabaseManager(database_path='sabs_rdbms.db')
    manager.drop_all()
    manager.create()
    manager.populate_compounds_table(all_data_file=all_data_file)
    manager.populate_assays_table(all_data_file=all_data_file)
    compound_data = manager.select_compounds(2000)

    # The first entry in each compound_data dict is the SMILES string
    smiles_dict = {key: compound_data[key][0]
                   for key in compound_data.keys()}

    # Create an instance of DistanceCalculator and write Tanimoto distances
    # to a csv file
    distance_calculator = DistanceCalculator(smiles_dict)
    distance_calculator.write_to_csv()

    # Create an instance of Cluster and perform clustering
    cluster = Cluster(smiles_dict, distance_file="tanimoto.csv",
                      num_clusters=5)
    cluster_ids, smiles_clusters = cluster.make_clusters(num_iterations=100)
    print(smiles_clusters)

    # Find silhouette score for clustering
    distance_matrix = pd.read_csv("tanimoto.csv")
    distance_matrix = distance_matrix.iloc[:, 1:]
    labels = [cluster_ids[molecule] for molecule in distance_matrix.columns]
    print(silhouette_score(distance_matrix, labels, metric="precomputed"))

    # Create a scatter plot using the relevant assays to see if structure
    # corresponds to assay
    compound_df = create_data_frame(compound_data, cluster_ids)

    scatter_plot(compound_df, 2000, 5)

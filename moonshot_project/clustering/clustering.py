import pandas as pd
import random


class Cluster:

    def __init__(self, smiles_dict: dict[str, str], distance_file: str,
                 num_clusters: int):
        self.smiles_dict = smiles_dict
        self.distance_file = distance_file
        self.num_clusters = num_clusters
        self.df = pd.read_csv(distance_file)
        self.df.index = smiles_dict.keys()

    def _sample_centroids(self, num_samples: int):
        """Function to sample `num_samples` data points from a dataset and return
        a reduced dataset

        :param num_samples: The number of samples to use
        :return: The reduced dataset (centroids)
        """

        country_centroids = random.sample(list(self.smiles_dict.keys()), num_samples)
        return country_centroids

    def _distance(self, mol1: str, mol2: str):
        """Returns distance between two molecules using df

        :param mol1: first molecule
        :param mol2: second molecule
        :return: distance
        """
        return self.df.at[mol1, mol2]

    def _distance_between_rows(self, centroid: str):
        return [self._distance(molecule, centroid)
                for molecule in self.smiles_dict.keys()]

    def _all_distances_between_centroids(self, centroids: list):
        distance_rows = [self._distance_between_rows(centroid) for centroid in centroids]
        distance_frame = pd.DataFrame(distance_rows,
                                      columns=list(self.smiles_dict.keys()))
        return distance_frame.transpose()

    def _assign_points_to_clusters(self, centroids: list):
        distance_frame = self._all_distances_between_centroids(centroids)
        cluster_ids = {}
        for index, datapoint in distance_frame.iterrows():
            minimum = min(datapoint)
            for j in range(distance_frame.shape[1]):
                if minimum == datapoint[j]:
                    cluster_ids[index] = j
                    break
        return cluster_ids

    def _new_centroid(self, current_cluster: list):
        """Here we assign the closest point we have to the centroid

        :param current_cluster: The cluster we are dealing with
        :return: A new centroid
        """
        sum_of_distances = {}
        for molecule in current_cluster:
            sum_of_distance = sum([self._distance(molecule, other)
                                  for other in current_cluster])
            sum_of_distances[molecule] = sum_of_distance
        if sum_of_distances:
            min_value = min(sum_of_distances.values())
            for key, value in sum_of_distances.items():
                if min_value == value:
                    return key
        else:
            return current_cluster[0]

    def _all_new_centroids(self, cluster_ids: dict, old_centroids: list):
        current_clusters = [[] for j in range(len(old_centroids))]
        for index, cluster_id in cluster_ids.items():
            current_clusters[cluster_id].append(index)
        new_centroids = [self._new_centroid(cluster) for cluster in current_clusters]
        return new_centroids

    def fit_model(self):
        starting_centroids = self._sample_centroids(self.num_clusters)
        for i in range(100):
            cluster_ids = self._assign_points_to_clusters(starting_centroids)
            new_centroids = self._all_new_centroids(cluster_ids, starting_centroids)
            starting_centroids = new_centroids

        clusters = {j: [] for j in range(1, self.num_clusters + 1)}
        for molecule in cluster_ids.keys():
            clusters[cluster_ids[molecule] + 1].append(self.smiles_dict[molecule])

        return cluster_ids, clusters

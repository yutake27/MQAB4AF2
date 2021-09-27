import pickle


class Cluster:
    bc40 = '../../../data/raw/bc-40.out'

    def __init__(self):
        cluster_list = []
        with open(self.bc40, 'r') as f:
            for line in f.readlines():
                cluster_list.append(line.split())
            self.origin_cluster_list = cluster_list

    def clustering(self, entry_set: set) -> list:
        extracted_cluster_list = []
        for cluster in self.origin_cluster_list:
            extracted_cluster = []
            for ent in cluster:
                if ent[: 4] in entry_set:
                    extracted_cluster.append(ent)
            if len(extracted_cluster) > 0:
                extracted_cluster_list.append(extracted_cluster)
        return extracted_cluster_list


def main():
    with open('../../../data/interim/extracted_entry.pkl', 'rb') as f:
        entry_list = pickle.load(f)
        entry_set = set(entry_list)
        cluster = Cluster()
        # print(cluster.origin_cluster_list)
        cluster_list = cluster.clustering(entry_set)
        print(cluster_list)
        print(len(cluster_list))


if __name__ == '__main__':
    main()

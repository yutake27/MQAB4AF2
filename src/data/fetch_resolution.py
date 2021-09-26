import requests


class GraphQL:
    api_url = 'https://data.rcsb.org/graphql?query='

    @classmethod
    def get(cls, query):
        res = requests.get(cls.api_url + query)
        return res.json()


def main():
    entry_list = ["4HHB", "12CA", "3PQR"]
    query = "{{entries(entry_ids: [{}]){{rcsb_entry_info{{resolution_combined}}}}}}" \
        .format(','.join(['\"' + ent + '\"' for ent in entry_list]))
    print(query)
    res_json = GraphQL.get(query)
    print(res_json)
    resolution_list = [ent['rcsb_entry_info']['resolution_combined'] for ent in res_json['data']['entries']]
    print(resolution_list)


if __name__ == '__main__':
    main()

import requests


class GraphQL:
    api_url = 'https://data.rcsb.org/graphql?query='

    @classmethod
    def get(cls, query):
        res = requests.get(cls.api_url + query)
        return res.json()

    @classmethod
    def post(cls, query: str, variables: list):
        res = requests.post(cls.api_url, json={'query': query, 'variables': variables})
        print(res.status_code)
        print(res.headers)
        print(res.text)
        assert res.status_code == 200
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

    query = 'query($ids: [String!]!){entries(entry_ids: $ids){rcsb_entry_info{resolution_combined}}}'
    variables = {'ids': entry_list}
    res_json = GraphQL.post(query, variables)
    print(res_json)


if __name__ == '__main__':
    main()

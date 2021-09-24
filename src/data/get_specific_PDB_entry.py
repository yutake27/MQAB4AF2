import json
import requests


class SearchAPI:
    api_url = 'https://search.rcsb.org/rcsbsearch/v1/query?json='

    @classmethod
    def post(cls, query):
        res = requests.get(cls.api_url + json.dumps(query))
        print(res.json())


if __name__ == '__main__':
    with open('./specific_PDB_entry_query.json', 'r') as f:
        query = json.load(f)
        SearchAPI.post(query)

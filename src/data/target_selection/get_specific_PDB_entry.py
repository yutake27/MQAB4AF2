import json
import requests
import pickle


class SearchAPI:
    api_url = 'https://search.rcsb.org/rcsbsearch/v1/query?json='

    @classmethod
    def get(cls, query):
        res = requests.get(cls.api_url + json.dumps(query))
        return res.json()

    @staticmethod
    def parseJson(res_json):
        entry_list = []
        for ent in res_json['result_set']:
            entry_list.append(ent['identifier'])
        return entry_list


def main():
    with open('./specific_PDB_entry_query.json', 'r') as f:
        query = json.load(f)
        res_json = SearchAPI.get(query)
        entry_list = SearchAPI.parseJson(res_json)
        print(entry_list)
        with open('../../../data/interim/extracted_entry.pkl', 'wb') as f:
            pickle.dump(entry_list, f)


if __name__ == '__main__':
    main()

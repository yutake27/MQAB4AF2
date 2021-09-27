import pickle


def main():
    with open('../../../data/interim/extracted_entry.pkl', 'rb') as f:
        entry_list = pickle.load(f)
        entry_set = set(entry_list)
    with open('../../../data/raw/pdb_seqres.txt', 'r') as f:
        output_line = []
        while True:
            header = f.readline()
            seq = f.readline()
            if not header:
                break
            if header[1: 5].upper() in entry_set:
                output_line.append(header)
                output_line.append(seq)
        print(output_line)


if __name__ == '__main__':
    main()

import gzip
import json
import pickle
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Convert buyables pickle file to json for seeding mongodb')
    parser.add_argument('pickle_file')
    parser.add_argument('json_file')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    with open(args.pickle_file, 'rb') as f:
        prices = pickle.load(f)
    prices_json = []
    for smiles, price in prices.items():
        prices_json.append({'smiles': smiles, 'ppg': price})
    prices_json_str = json.dumps(prices_json) + '\n'
    with gzip.open(args.json_file, 'w') as f:
        f.write(prices_json_str.encode('utf-8'))

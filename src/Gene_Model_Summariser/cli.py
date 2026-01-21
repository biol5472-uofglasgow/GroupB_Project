import argparse


def app():
    parser = argparse.ArgumentParser(description="A simple CLI tool.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--help', action='help', help='Show this help message and exit.')
    print('hello world')
    
    return parser.parse_args()
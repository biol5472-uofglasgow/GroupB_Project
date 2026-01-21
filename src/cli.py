import argparse

parser = argparse.ArgumentParser(description="A simple CLI tool.")
parser.add_argument('--help', action='help', help='Show this help message and exit.')

def app():
    print('hello world')
    # Application logic goes here
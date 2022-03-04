import argparse, pathlib

parser = argparse.ArgumentParser()
parser.add_argument('file', type=pathlib.Path, help="mettre un file")
parser.add_argument("-l", "--launch", help="Launch MD", action="store_true")
args = parser.parse_args()

with args.file.open('w') as file:
    file.write("coucou")

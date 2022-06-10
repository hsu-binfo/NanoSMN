import argparse
from ontSMNtools import ontSMNtools

parser = argparse.ArgumentParser(description='Draw SMN CNV plot after graph alignment.', usage='%(prog)s [options]')
parser.add_argument('--path', dest='path', help='path to the location of .gaf file', required=True)
parser.add_argument('--sample', dest='sample', help='name of the sample', required=True)
args = parser.parse_args()
# parser.print_help()

ost = ontSMNtools(args.path)
ost.generate_report(args.sample)
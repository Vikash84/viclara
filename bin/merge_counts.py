#!/usr/bin/env python3
"""
merge mpa-style kraken tables into a single table with clades and relative abundance
"""

# --- standard python imports ---#

import os
import re
import sys
import logging
import argparse
from textwrap import dedent
from itertools import takewhile

# --- third-party imports ---#
import pandas as pd

# --- project specific imports ---#
from utils import mkdir

logging.getLogger().setLevel(logging.DEBUG)
logFormatter = logging.Formatter('[%(asctime)s] - %(''levelname)s - %(message)s',
                                 datefmt='%Y-%m-%d %I:%M:%S')

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logging.getLogger().addHandler(consoleHandler)


def parse_args():
    """
    command line options
    :return:
    """
    parser = argparse.ArgumentParser(prog="merge_mpa_tables.py",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     argument_default=argparse.SUPPRESS, prefix_chars='--',
                                     description=dedent('''Performs a table join on one or more kraken2 mpa-styled 
                                     output files.''')
                                     )

    helpstr = """python merge_mpa_tables.py [options]"""
    required_group = parser.add_argument_group(dedent('''INPUT'''))
    required_group.add_argument("-i", "--input", metavar="<input.txt>", nargs="+", dest="input",
                                help="One or more tab-delimited text tables to join")
    parser.add_argument("-o", "--output", metavar="<output.txt>", dest="output", default=sys.stdout,
                        help="Name of output file in which joined tables are saved")

    return parser


def merge(input, output):
    """
    outputs the table join of the given pre-split string collection.
    :param input: <str> collection of collections of string collections
    :param output: <str> output stream to which matched rows are written.
    :return: output
    """

    merged_tables = pd.DataFrame()
    names = []

    for f in input:
        names.append(os.path.splitext(os.path.basename(f))[0])
        with open(f) as fin:
            headers = [x.strip() for x in takewhile(lambda x: x.startswith('#'), fin)]
            if len(headers) == 0:
                headers = []
                index_col = [0]

        df = pd.read_csv(f,
                         sep='\t',
                         skiprows=1,
                         names=names,
                         index_col=index_col
                         )

        

        if merged_tables.empty:
            merged_tables = df.iloc[:, 0].rename(os.path.splitext(os.path.basename(f))[0]).to_frame()
        else:
            merged_tables = pd.merge(
                df.iloc[:, 0].rename(os.path.splitext(os.path.basename(f))[0]).to_frame(),
                merged_tables,
                how='outer',
                left_index=True,
                right_index=True
            )
    df = merged_tables.reset_index().set_index('index', drop=False)
    df = df.rename(columns={'index': 'sample'})
    df.to_csv(output, index=False, sep="\t")
    return output


def main():
    parser = parse_args()
    args = parser.parse_args()
    logging_level = logging.WARN
    consoleHandler.setLevel(logging_level)

    try:
        mkdir(os.path.dirname(args.output))
    except Exception as e:
        logging.info("redirecting standard output")
    finally:
        merge(input=args.input, output=args.output)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logging.error(str(e))
        logging.debug('', exc_info=True)
        try:
            sys.exit(e.errno)
        except AttributeError:
            sys.exit(1)

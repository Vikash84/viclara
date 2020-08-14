#!/usr/bin/env python

"""
count reads in raw fastq files
"""

import os
import re
import sys
import argparse
from textwrap import dedent
import pandas as pd
from pathlib import Path

from utils import mkdir
import subprocess
from utils import start_time
from utils import end_time
from utils import run_time
from utils import read_type
from utils import find_executable
from utils import run_shell_command

import logging

logfile = os.path.splitext(os.path.basename(os.path.abspath(__file__)))[0] + '.log'
logging.getLogger().setLevel(logging.DEBUG)
logFormatter = logging.Formatter('[%(asctime)s] {%(filename)s:%(lineno)d} %('
                                 'levelname)s - %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p')

fileHandler = logging.FileHandler(filename=logfile)
fileHandler.setFormatter(logFormatter)
logging.getLogger().addHandler(fileHandler)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logging.getLogger().addHandler(consoleHandler)

f_out = open(logfile, 'w')


def count_reads(infile, outfile):
    """
    Perform a read count on FASTQ files in a directory, return dict of counts (int)

    :param infile: <str> path to the input file in fastq format
    :param outfile: <str> file to write the summary
    :return: <str> log file
    """

    # locate the executable
    reformat = find_executable(["reformat.sh"])

    summary_dict = dict()
    # loop through the dict items
    # for sample, reads in read_type(datadir).items():
    #     lines = []
    #     for read in reads:
    #print("\n# sample: {}\n".format(read), sep=' ', end='\n', file=f_out, flush=True)
    call = ["{} in={} -Xmx4G".format(reformat, infile)]
    cmd = " ".join(call)
    try:
        logging.info("[counting reads]\n\t" + cmd +
                        "\n\tBrian Bushnell (2017)."
                        "\n\tBBTools: a suite of fast, multithreaded bioinformatics tools designed for "
                        "analysis of DNA and RNA sequence data. "
                        "\n\thttps://jgi.doe.gov/data-and-tools/bbtools//\n "
                        )
        process = run_shell_command(cmd=cmd, logfile=f_out, raise_errors=True)
        if process:
            with open(logfile, 'r') as infile:
                for line in infile:
                    m1 = re.match(r"^java.*$", line.strip())
                    m2 = re.match(r"^Input:.*$", line.strip())
                    if m1:
                        read = os.path.basename(m1.group(0).split()[-2]).strip("in=")
                        summary_dict[read] = ''
                    if m2:
                        read_count = int(m2.group(0).split()[1])
                        if read in summary_dict:
                            summary_dict[read] = read_count
    except Exception:
        raise Exception("ERROR: COUNTING READS FAILED")

    df = pd.DataFrame(summary_dict.items(), columns=['sample', 'reads'])
    df.to_csv(outfile.name, index=False, sep="\t")
    return outfile


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        argument_default=argparse.SUPPRESS, prefix_chars='--',
        description=dedent('''count reads in fastq file''')
    )
    helpstr = """python count_reads.py [options]"""
    required_group = parser.add_argument_group(dedent('''INPUT'''))

    required_group.add_argument('--infile', type=Path, metavar="<path>",
                                dest='infile',
                                help="path to the file in fastq format"
                                )
    parser.add_argument('--output', type=argparse.FileType('w'), metavar="<file>",
                        dest='output', required=True,
                        help="file to write summary results"
                        )

    args = parser.parse_args()

    log = count_reads(infile=args.infile, outfile=args.output)


if __name__ == '__main__':
    main()

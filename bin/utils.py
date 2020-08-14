#!/usr/bin/env python
"""
utility functions
"""
import argparse
import os
import re
import sys
import errno
import shutil
import tempfile
import fileinput
import contextlib
from pathlib import Path
from itertools import groupby
import pandas as pd

import subprocess
from textwrap import dedent
from datetime import datetime


def start_time():
    """
    determine the start time of a process
    :return:
    """
    # set start time format
    format = "%a %b %d %H:%M:%S %Y"
    stime = datetime.now()
    s = stime.strftime(format)
    print("started:\t{}".format(s))
    return stime


def end_time(stime):
    """

    :param stime:
    :return:
    """

    # set end time format
    format = "%a %b %d %H:%M:%S %Y"
    etime = datetime.now()
    e = etime.strftime(format)
    print("\ncompleted:\t {}".format(e))
    return etime


def run_time(stime, etime):
    """
    determine the total time to run a process
    :param stime:
    :param etime:
    :return:
    """

    tdelta = etime - stime

    # format the time delta object to human readable form
    d = dict(days=tdelta.days)
    d['hrs'], rem = divmod(tdelta.seconds, 3600)
    d['min'], d['sec'] = divmod(rem, 60)

    if d['min'] == 0:
        fmt = '{sec} sec'
    elif d['hrs'] == 0:
        fmt = '{min} min {sec} sec'
    elif d['days'] == 0:
        fmt = '{hrs} hr(s) {min} min {sec} sec'
    else:
        fmt = '{days} day(s) {hrs} hr(s) {min} min {sec} sec'
    print("\n[ALL done] Runtime: " + '\t' + fmt.format(**d))


def find_executable(names, default=None):
    """
    find an executable PATH from the given list of names.
    Raises an error if no executable is found.  Provide a
    value for the "default" parameter to instead return a value.

    :param names: list of given executable names
    :param default:
    :return: <str> path to the first executable found in PATH from the given list of names.
    """
    exe = next(filter(shutil.which, names), default)

    if exe is None:
        print("Unable to find any of {} in PATH={}".format(names, os.environ["PATH"]))
        print("\nHint: You can install the missing program using conda or homebrew or apt-get.\n")
        raise Exception
    return exe


def run_shell_command(cmd, logfile, raise_errors=False, extra_env=None):
    """
    run the given command string via Bash with error checking

    :param cmd: command given to the bash shell for executing
    :param logfile: file object to write the standard errors
    :param raise_errors: bool to raise error if running command fails/succeeds
    :param extra_env: mapping that provides keys and values which are overlayed onto the default subprocess environment.
    :return:
    """

    env = os.environ.copy()

    if extra_env:
        env.update(extra_env)

    try:

        subprocess.check_call("set -euo pipefail; " + cmd,
                                  shell=True,
                                  stderr=logfile,
                                  universal_newlines=True,
                                  executable="/bin/bash"
                                  )
    except (subprocess.CalledProcessError, OSError) as error:
        # except subprocess.CalledProcessError as error:
        rc = error.returncode
        if rc == 127:
            extra = "Are you sure this program is installed?"
        else:
            extra = " "
        print("Error occurred: shell exited with return code: {}\ncommand running: {}\n{}".format(
            error.returncode, cmd, extra), file=sys.stderr
        )
        if raise_errors:
            raise
        else:
            return False

    except FileNotFoundError as error:
        print("Unable to run shell command using {}! tool requires {} to be installed.".format(
            error.filename, error.filename), file=sys.stderr
        )
        if raise_errors:
            raise
        else:
            return False
    else:
        return True


# def print_error(message, **kwargs):
#     """
#     Formats *message* with *kwargs* using :meth:`str.format` and
#     :func:`textwrap.dedent` and uses it to print an error message to
#     ``sys.stderr``.
#     """
#     print(**kwargs)
#     print("\nERROR: " + dedent(message.format(**kwargs)).lstrip("\n") + "\n", file=sys.stderr)


# def run_shell_command(cmd, logfile, raise_errors=False, extra_env=None):
#     """
#     run the given command string via Bash with error checking
#
#     :param cmd: command given to the bash shell for executing
#     :param logfile: file object to write the standard errors
#     :param raise_errors: bool to raise error if running command fails/succeeds
#     :param extra_env: mapping that provides keys and values which are overlayed onto the default subprocess environment.
#     :return:
#     """
#
#     env = os.environ.copy()
#
#     if extra_env:
#         env.update(extra_env)
#
#     try:
#         # Use check_call() instead of run() since the latter was added only in Python 3.5.
#         p = subprocess.check_call(
#             "set -euo pipefail; " + cmd,
#             stderr=logfile,
#             shell=True,
#             universal_newlines=True,
#             executable="/bin/bash",
#             env=env
#         )
#
#     except subprocess.CalledProcessError as error:
#         print_error(
#             "shell exited {rc} when running: {cmd}{extra}",
#             rc=error.returncode,
#             cmd=cmd,
#             extra="\nAre you sure this program is installed?" if error.returncode == 127 else "",
#         )
#         if raise_errors:
#             raise
#         else:
#             return False
#
#     except FileNotFoundError as error:
#         print_error(
#             """
#             Unable to run shell commands using {shell}!
#             tool requires {shell} to be installed.
#             """,
#             shell=error.filename
#         )
#         if raise_errors:
#             raise
#         else:
#             return False
#     else:
#         return True
#
#
# def print_error(message, **kwargs):
#     """
#     Formats *message* with *kwargs* using :meth:`str.format` and
#     :func:`textwrap.dedent` and uses it to print an error message to
#     ``sys.stderr``.
#     """
#     print(**kwargs)
#     print("\nERROR: " + dedent(message.format(**kwargs)).lstrip("\n") + "\n", file=sys.stderr)


def first_line(text):
    """
    Returns the first line of the given text, ignoring leading and trailing
    whitespace.
    """
    return text.strip().splitlines()[0]


def mkdir(directory):
    """
    recursivley create a directory if it does not exist

    :param directory: path to the directory to be created
    :return: directory
    """
    directory = os.path.abspath(directory)
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise  # raises the error again
    return directory


def available_cpu_cores(fallback: int = 1) -> int:
    """
    Returns the number (an int) of CPU cores available to this **process**, if
    determinable, otherwise the number of CPU cores available to the
    **computer**, if determinable, otherwise the *fallback* number (which
    defaults to 1).
    """
    try:
        # Note that this is the correct function to use, not os.cpu_count(), as
        # described in the latter's documentation.
        #
        # The reason, which the documentation does not detail, is that
        # processes may be pinned or restricted to certain CPUs by setting
        # their "affinity".  This is not typical except in high-performance
        # computing environments, but if it is done, then a computer with say
        # 24 total cores may only allow our process to use 12.  If we tried to
        # naively use all 24, we'd end up with two threads across the 12 cores.
        # This would degrade performance rather than improve it!
        return len(os.sched_getaffinity(0))
    except:
        # cpu_count() returns None if the value is indeterminable.
        return os.cpu_count() or fallback


def locate_picard(tool, java_file, versions):
    """
    get the full path to the picard java executable jar file (picard.jar)

    :param tool:
    :param java_file:
    :param versions: list of version numbers to form a regular expression to search
    :return:
    """
    # locate the executable
    picard = find_executable([tool])

    # get the path(s) to the executable java file (picard.jar)
    jar_file = None
    base_paths = os.environ['PATH'].split(':')
    for p in base_paths:
        match = re.search(picard, p)
        if match is not None:
            for root, dirs, fns in os.walk(p):
                for fn in fns:
                    if fn in ['picard.jar']:
                        jar_file = os.path.join(root, fn)

    jar = None
    # versions_regex = '(?:%s)' % '|'.join(versions)

    # base_path = Path(os.path.dirname(os.environ['PATH'].split(':')[1]))
    p = os.path.dirname(os.environ['PATH']).split(':')
    for f in p:
        match = re.search('picard', f)
        if match is not None:
            for root, dirs, fns in os.walk(p):
                for fn in fns:
                    if fn in ['picard.jar']:
                        jar = os.path.join(root, fn)
    print(jar)

    # for root, dirs, fnames in os.walk(base_path):
    #     for filename in fnames:
    #         if filename == java_file:
    #             jar_file = os.path.join(root, filename)
    #             if os.path.isfile(jar_file):
    #                 _version = os.path.basename(os.path.dirname(os.path.normpath(jar_file)))
    #                 match = re.search(versions_regex, _version)
    #                 if match is not None:
    #                     jar = os.path.join(root, filename)
    return jar_file


def nthreads(value):
    """
    Argument value validation and casting function for --nthreads.
    """

    if value.lower() == 'auto':
        return available_cpu_cores()

    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("{} is not an integer or the word 'auto'".format(value)) from None


def read_type(datadir):
    """
    determine reads type (either as single-end or paired-end reads)
    this assumes that files are named according to illumina file naming
    convention: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm


    :param datadir: <str> path to the directory having the raw reads
    :return: reads_dict: dict
        dict object having sample id and paths to the read files
        {'B66': ['B66_S9_L001_R1_001.fastq.gz', 'B66_S9_L001_R2_001.fastq.gz']}
    """

    ext = "fastq fq".split()
    samples_dict = dict()

    if not os.path.exists(datadir):
        print("data directory does not exist")
    # elif os.path.exists(datadir) and len(os.listdir(datadir)) == 0:
    #     print("data directory is empty")
    else:
        # traverse the directory to look for patterns: _R1, _R2, _1, _2
        for dirname, dirnames, filenames in os.walk(os.path.abspath(datadir)):
            filenames.sort()
            for filename in filenames:
                if ext[0] in filename or ext[1] in filename:
                    sampleid = filename.split('.fastq')[0].split('.fq')[0]
                    sampleid = sampleid.replace("_R1", "").replace("_r1", "").replace("_R2", "").replace("_r2", "")
                    sampleid = sampleid.replace('.R1', '_R1').replace('.R2', '_R2').replace("_R1", "").replace("_R2",
                                                                                                               "")
                    sampleid = re.sub(r'_S(\d+)_L\d{3}_\d{3}', '', sampleid)
                    sampleid = sampleid.replace(".", "_").replace(" ", "_")
                    sampleid = re.sub(r'\_.*', '', sampleid)

                    if os.path.basename(dirname) == 'contaminants':
                        continue
                    else:
                        filepath = os.path.join(dirname, filename)
                        if not sampleid in samples_dict:
                            samples_dict[sampleid] = [filepath]
                        else:
                            samples_dict[sampleid].append(filepath)
        # print("\n{:<30}\t{:>20}".format("Sample", "Files"))
        # for sample, reads in samples_dict.items():
        #     print("{:<30}\t{:>20}\n".format(sample, '\n\t\t\t\t'.join(reads)))
    return samples_dict


def compress_fq(filename, suffix=".gz"):
    """
    compress fastq format files to gz

    :param filename: file path
    :param suffix: suffix for the file extension
    :return: filename: file path
    """

    match = re.findall(r"(\w+.)", os.path.basename(filename))

    if match[-1] == 'fastq' or match[-1] == 'fq':
        outfile = filename + suffix
        call = ['gzip --suffix={} \\\n\t {}\n'.format(suffix, filename)]
        cmd = ' '.join(call)
        print("\ncompressing via: {}".format(cmd))
        try:
            subprocess.check_call(cmd, shell=True)
        except:
            print("ERROR: FILE COMPRESSION FAILED")
        return outfile


def uncompress_fq(filename, suffix=".fastq"):
    """
    compress fastq format files to gz

    :param filename: file path
    :param suffix: suffix for the file extension
    :return: filename: file path
    """

    match = re.findall(r"(\w+.)", os.path.basename(filename))
    if match[-2] == 'fastq.' or match[-2] == 'fq.' and match[-1] == 'gz':
        outfile = os.path.join(os.path.dirname(filename), ''.join(match[:-1]).strip('.'))
        if Path(outfile).is_file():
            print("[decompression done], check file: \n\t {}".format(outfile))
        else:
            call = ['gunzip --keep \\\n\t --force \\\n\t --suffix={} \\\n\t {}\n'.format(suffix, filename)]
            cmd = ' '.join(call)
            print("\nuncompressing via: {}".format(cmd))
            try:
                subprocess.check_call(cmd, shell==True)
            except:
                print("ERROR: FILE UNCOMPRESSION FAILED")
        return outfile


def uncompress_fasta(filename, suffix=".fasta"):
    """

    :param filename:
    :param suffix:
    :return:
    """

    match = re.findall(r"(\w+.)", os.path.basename(filename))
    if match[-1] == 'gz':
        out = os.path.splitext(filename)[0]
        print(out)
        ext = os.path.splitext(out)[1]
        outfile = os.path.join(os.path.dirname(filename), ''.join(match[:-1]).strip(ext)) + suffix

        if Path(outfile).is_file():
            print("[decompression done], check file: \n\t {}".format(outfile))
        else:
            cmd = 'gzip -dc <{}> {}'.format(filename, outfile)
            try:
                print("[uncompressing file] {}\n\t".format(cmd))
                p = subprocess.check_call(cmd, shell=True)
                if p == 0:
                    os.remove(os.path.abspath(filename))
            except Exception as error:
                print("Error: {}".format(error))
        return outfile


def _out(read, prefix):
    """
    get output names for bbmap executables output files

    :param read: basename of the read file
    :param prefix: string to include in the output name
    :return: renamed basename read
    """
    print(os.path.splitext(read)[1])
    M = re.match(r"(\w+.*?)_R(\d).fastq(?:.gz|)$", read)
    out = M.group(1) + prefix + M.group(2) + '.fastq.gz'
    return out


def str2bool(v):
    """
    convert string to a boolean

    :param v: string value either true or false
    :return: False or True
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ['true', 't']:
        return True
    elif v.lower() in ['false', 'f']:
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


@contextlib.contextmanager
def tempdir(cleanup=True):
    """
    create a temporary directory and cleanup

    :param cleanup: boolean, if True, the temp directory will be cleaned
    :return:
    """
    dirpath = tempfile.mkdtemp()
    try:
        yield dirpath
    finally:
        if cleanup:
            shutil.rmtree(dirpath)
        else:
            return dirpath


def copy_file(src, dest):
    """
    copy file from one directory to another

    :param src: path of the source file
    :param dest: path of the destination
    :return:
    """
    # copy to directory
    try:
        shutil.copy(src, dest)

    # source and destination are same
    except shutil.SameFileError:
        print("source\t{} and destination\t{} represents the same file.".format(src, dest))

    # permission error
    except PermissionError:
        print("Permission denied!")

    # any other errors
    except Exception as error:
        print(error)
    return dest


def check_contigs(reference):
    """
    traverse through the samples assembly directories to find contigs
    {megahit: '.contigs.fa', metaspades: 'contigs.fasta', ray-meta: Contigs.fa}

    returns: refs_dict <dict> a dictionary with sampleids and reference files and their basenames

    :param reference: <str> path to the assembly output directory (from megahit, metaspades or ray-meta)
    :return refs_dict: <dict> a dictionary with sampleid as key and reference file and the reference index basename as
                                values
    """

    is_contig = False
    refs_dict = dict()
    if Path(reference).is_file():
        print("please specify the path to the assembly output containing sample subdirectories")
    else:
        while not is_contig:
            for dirname, dirnames, filenames in os.walk(reference):
                filenames.sort()
                for f in filenames:
                    match = re.match(r'^(^(?!.*k(\d)).*\.contigs.fa$|^contig.fa$|contigs.fasta$|^Contigs.fasta$)', f)
                    if match is None:
                        pass
                    else:
                        reference = os.path.join(dirname, match.group(0))
                        # create an index basename (same as the fasta name)
                        index_base = os.path.splitext(reference)[0]
                        sample = os.path.basename(dirname)
                        if not sample in refs_dict:
                            refs_dict[sample] = [reference, index_base]
                    is_contig = True
    return refs_dict


def fasta_iterator(fasta_name):
    """
    Given a fasta file. yield tuples of header, sequence
    Author: Brent Pedersen
    https://www.biostars.org/p/710/
    """
    fa_fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fa_fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq
    fa_fh.close()


def fix_ftp_paths(assembly_text_file):
    """

    :param assembly_text_file:
    :return:
    """
    bak = os.path.join(os.path.dirname(assembly_text_file), os.path.basename(assembly_text_file))
    os.chmod(assembly_text_file, 0o777)
    df = pd.read_csv(assembly_text_file, sep="\t", skiprows=1)
    df = df[df.ftp_path.str.startswith('ftp')]
    os.remove(assembly_text_file)
    df.to_csv(bak, sep="\t", index=False)
    return bak


def ref_strings():
    """
    references to the tools used
    :return:
    """

    tools_refs = {
        'fastqc': "\nAndrews, S: FastQC - A high throughput sequence QC analysis tool "
                  "\nhttps://www.bioinformatics.babraham.ac.uk/projects/fastqc/\n",

        'bbmap': "\nBrian Bushnell (2017)."
                 "\nBBTools: a suite of fast, multithreaded bioinformatics tools designed for "
                 "analysis of DNA and RNA sequence data. "
                 "\nhttps://jgi.doe.gov/data-and-tools/bbtools/\n",

        'trimmomatic': "\nBolger, A. M., Lohse, M., & Usadel, B. (2014)."
                       "\nTrimmomatic: A flexible trimmer for Illumina Sequence Data. "
                       "\nBioinformatics, Volume 30, Issue 15, 1 August 2014, Pages 2114–2120 "
                       "\nhttps://doi.org/10.1093/bioinformatics/btu170\n",

        'bowtie2': "\nLangmead, B., & Salzberg, S. L. (2012). "
                   "\nFast gapped-read alignment with Bowtie 2. "
                   "\nNature methods, 9(4), 357–359. https://doi.org/10.1038/nmeth.1923\n",

        'bwa': "\nLi, H., Durbin, D. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform, "
               "\nBioinformatics, Volume 25, Issue 14, 15 July 2009, Pages 1754–1760, "
               "\nhttps://doi.org/10.1093/bioinformatics/btp324\n",

        'samtools': "\nLi, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., "
                    "\nDurbin, R., & 1000 Genome Project Data Processing Subgroup (2009). "
                    "\nThe Sequence Alignment/Map format and SAMtools. "
                    "\nBioinformatics (Oxford, England), 25(16), 2078–2079. "
                    "\nhttps://doi.org/10.1093/bioinformatics/btp352\n",

        'kraken2': "\nWood, D.E., Lu, J. & Langmead, B. (2019) "
                   "\nImproved metagenomic analysis with Kraken 2. "
                   "\nGenome Biol 20, 257 "
                   "\nhttps://doi.org/10.1186/s13059-019-1891-0\n",

        'centrifuge': "\nKim D, Song L, Breitwieser FP, and Salzberg SL. (2016). "
                      "\nCentrifuge: rapid and sensitive classification of metagenomic sequences. "
                      "\nGenome Research 2016 "
                      "\nhttps://genome.cshlp.org/content/early/2016/11/16/gr.210641.116\n",

        'last': "\nKiełbasa SM, Wan R, Sato K, Horton P, Frith MC. (2011) "
                "\nAdaptive seeds tame genomic sequence comparison.  "
                "\nGenome Res. 21(3):487-93.\n",

        'picard': "\nBroad Institute "
                  "\nA set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and "
                  "\nformats such as SAM/BAM/CRAM and VCF. "
                  "\nhttps://broadinstitute.github.io/picard/\n",

        'trinity': "\nGrabherr, M. G., Haas, B. J., Yassour, M., Levin, J. Z., Thompson, D. A., Amit, I., Adiconis, "
                   "X., "
                   "\nFan, L., Raychowdhury, R., Zeng, Q., Chen, Z., Mauceli, E., Hacohen, N., Gnirke, A., Rhind, N., "
                   "\ndi Palma, F., Birren, B. W., Nusbaum, C., Lindblad-Toh, K., Friedman, N., … Regev, A. (2011). "
                   "\nFull-length transcriptome assembly from RNA-Seq data without a reference genome. "
                   "\nNature biotechnology, 29(7), 644–652. "
                   "\nhttps://doi.org/10.1038/nbt.1883\n",
        'mummer': "\nKurtz, S., Phillippy, A., Delcher, A.L. et al. (2004) "
                  "\nVersatile and open software for comparing large genomes. "
                  "\nGenome Biol 5, R12 "
                  "\nhttps://doi.org/10.1186/gb-2004-5-2-r12\n",

        'gatk': "\nMcKenna A, Hanna M, Banks E, et al. "
                "\nThe Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing "
                "data."
                "\nGenome Res. 2010;20(9):1297‐1303. "
                "\ndoi:10.1101/gr.107524.110 ",

        'bedtools': "\nQuinlan Aaron, Hall Ira. (2010) "
                    "\nBEDTools: a flexible suite of utilities for comparing genomic features."
                    "\nBioinformatics, Volume 26, Issue 6, 15, Pages 841–842, "
                    "\nhttps://doi.org/10.1093/bioinformatics/btq033",

        'ncbi-genome-download': "\nBli, Kai "
                                "\nDownload genome files from the NCBI FTP server. "
                                "\nhttps://pypi.org/project/ncbi-genome-download/"

    }

    return tools_refs

# if __name__ == '__main__':
#     main()

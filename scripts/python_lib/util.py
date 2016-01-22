__author__ = 'Jacques van Helden'
__license__ = "GPL"
__version__ = "0.01"
__maintainer__ = "Jacques van Helden"
__email__ = "Jacques.van-Helden@univ-amu.fr"
__status__ = "In development"

import pandas as pd

def read_sample_ids(file_path:str, base_dir:str=".", column:int=1, verbosity:int=0) -> list:
    """Read sample descriptions from a tab-delimited file. 

    The first column of the sample description file contains the
    sample ID. The other columns are currently ignored.

    Note: this function has been renamed read_column_from_tab, because
    it is generic. The function read_sample_id() is maintained for
    backward compatibility.

    """
    if verbosity >= 4:
        print ("read_sample_ids()\t" +"column\t" + str(column) +"\tfile_path\t" + file_path)
    sample_ids = read_column_from_tab(file_path, column, verbosity)
    return(sample_ids)


def read_column_from_tab(file_path:str, column:int=1, verbosity:int=0) -> list:
    """Read a column from a tab-delimited file, and returns the values in
    a list.

    Lines starting with a semicolumn (;) are considered as comment
    lines, and thus ignored.

    A line starting with a # can optionally be used at the beginning
    of the file to specify column headers (still not taken in
    consideration, but should be treated in future versions of this
    function).

    :param file_path: path to a tab-delimited text file
    :type file_path: string
    :param column: number of the column to read (Default: 1)
    :type column: Integer >= 1
    :param verbosity: verbosity level
    :type verbosity: Integer
    :return: the list of values from the specified columns.
    :rtype: list of strings
    """
    if verbosity >= 4:
        print ("read_column_from_tab()\t" +"column\t" + str(column) +"\tfile_path\t" + file_path)
    values = []

    f = open(file_path, "r")
    for line in f.readlines():
        line = line.rstrip("\n") ## Suppress carriage.return
        if line[0:1] != ';': ## Skip comment lines (starting by "--")
            if line[0:1] == '#': ## Line containing column headers
                line = line.lstrip("#")
                keys = line.split("\t")
#                print("\t".join(keys))
            else:
                fields = line.split("\t")
                current_value = fields[column -1]
#                current_value = fields[0]
                values.append(current_value)
    if verbosity >= 5:
        print("\t".join(["read_sample_ids()", "Result:", str(len(values)), "sample IDs"]))
    return(values)
    
def read_table(file:str, verbosity:int=0, header:int=0, skip_blank_lines=True, comment=';') -> pd.DataFrame:
    """Read a tab-delimited text file and return the content as a data
    frame (object of the class panda.DataFrame).

    This is a simple wrapper around panda.read_csv, with the appropriate
    options for our customized tab-separated files (tsv): comment lines
    start with ';', the first non-comment line contains the
    header. These options can be overwritten in the function call.

    """
    if verbosity >= 3:
        print ("read_table()\t" + file)
    df = pd.read_csv(file, sep="\t", 
                     header=header, 
                     skip_blank_lines=skip_blank_lines,
                     comment=comment)
    if verbosity >= 4:
        print("\tColumns:\t" + ";".join(list(df.columns)))
    return(df)
    
# ################################################################
# def read_sample_ids(file, base_dir=".", column=1, verbosity=0):
#     """Read sample descriptions from a tab-delimited file. 
#
#     The first column of the sample description file contains the
#     sample ID. The other columns are currently ignored.
#
#     Lines starting with a semicolumn (;) are considered as comment
#     lines, and thus ignored.
#
#     A line starting with a # can optionally be used at the beginning
#     of the file to specify column headers (still not taken in
#     consideration, but should be treated in future versions of this
#     function).
#
#     :param file: path to the sample description file
#     :type file: string
#     :param column: number of the column containing the sample ID (Default: 1)
#     :type column: Integer >= 1
#     :param verbosity: verbosity level
#     :type verbosity: Integer
#     :return: a list of sample IDs, taken from the first column of the sample file.
#     :rtype: list of strings
#
#     """
#     if verbosity >= 1:
#         print ("read_sample_ids()\t" +"Reading sample IDs from file\t" + file)
#     samples = []
#
#     f = open(file, "r")
#     for line in f.readlines():
#         line = line.rstrip("\n") ## Suppress carriage.return
#         if line[0:1] != ';': ## Skip comment lines (starting by "--")
#             if line[0:1] == '#': ## Line containing column headers
#                 line = line.lstrip("#")
#                 keys = line.split("\t")
# #                print("\t".join(keys))
#             else:
#                 fields = line.split("\t")
#                 samples.append(fields[column-1])
# #                print(fields[column -1])
#     if verbosity >= 1:
#         print("\t".join(["read_sample_ids()", "Result:", str(len(samples)), "sample IDs"]))
#     return(samples)
    
################################################################
def read_chipseq_design(file, test_column=1, input_column=2, verbosity=0):
    """Read ChIP-seq analysis design from a tab-delimited file. 

    Each row describes one analysis, which requires to define two
    samples: (1) test sample (generally, the ChIP result); (2) input
    sample (can be of different types: mock, genomic input, a ChIP
    result in control conditions, or any other relevant input).

    Lines starting with a semicolumn (;) are considered as comment
    lines, and thus ignored.

    A line starting with a # can optionally be used at the beginning
    of the file to specify column headers (still not taken in
    consideration, but should be treated in future versions of this
    function).

    :param file: path to the sample description file
    :type file: string
    :param test_column: number of the column containing the test ID (Default: 1)
    :type test_column: Integer >= 1
    :param input_column: number of the column containing the input ID (Default: 1)
    :type input_column: Integer >= 1
    :param verbosity: verbosity level
    :type verbosity: Integer
    :return: two lists giving respectively test and input sample IDs.
    :rtype: lists of strings

    """
    if verbosity >= 4:
        print ("read_sample_ids()\t" +"Reading ChIP-seq design from file\t" + file)
    test_ids = []
    input_ids = []

    f = open(file, "r")
    for line in f.readlines():
        line = line.rstrip("\n") ## Suppress carriage.return
        if line[0:1] != ';': ## Skip comment lines (starting by "--")
            if line[0:1] == '#': ## Line containing column headers
                line = line.lstrip("#")
                keys = line.split("\t")
#                print("\t".join(keys))
            else:
                fields = line.split("\t")
                test = fields[test_column-1]
                test_ids.append(test)
                input = fields[input_column-1]
                input_ids.append(input)
                if verbosity >= 5:
                    print("\t".join(["\tanalysis", test, "versus", input]))
    if verbosity >= 5:
        print("\t".join(["read_chipseq_design()", "Result:", str(len(test_ids)), "analyses"]))
    return(test_ids, input_ids)
    

################################################################
def glob_multi_dir(dir_list, pattern="", base_dir=".", ext="", verbosity=0):
    """Given a list of directories, returns the list of files matching a
    given pattern. In addition to the list of matching files, this
    function returns the list of directories associated to each file,
    and a list of basenames. An extension can optionally be removed
    from the basename list.

    :param dir_list: list of directories
    :type dir_list: list of strings
    :param pattern: matching pattern for file names (default: ""). If not specified, all files found in the directories are returned.
    :type pattern: string
    :param base_dir: Base directory (default: "."). All directories of dir_list are computed relative to base_dir.
    :type base_dir: string
    :param ext: File extension to suppress in basename list (default: "").
    :type ext: string
    :return: thre lists of the same length, indicating the file paths, associated directories, and basenames, respectively.
    :rtype: list of strings

    """
    files = [] ## List of paths to the files
    file_dirs = [] ## List of directories associated to each file (same length as the "files" list)
    basenames = [] ## List of file basenmes (same length as the "files" list
    if verbosity >= 4:
        print("\t".join(["glob_multi_dir()", "Listing files in", str(len(dir_list)), "directories", "pattern: " + pattern]))
    import glob
    import os
    for current_dir in dir_list:
        ## Find files matching the pattern in the current directory
        glob_pattern = base_dir + "/" + current_dir + "/" + pattern
        files_in_dir = glob.glob(glob_pattern)
        files = files + files_in_dir

        ## Handle the list of attributes of the matched files
        for filename in files_in_dir:
            ## Append current directory to the list of file-associated directories
            file_dirs.append(current_dir)

            basename = os.path.basename(filename)
            if ext != "":
                basename = basename.replace(ext, "")
            basenames.append(basename)
#            basenames = basenames + [basename]
#        print ("files_in_dir:\t"+"; ".join(files_in_dir))
    if verbosity >= 5:
        print("\t".join(["glob_multi_dir()", "Result:", str(len(files)), "files"]))
    return(files, file_dirs, basenames)


################################################################
def report_numbered_list(list):
    """
    Taking as input a list of strings, return a numbered list in
    reStructured text for the snakemake report.

    :param list: List of items to be numbered.
    :type list: A list of strings.
    :return: A numbered list in reStructured format
    :rtype: String

    """
    result = ""
    n = 0
    for element in list:
        n = n + 1
        result = result + str(n) + ". " + element + "\n"
    return(result)

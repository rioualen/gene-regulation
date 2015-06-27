__author__ = 'Jacques van Helden'
__license__ = "GPL"
__version__ = "0.01"
__maintainer__ = "Jacques van Helden"
__email__ = "Jacques.van-Helden@univ-amu.fr"
__status__ = "Embryonic"

def read_sample_ids(file):
    """Read sample descriptions from a tab-delimited file. 

    The first column of the sample description file contains the
    sample ID. The other columns are currently ignored.

    Lines starting with a semicolumn (;) are considered as comment
    lines, and thus ignored.

    A line starting with a # can optionally be used at the beginning
    of the file to specify column headers (still not taken in
    consideration, but should be treated in future versions of this
    function).

    :param file: path to the sample description file
    :type file: string
    :return: a list of sample IDs, taken from the first column of the sample file.
    :rtype: list of strings

    """
    print ("read_sample_ids()\t" +"Reading sample IDs from file\t" + file)
    samples = []
    f = open(config["dir"]["base"] + "/" + file, "r")
    for line in f.readlines():
        line = line.rstrip("\n") ## Suppress carriage.return
        if line[0:1] != ';': ## Skip comment lines (starting by "--")
            if line[0:1] == '#': ## Line containing column headers
                line = line.lstrip("#")
                keys = line.split("\t")
#                print("\t".join(keys))
            else:
                fields = line.split("\t")
                samples.append(fields[0])
#                print(fields[0])
    print("\t".join(["read_sample_ids()", "Result:", str(len(samples)), "sample IDs"]))
    return(samples)
    

def glob_multi_dir(dir_list, pattern="", base_dir="."):
    """Given a list of directories, returns the list of files matching a
    given pattern.

    :param dir_list: list of directories
    :type dir_list: list of strings
    :param pattern: matching pattern for file names (default: ""). If not specified, all files found in the directories are returned.
    :type pattern: string
    :param base_dir (default: "."). Base directory. All directories of dir_list are computed relative to base_dir.
    :type base_dir: string
    :return: a list of file paths, defined relative to base_dir.
    :rtype: list of strings

    """
    files = []
    print("\t".join(["glob_multi_dir()", "Listing files in", str(len(dir_list)), "directories", "pattern: " + pattern]))
    import glob
    for dir in dir_list:
        glob_pattern = base_dir + "/" + dir + "/" + pattern
#        print ("glob pattern:\t"+glob_pattern)
        dir_files = glob.glob(glob_pattern)
#        print ("dir_files:\t"+"; ".join(dir_files))
        files = files + dir_files
    print("\t".join(["glob_multi_dir()", "Result:", str(len(files)), "files"]))
    return(files)


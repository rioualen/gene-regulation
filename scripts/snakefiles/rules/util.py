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
    

def glob_multi_dir(dir_list, pattern="", base_dir=".", ext=""):
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
    print("\t".join(["glob_multi_dir()", "Listing files in", str(len(dir_list)), "directories", "pattern: " + pattern]))
    import glob
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
    print("\t".join(["glob_multi_dir()", "Result:", str(len(files)), "files"]))
    return(files, file_dirs, basenames)


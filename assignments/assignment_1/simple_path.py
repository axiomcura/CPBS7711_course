import os
import sys
import glob
import argparse
from datetime import datetime
from collections import defaultdict
from itertools import permutations

#--------------------
# Classes
#--------------------
class StringDB:
    """ Creates a simple python object that handles lookups into the
    STRING.txt file.

    Arguments
    ---------
    fname : str
        path to database

    Returns
    -------
    StringDB Object
        Small container where the the String.txt is stored
    """

    def __init__(self, fname):
        db = self._parse_string_db(fname)
        self._db = db
        self.size = "{} MB".format(round(sys.getsizeof(db)/1024**2, 4))


    def _parse_string_db(self, fname: str) -> dict:
        """ Parses the STRING.txt file and converts it into a dictionary

        Arguments
        fname : str
            path to STRING.txt file

        Returns
        -------
        dict
            dictionary containing protein pairs as key and scores as its value

        # Example:
        >>> results = {"PROT1 PROT2" : score, ...}
        """

        # Reads the STRING.txt file and iterates over all lines
        # each line is split into two parts
        # -= protein pairs (key)
        # -- score (value)
        # then stored into dictionary
        db_contents = defaultdict(lambda: None)
        with open(fname, "r") as infile:
            record_contents = infile.readlines()
            for records in record_contents:
                data = records.replace("\n", "").split("\t")
                protein_pair = "{} {}".format(data[0], data[1])
                score = data[2]
                db_contents[protein_pair] = score
        return db_contents


    def get_pair_score(self, locus_name: str, proteins: list, interaction_type="pp") -> list:
        """ Accepts a list of genes and queries to the database

        Summary:
        -------
        Creates all possible permutations of protein interactions within the protein list.
        Each pair will be queried into database and returns a score. The database will return
        'None' if the score is not found and will not be recorded into the results. In addition,
        a "tracking" list is also implemented to prevent repetitive query. This means that
        reversed interaction queries will be ignored if the original query has been recorded.

        Argument
        -------
        genes : list
            list of genes found in the locus

        interaction_type : str (choices=["pp", "pd", "pr", "rc", "cr", "gl", "pm", "mp"])
            Interaction type of both genes/molecules. Supported interaction types are:
            - p:protein - protein interaction
            - pd: protein -> DNA
            - pr: protein -> reaction
            - rc: reaction -> compound
            - cr: compound -> reaction
            - gl: genetic lethal relationship
            - pm: protein-metabolite interaction
            - mp: metabolite-protein interaction

        More information about compatibility found here:
        http://manual.cytoscape.org/en/stable/Supported_Network_File_Formats.html#sif-format

        Returns
        -------
        list
            Contains a list of strings that describes the interaction type
            between two genes and its score. This cotnents is what is going
            to be used to produce the sif file

        """

        # type checking
        if not isinstance(proteins, list):
            genes = [proteins]

        # getting all possible combinations
        # -- query all combinations into the database
        # -- reversed queries are ignored if original query is recorded
        pairs = permutations(proteins, 2)

        results = []
        searched = []
        for gene1, gene2 in pairs:
            query = "{} {}".format(gene1, gene2)
            query_rev = "{} {}".format(gene2, gene1)
            score = self._db[query]
            if score == None:
                continue
            if  query_rev in searched:
                continue
            result = "{} {} {} {}".format(gene1, interaction_type, gene2, score)
            results.append(result)
            searched.append(query)
        return results


#--------------------
# Functions
#--------------------
def parse_input(input_file : str) -> dict:
    """ Parses input file and converts it into a dictionary

    Arguments
    ---------
    input_file : str
        path to input file

   Returns
   -------
   dict
        dictionary containing the locus name as the key and all
        the genes within the locus as value

   Example
   --------
    >>> # this is an example output of this function.
    >>> parsed_input = {"locus_name1" : ["gene1", "gene2", "gene3", "geneN"],
    >>>                 "locus_name2" : ["gene1", "gene2", "gene3", "geneN"]}
    """

    locus_genes = defaultdict(lambda: None)
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        for line in lines:
            data = line.strip("\n").split("\t")
            locus = data[1].split()[-1]
            genes = data[2:]
            locus_genes[locus] = genes

    return locus_genes


def _flatten_data(data : dict) -> list:
    """ Flattens data into 1D array

    This is usefull and cleaner processing for the data. embeded lists
    is a result of multple loci being present in the input. This will
    flatten the list of list into one single list.

    Arguments:
    ---------
    data : dict
        Labled data containing the locus name paried with all protein protein
        interaction scores.

    Returns:
    --------
    list
        a flatten list conining all proteins pair interactions scores
    """

    flatten_data = []
    for gene_list in data.values():
        flatten_data += gene_list

    return flatten_data


def save_as_sif(data : dict, out_dir="Simple_out") -> None:
    """ Converts data into sif format.

    Generates one SIF file per locus and also all loci in one SIF. The outputed files
    will contains a unique id that pertains to `Month Day Year - Hours Minutes Seconds'
    format. This allows the user to not only keep track when these files were generated
    but also to prevent overwriting.

    Arguments
    ---------
    data : dict
        Dictioanry containing the locus name and interaction
        scores

    out_dir : str (default="Simple_out")
        Name of the directory where all the generated SIF
        files will be stored

    Returns
    -------
    None
        Writes SIF files
    """
    # NOTE: add comment here 
    unique_id = datetime.today().strftime("%m%d%y-%H%M%S")
    for locus_name, interactions in data.items():
        outname = "{}-{}.sif".format(locus_name, unique_id)
        with open(outname, "w") as sifile:
            sifile.write("gene1 interaction gene2 score\n")
            for interaction in interactions:
                sifile.write("{}\n".format(interaction))

    # NOTE: add comment here how this is different from 
    # here the results are flatten and generates a single SIF file
    interactions_flatten = _flatten_data(data)
    global_outname = "all_nodes-{}.sif".format(unique_id)
    with open(global_outname, "w") as globalout:
        globalout.write("gene1 interaction gene2 score\n")
        for interaction in interactions_flatten:
            globalout.write("{}\n".format(interaction))


    # storing SIF files into directory in the format of
    # "outname-uniqueID"
    dir_name = "{}_{}".format(out_dir, unique_id)
    os.mkdir(dir_name)
    all_sifs = glob.glob("./*-{}.sif".format(unique_id))
    cmd = "mv {} {}".format(" ".join(all_sifs), dir_name)
    print(cmd)
    os.system(cmd)


if __name__ == "__main__":

    # CLI arguments
    description = "simple program for generating SIF to full understand relationships between biological processes."
    parser = argparse.ArgumentParser(description=description)
    required = parser.add_argument_group("Required Arguments")
    required.add_argument("-i", "--input", type=str, metavar="INPUT",
                        help="input file")
    parser.add_argument("-o", "--output", type=str, required=False, default="Simple_out",
                        metavar="PARAMETER", help="name of the outfile default='Simple_out'")
    required.add_argument("-t", "--interaction_type", type=str, metavar="PARAMETER",
                        choices=["pp", "pd", "pr", "rc", "cr", "gl", "pm", "mp"],
                        help="Interaction type")
    parser.add_argument('-db', "--database", type=str, metavar="FILE",
                        help="path to database. Default path is `./Data/STRING.txt",
                        default="./Data/STRING.txt")
    args = parser.parse_args()

    # parsing input and instantiating database object "connecting"
    parsed_input = parse_input(args.input)
    db = StringDB(args.database)

    # collecting all data
    results = defaultdict(lambda: None)
    for locus_name, gene_list in parsed_input.items():

        # each genes list within the parsed data will be send to the database
        # -- scores will be recorded into a dictioanry along with locus name
        result = db.get_pair_score(locus_name, gene_list, interaction_type="pp")
        results[locus_name] = result

    # saving results into SIF files
    save_as_sif(results, out_dir=args.output)

import warnings 
from collections import defaultdict
import argparse
import pandas as pd
import numpy as np

class StringDB:
    def __init__(self, fname):

        # automatically loading in database when instantiating StringDB object
        self.fname = fname
        db = self._load_string_data(fname)
        self.db = db


    def _load_string_data(self, f_path: str) -> pd.DataFrame:
        """ Loads in the STRING database and converts it into a pandas dataframe object"""
        string_df = pd.read_csv(f_path, sep="\t")
        string_df.columns = ["gene1", "gene2", "score"]
        return string_df


    def _cross_ref(self, locus: str, target: list, reference: list) -> None:
        """Cross references matches with initial inputs to see which pairs
        were not found"""

        target_set = set(target)
        ref_set = set(reference)

        
        missing_set = ref_set - target_set
        for missing in missing_set:
            msg = "Not found {} - {}".format(locus, missing)
            warnings.warn(msg)
    

    def _to_adjacency_dict(self, locus, selected_pairs):
        """ converts the selected pairs into a adjacency dict
        
        Arguments
        ---------
        locus: str
            Targeted gene

        selected_pairs: list, np.ndarray
            An array of genes that were identified in the STRING database 
        
        Returns
        -------
        adjacency_dict: dict {str : {str : float}}
            Returns an adjacency dict where the locus is the main key and the
            value is a sub dictionary containing the protein gene and interaction
            score (float)

        >>> # Example result
        >>> adj_dict = {"locus gene":{"matched_gene1":score, "matched_gene2":score, "matched_gene3":score}}
        
        """

        main_result = {} # stores locus and all genes and scores
        gene_score = {} # --> stores gene and score
        for idx in range(len(selected_pairs.index.tolist())):
            data = selected_pairs.iloc[idx]
            gene, score = (data["gene2"], data["score"])
            gene_score[gene] = score

        main_result[locus] = gene_score
        return main_result
    

    def find_pairs(self, locus, genes):
        """ Attempts to find all pairs with a given locus.
        
        arguments
        ----------
        locus: str
            Main gene that will be used to compare all genes
        genes: list, np.array
            An array genes names

        returns
        -------
        dict
            A adjacency_dict containing the locus as the main key and the sub dictionary 
            containing the queryed gene along with its score. If the score is not found,
            then it will not be included in the adjacency_dict
        
        """
        # data type checking
        if not isinstance(locus, str):
            raise TypeError("locus must be a string. you have provided {}". type(locus))
        if not isinstance(genes, list) and not isinstance(genes, np.ndarray):
            raise TypeError("'genes' data must be a list or numpy array, you have provided {}".format(type(genes)))

        # query searches all given genes with one locus
        # -- this process is vectorized does not use a for loop to find every single match
        # -- pairs that are NOT found will not be included in the results
        # -- -- We use the _cross_ref() function to let the user know which pairs where not found
        query = self.db.loc[(self.db["gene1"] == locus) & (self.db["gene2"].isin(genes))] 
        selected_genes = query["gene2"].values.tolist()
        
        # checking for missing pairs
        self._cross_ref(locus, selected_genes, genes)
        
        # return results
        results = self._to_adjacency_dict(locus, query)
        return results


# functions
def save_as_sif(adjacency_dict, interaction, outname, path="."):
    """ Converts adjacency_list into SIF file
    
    Argument
    -------
    adjacency_dict : dict
        Contains all adjacency_dict pathways 
    interaction : str (choices=["pp", "pd", "pr", "rc", "cr", "gl", "pm", "mp"])
        Interaction type of both genes/molecules. Supported interaction types are:
        - p:protein - protein interaction 
        - pd: protein -> DNA
        - pr: protein -> reaction
        - rc: reaction -> compound
        - cr: compound -> reaction
        - gl: genetic lethal relationship
        - pm: protein-metabolite interaction
        - mp: metabolite-protein interaction

    outname : str
        name of the output SIF file
    path : string, optional (default=".")
        Path where the SIF file is going to be written. By default it will be 
        created at current dictory 
    optional

    Returns:
    -------
    ValueError 
        Raised when an incorrect interaction is provided
    File
        SIF file written at provided path (default: current path)
    """
    # type checking 
    known_interaction_types = ["pp", "pd", "pr", "rc", "cr", "gl", "pm", "mp"]
    if interaction not in known_interaction_types:
        raise ValueError("'{}' is an unsupported interaction type. Supported interaction: {}".format(interaction, ", ".join(known_interaction_types)))

    # writing out SIF output file
    full_path = "{}/{}".format(path, outname)
    with open(full_path, "w") as outfile:
        for locus, matches in adjacency_dict.items():
            for gene, score in matches.items():
                result = "{} {} {} {}\n".format(locus, interaction, gene, score)
                outfile.write(result)

def parse_input(input_file):
    locus_genes = defaultdict(None)
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        for line in lines:
            data = line.split("\t")
            locus = data[1].split()[-1]
            genes = data[2:]
            locus_genes[locus] = genes

    return locus_genes


if __name__ == "__main__":

    # CLI arguments
    description = "No description specified"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input", type=str, metavar="INPUT",
                        help="input file")
    parser.add_argument("-o", "--output", type=str, metavar="PARAMETER",
                        help="name of the outfile")
    parser.add_argument("-t", "--interaction_type", type=str, metavar="PARAMETER",
                        help="Interaction type")
    parser.add_argument('-db', "--database", type=str, metavar="FILE", 
                        help="path to database. Default path is `./Data/STRING.txt", 
                        default="./Data/STRING.txt")
    args = parser.parse_args()

    # parsing input and instantiating database object "connecting"
    parsed_input = parse_input(args.input)
    string_db = StringDB(args.database)

    # collecting 
    print(parsed_input)
    
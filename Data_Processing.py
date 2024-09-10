
import os
import xml.etree.ElementTree as ET
import pandas as pd
import pickle
import numpy as np
import pickle
from goatools.obo_parser import GODag

def find_gene_name(gene_names, gene_synonyms):
    """
    Finds the first gene name from a list of gene names that exists in the gene_synonyms dictionary.
    If no direct match is found, it checks for partial matches based on splitting names by '-'.
    
    Parameters:
    gene_names (list): List of possible gene names.
    gene_synonyms (dict): Dictionary of gene names and their synonyms.
    
    Returns:
    str: The key (primary gene name) if a match is found, otherwise None.
    """

    for name in gene_names:
        for key, synonyms in gene_synonyms.items():
            if name in synonyms:
                return key

    for name in gene_names: 
        #Special cases
        possible_names = name.split('-')       
        for possible_name in possible_names:
            for key, synonyms in gene_synonyms.items():
                if possible_name in synonyms:
                    return key

    return None

def process_kegg_gaf(kegg_dir, gaf_file):

    """
    Processes KEGG XML files and a GAF file to extract gene-related information and relations.
    
    Parameters:
    kegg_dir (str): Directory containing KEGG XML files.
    gaf_file (str): Path to the GAF file.
    
    Returns:
    tuple: Contains dictionaries for gene synonyms, disease-gene relationships, graphs of gene relations,
           a gene-GO matrix, list of GO terms, and an array of gene names.
    """

    # Define the column names as per the GAF 2.1 specification
    gaf_columns = [
        "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
        "DB_Reference", "Evidence_Code", "With_From", "Aspect",
        "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type",
        "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"
    ]

    # Read the GAF file, skipping lines that start with '!'
    gaf_df = pd.read_csv(gaf_file, sep='\t', comment='!', header=None, names=gaf_columns, dtype=str)

    # Create a dictionary of gene synonyms
    gene_synonyms = {}

    for gene_name, synonyms_str in zip(gaf_df['DB_Object_Symbol'], gaf_df['DB_Object_Synonym']):
        # Split synonyms if available, otherwise set an empty list
        synonyms = synonyms_str.split('|') if pd.notna(synonyms_str) else []
        synonyms.append(gene_name)  # Add the gene name itself
        gene_synonyms[gene_name] = list(set(synonyms))  # Remove duplicates

    # Dictionary to store gene relations for each disease
    graphs = {}
    # Dictionary to store genes for each disease
    disease_gene_dict = {}

    # Parse each KEGG XML file to extract gene entries and relations
    for kegg_file in os.listdir(kegg_dir):
        if not kegg_file.endswith(".xml"):
            continue

        tree = ET.parse(os.path.join(kegg_dir, kegg_file))
        root = tree.getroot()

        # Get pathway title
        pathway_title = root.attrib.get('title', 'Unknown Pathway')
        disease_gene_dict[pathway_title] = {}

        # Extract gene entries and relations
        for entry in root.findall('entry'):
            if entry.attrib.get('type') == 'gene':
                graphics = entry.find('graphics')
                if graphics is not None and 'name' in graphics.attrib:
                    gene_names = graphics.attrib['name'].replace("...", "").split(',')  # List of possible names
                    valid_gene_name = find_gene_name(gene_names, gene_synonyms)
                    if valid_gene_name:
                        # Accumulate gene occurrences
                        disease_gene_dict[pathway_title][valid_gene_name] = disease_gene_dict[pathway_title].get(valid_gene_name, 0) + 1

        # Graph dictionary for the current disease
        relation_dict = {}

        # Collect relation data
        for relation in root.findall('relation'):
            entry1 = relation.attrib['entry1']
            entry2 = relation.attrib['entry2']
            rel_type = relation.attrib['type']
            subtypes = relation.findall('subtype')

            # Collect edge types
            edge_type = [(rel_type, subtype.attrib['name']) for subtype in subtypes] or [(rel_type, None)]

            # Get node names safely
            entry1_element = root.find(f".//entry[@id='{entry1}']")
            entry2_element = root.find(f".//entry[@id='{entry2}']")
            
            if entry1_element.attrib['type'] != 'group':

                graphics1 = entry1_element.find('graphics')
                node1_names = graphics1.attrib['name'].replace("...", "").split(',')
                node1_name = find_gene_name(node1_names, gene_synonyms) if entry1_element.attrib['type'] == 'gene' else node1_names[0]
            else: 
                node1_name = 'Undefined'

            if entry2_element.attrib['type'] != 'group':

                graphics2 = entry2_element.find('graphics')
                node2_names = graphics2.attrib['name'].replace("...", "").split(',')
                node2_name = find_gene_name(node2_names, gene_synonyms) if entry2_element.attrib['type'] == 'gene' else node2_names[0]
            else: 
                node2_name = 'Undefined'

            relation_dict[(node1_name, node2_name)] = edge_type

        graphs[pathway_title] = relation_dict

    # Create a gene-GO matrix
    genes = set(gene for gene_dict in disease_gene_dict.values() for gene in gene_dict.keys())
    go_terms = gaf_df['GO_ID'].unique()

    # Initialize gene-GO matrix with zeros
    gene_go_matrix = np.zeros((len(genes), len(go_terms)))

    gene_idx = {gene: i for i, gene in enumerate(genes)}
    go_idx = {go: i for i, go in enumerate(go_terms)}

    for gene, go_ids in gaf_df.groupby('DB_Object_Symbol')['GO_ID']:
        if gene in gene_idx:
            for go_id in go_ids.unique():
                if go_id in go_idx:
                    gene_go_matrix[gene_idx[gene], go_idx[go_id]] = 1

    return gene_synonyms, disease_gene_dict, graphs, gene_go_matrix, go_terms, np.array(list(genes))



def parse_obo_to_dict(obo_file):
    """
    Parses an OBO file to extract GO term definitions.
    
    Parameters:
    obo_file (str): Path to the OBO file.
    
    Returns:
    dict: A dictionary with GO term IDs as keys and their definitions as values.
    """
    go_dict = {}
    with open(obo_file, 'r') as f:
        inside_term = False
        current_go_id = None
        for line in f:
            line = line.strip()
            
            # Detect the start of a new term
            if line == "[Term]":
                inside_term = True
                current_go_id = None
            
            # When inside a term, find the GO ID
            if inside_term and line.startswith("id:"):
                current_go_id = line.split("id: ")[1]
            
            # When inside a term, find the definition
            if inside_term and line.startswith("def:"):
                go_definition = line.split("def: ")[1].strip('"')
                if current_go_id:
                    go_dict[current_go_id] = go_definition
                    inside_term = False  # Finish processing this term

    return go_dict



# Paths
kgml_dir = "data/KEGG_data/KGML"
gaf_path = "data/GOA_human/goa_human.gaf"
obo_path = "data/GOA_human/go-basic.obo"

# Process files
gene_synonyms, disease_gene_dict, graphs, gene_go_matrix, go_terms, genes = process_kegg_gaf(kgml_dir, gaf_path)
go_definitions = parse_obo_to_dict(obo_path)

# Create a folder to save the data structures
data_structures = "data/data_structures"
os.makedirs(data_structures, exist_ok=True)

# Save data structures to disk
with open(f"{data_structures}/gene_synonyms.pkl", "wb") as f:
    pickle.dump(gene_synonyms, f)

with open(f"{data_structures}/disease_gene_dict.pkl", "wb") as f:
    pickle.dump(disease_gene_dict, f)

with open(f"{data_structures}/graphs.pkl", "wb") as f:
    pickle.dump(graphs, f)

with open(f"{data_structures}/go_definitions.pkl", "wb") as f:
    pickle.dump(go_definitions, f)

np.save(f"{data_structures}/gene_go_matrix.npy", gene_go_matrix)
np.save(f"{data_structures}/gene_go_matrix.npy", gene_go_matrix)
np.save(f"{data_structures}/go_terms.npy", go_terms)
np.save(f"{data_structures}/genes.npy", genes)

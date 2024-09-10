from langchain.agents import tool
import streamlit as st
import os 
from langchain_openai import ChatOpenAI
import pickle
import numpy as np
from goatools.obo_parser import GODag
import networkx as nx
from config import OPENAI_API_KEY

llm = ChatOpenAI(model='gpt-4o', temperature=0, api_key=OPENAI_API_KEY)

data_directory = 'data/data_structures'

# Define file paths
gene_synonyms_file = os.path.join(data_directory, 'gene_synonyms.pkl')
disease_gene_dict_file = os.path.join(data_directory, 'disease_gene_dict.pkl')
graphs_file = os.path.join(data_directory, 'graphs.pkl')
gene_go_matrix_file = os.path.join(data_directory, 'gene_go_matrix.npy')
go_terms_file = os.path.join(data_directory, 'go_terms.npy')
genes_file = os.path.join(data_directory, 'genes.npy')

# Load the data from pickle files (gene synonyms, disease-gene relationships, graphs)
with open(gene_synonyms_file, 'rb') as f:
    gene_synonyms = pickle.load(f)

with open(disease_gene_dict_file, 'rb') as f:
    disease_gene_dict = pickle.load(f)

with open(graphs_file, 'rb') as f:
    graphs = pickle.load(f)

with open('data/data_structures/go_definitions.pkl', 'rb') as f:
    go_def = pickle.load(f)

# Load numerical data from numpy files (gene-GO matrix, GO terms, and genes)
gene_go_matrix = np.load(gene_go_matrix_file, allow_pickle=True)
go_terms = np.load(go_terms_file, allow_pickle=True)
genes = np.load(genes_file, allow_pickle=True)

# Build directed graphs for each disease using NetworkX
nx_graphs = []  
for disease in disease_gene_dict.keys():
    G = nx.DiGraph() 
    for (gene1, gene2), relations in graphs[disease].items():
        for relation in relations:
            # Add edge for each relation 
            G.add_edge(gene1, gene2, relation=relation)
    nx_graphs.append(G)

# Load Gene Ontology DAG structure
godag = GODag("data/GOA_human/go-basic.obo")

# Function to find the correct gene name using a gene synonyms dictionary
def find_gene(gene_name, gene_synonyms):
    for key, synonyms in gene_synonyms.items():
        if gene_name in synonyms:
            return key
    return None

@tool
def find_diseases_for_gene(gene_name):
    """
    This function takes a gene name or its synonym, queries the dictionary of gene synonyms, 
    and returns a list of diseases (pathways) that the gene is associated with.

    Parameters:
    gene_name (str): The name of the gene to search for.

    Returns:
    str: A list of diseases the gene is involved in, or a message if no disease is found.    
    """
    # Find the correct gene name using the synonym dictionary
    gene_name_lower = gene_name.lower()
    correct_gene_name = None

    # Check if the input gene name matches either the key or any of the synonyms
    for key, synonyms in gene_synonyms.items():
        if gene_name_lower == key.lower() or any(gene_name_lower == synonym.lower() for synonym in synonyms):
            correct_gene_name = key
            break
    
    # If the gene is found, look for it in disease-gene dictionary
    if correct_gene_name is not None:
        diseases = [disease for disease, gene_dict in disease_gene_dict.items() if correct_gene_name in gene_dict]

        if len(diseases) > 0:
            return ', '.join(diseases)
    
    return 'No disease is associated with this gene in the database.'

@tool
def downstream_analysis(disease, start_gene):

    """
    Perform downstream analysis on a given gene and a specific disease. 
    This function finds all paths from the start gene to the leaf nodes in the disease-specific pathway graph, and 
    returns the most common biological properties and molecular functions shared across the genes in each path.

    Parameters:
    disease (str): The disease (pathway) to analyze.
    start_gene (str): The gene from which to start the analysis.

    Returns:
    str: A detailed report of the paths and common biological properties between the genes.
     """

    # Find the correct gene name using the synonym dictionary
    correct_gene_name = find_gene(start_gene, gene_synonyms)
    if correct_gene_name is None: return  f"Gene {start_gene} not found in the data."
    # Check if the graph contains the start gene
    if correct_gene_name not in disease_gene_dict[disease]:
        return f"Gene {start_gene} is not associated with the disease {disease}."

    disease_index = list(disease_gene_dict.keys()).index(disease)

    # Find all paths from the input gene to the leaf nodes
    all_paths = []
    graph = nx_graphs[disease_index]
    
    # Explore all paths
    for end_gene in graph.nodes():
        if graph.out_degree(end_gene) == 0:  # Leaf node (no outgoing edges)
            # Find all paths from start_gene to the leaf end_gene
            paths = list(nx.all_simple_paths(graph, source=correct_gene_name, target=end_gene))
            all_paths.extend(paths)
    
    # Analyze each path and retrieve the most common biological properties
    ans = ''
    for i, path in enumerate(all_paths):
        ans += f"Path {i + 1}: \n"
        gene_indexes = []
        # Iterate over each pair of consecutive genes in the path
        for i in range(len(path) - 1):
            node1 = path[i]
            node2 = path[i + 1]
            # Get the edge data between node1 and node2
            edge_data = graph.get_edge_data(node1, node2)['relation']
            # Add the first node and edge information to the string
            ans += f"{node1} -- {edge_data} --> "

            # Add the gene indexes for the gene-GO matrix
            if node1 in genes: gene_indexes.append(np.where(genes == node1)[0][0])
            if node2 in genes: gene_indexes.append(np.where(genes == node2)[0][0])

        # Add the final gene in the path (leaf node)
        ans += f'{path[-1]} \n'

        # Analyze GO terms for genes in the path
        go_term_count = np.sum(gene_go_matrix[gene_indexes, :], axis=0)
        sorted_go_idx = np.argsort(go_term_count)

        ans += f'The most common biological properties and molecular functions shared across genes in this path are: \n'
        i, c = 1, 0
        done = False
        while not done:
            go_id = go_terms[sorted_go_idx[-i]]
            if godag[go_id].namespace != 'cellular_component':
                ans += f'{godag[go_id].name} ({godag[go_id].namespace}). Definition: {go_def[go_id]}  \n'
                c += 1

            if c == 5 or i >= 50 or go_term_count[sorted_go_idx[-i]] == 0: 
                done = True
                if c == 1: ans += 'No common biological properties and molecular functions found. \n\n'

            i += 1
        
    return ans

@tool
def collective_impact_analysys(gene_names_string):
    """
    Retrieves the collective impact of genes in biological processes and disease outcomes.
    This function takes a string of gene names and for each disease finds all 
    paths where all the genes are present. For each path, it retrieves the most common biological 
    properties and molecular functions shared across the genes in the path.
    
    Parameters:
    gene_names_string (str): A comma-separated string of gene names.
    
    Returns:
    str: A report of the paths containing all genes and their shared biological properties.
    """
    
    # Convert the input string into a list of gene names
    gene_names = [gene.strip() for gene in gene_names_string.split(',')]
    
    # Map gene names to their correct keys using the gene_synonyms dictionary
    correct_gene_names = []
    for gene_name in gene_names:
        correct_gene = find_gene(gene_name, gene_synonyms)
        if correct_gene is None:
            return f"Gene {gene_name} not found in the data."
        correct_gene_names.append(correct_gene)
    
    result = ""
    
    # For each disease, check if the paths contain all the genes from the list
    for disease, gene_dict in disease_gene_dict.items():
        result += f"\nDisease: {disease}\n"
        disease_index = list(disease_gene_dict.keys()).index(disease)
        graph = nx_graphs[disease_index]
        
        all_paths_with_genes = []
        
        # Find all paths in the disease graph that contain all the genes in the list
        for start_gene in graph.nodes():
            if start_gene in correct_gene_names:
                for end_gene in graph.nodes():
                    if graph.out_degree(end_gene) == 0:  # Leaf node (no outgoing edges)
                        # Find all paths from start_gene to the leaf end_gene
                        paths = list(nx.all_simple_paths(graph, source=start_gene, target=end_gene))
                        
                        # Check if all genes in the input list are present in the path
                        for path in paths:
                            if all(gene in path for gene in correct_gene_names):
                                all_paths_with_genes.append(path)
        
        # Analyze each path and retrieve the most common biological properties
        for i, path in enumerate(all_paths_with_genes):
            result += f"Path {i + 1}: \n"
            gene_indexes = []
            
            # Iterate over each pair of consecutive genes in the path
            for j in range(len(path) - 1):
                node1 = path[j]
                node2 = path[j + 1]
                edge_data = graph.get_edge_data(node1, node2)['relation']
                result += f"{node1} -- {edge_data} --> "
                
                if node1 in genes: gene_indexes.append(np.where(genes == node1)[0][0])
                if node2 in genes: gene_indexes.append(np.where(genes == node2)[0][0])

            # Add the final gene in the path (leaf node)
            result += f'{path[-1]} \n'

            # Analyze GO terms for genes in the path
            go_term_count = np.sum(gene_go_matrix[gene_indexes, :], axis=0)
            sorted_go_idx = np.argsort(go_term_count)

            result += 'The most common biological properties and molecular functions shared across genes in this path are: \n'
            c = 0
            for k in range(1, 51):  # Consider the top 50 GO terms
                go_id = go_terms[sorted_go_idx[-k]]
                if godag[go_id].namespace != 'cellular_component':
                    result += f'{godag[go_id].name} ({godag[go_id].namespace}). Definition: {go_def[go_id]} \n'
                    c += 1
                if c == 5:  # Limit to the top 5 biological properties
                    break

            if c == 0:
                result += 'No common biological properties and molecular functions found.\n'
            result += '\n'

        if len(all_paths_with_genes) == 0:
            result += "No paths found containing all the genes in this disease graph.\n"
    
    return result


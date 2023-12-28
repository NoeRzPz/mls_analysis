from bioservices import QuickGO
import csv
import os

# Define your data directories
dataDir = '~/mls_analysis/data/'

# Expand the user's home directory in the directories path
expanded_dataDir = os.path.expanduser(dataDir)

go_terms = {
    "Genomic Instability": [("DNA_repair", "GO:0006281"), ("Nuclear_lamina", "GO:0005652")],
    "Telomere Attrition": [("Telomere_maintenance", "GO:0000723"), ("Telomere_capping", "GO:0016233")],
    "Epigenetic Alterations": [("Histone modification", "GO:0016570"), ("DNA methylation", "GO:0006306"), ("Chromatine remodeling", "GO:0006338")],
    "Loss of Proteostasis ": [("Chaperone-mediated protein folding", "GO:0061077"), ("Autophagy", "GO:0006914"), ("Ubiquitin-proteasome system", "GO:0043161")],
    "Deregulated Nutrient Sensing": [("Insulin receptor signaling Pathway", "GO:0008286"), ("TOR signaling", "GO:0031929")],
    "Mitochondrial Dysfunction": [("Response to ROS ", "GO:0000302"), ("Mitochondrial genome maintenance", "GO:0000002")],
    "Cellular senescence": [("Cellular senescence", "GO:0090398")],
    "Stem Cell Exhaustion": [("Stem cell proliferation", "GO:0072089")],
    "Altered Intracellular Communication": [("Inflammatory response", "GO:0006954"), ("Inflammasome complex", "GO:0061702")]
}

#curl -s http://www.ebi.ac.uk/QuickGO/GAnnotation?tax=9606&relType=IP=&goid=GO:0003015&format=tsv

#If you require a list of the proteins used in these annotations, you can use the following URL;

#curl -s http://www.ebi.ac.uk/QuickGO/GAnnotation?tax=9606&relType=IP=&goid=GO:0003015&format=proteinList

def get_genes_for_go_term(go_term):
    quickgo = QuickGO()
    response = quickgo.Annotation(goId=go_term, taxonId="9606") # 9606 is the taxon ID for Homo sapiens
    
    # Check if response is valid and contains data
    if 'results' in response:
        genes = [result['symbol'] for result in response['results'] if 'symbol' in result]
        return list(set(genes))  # Remove duplicates
    else:
        return []  # Return an empty list if no results are found


gmt_file_path = os.path.join(expanded_dataDir, "agein.Hs.symbols.gmt")

with open(gmt_file_path, "w", newline="") as file:
    writer = csv.writer(file, delimiter='\t')
    for category, terms in go_terms.items():
        for term_name, go_id in terms:
            genes = get_genes_for_go_term(go_id)
            #writer.writerow([term_name, category] + genes)
            writer.writerow([term_name] + genes)
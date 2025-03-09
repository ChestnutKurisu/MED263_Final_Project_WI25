import pandas as pd
import csv
import os

# -------------------------------
# PARAMETERS & SETTINGS
# -------------------------------
# Expanded gene list: original core genes plus additional high-value candidates.
genes_of_interest = [
    "Unnamed: 0",  # This column holds the cell-line identifier (ModelID)
    "TP53",       # Tumor suppressor ("guardian of the genome"); mutated in many cancers (breast, lung, colon, etc.)
    "EGFR",       # Receptor tyrosine kinase; often overexpressed/mutated in non-small cell lung cancer (NSCLC), glioblastoma, and colon cancer
    "KRAS",       # Oncogene; commonly mutated in pancreatic, colorectal, and lung cancers
    "BRCA1",      # Tumor suppressor involved in DNA repair; mutations predispose to hereditary breast and ovarian cancer
    "BRCA2",      # Similar to BRCA1; mutations increase risk for breast, ovarian, and pancreatic cancers
    "CDKN2A",     # Tumor suppressor that regulates the cell cycle; frequently altered in melanoma and pancreatic cancer
    "PTEN",       # Negative regulator of the PI3K/AKT pathway; loss-of-function mutations seen in glioblastoma, endometrial, and prostate cancers
    "MYC",        # Oncogene regulating cell proliferation; amplification occurs in cancers like Burkitt lymphoma, breast, and lung cancers
    "PIK3CA",     # Catalytic subunit of PI3K; mutations drive oncogenesis in breast, colorectal, and endometrial cancers
    "RB1",        # Tumor suppressor; its inactivation is central to retinoblastoma and contributes to other cancers (e.g., osteosarcoma)
    "BRAF",       # Oncogene; the V600E mutation is a key driver in melanoma, thyroid, and colorectal cancers
    # Additional cancer-implicated genes:
    "APC",        # Tumor suppressor; mutations are the hallmark of familial adenomatous polyposis (FAP) and sporadic colorectal cancers
    "CTNNB1",     # Encodes β-catenin; involved in the Wnt signaling pathway and implicated in hepatocellular carcinoma, colorectal, and endometrial cancers
    "SMAD4",      # Tumor suppressor in the TGF-β signaling pathway; frequently inactivated in pancreatic ductal adenocarcinoma and other cancers
    "ERBB2",      # Also known as HER2; receptor tyrosine kinase whose amplification is associated with aggressive breast and gastric cancers
    "ALK",        # Oncogene; rearrangements and fusions (e.g., EML4-ALK) drive tumorigenesis in non-small cell lung cancer and anaplastic large cell lymphoma
    "MET",        # Receptor tyrosine kinase; mutations/amplifications implicated in lung, kidney, and other cancers
    "NOTCH1",     # Part of the Notch signaling pathway; mutations are linked to T-cell acute lymphoblastic leukemia (T-ALL) and other malignancies
    "NOTCH2",     # Notch family member; aberrations have been observed in certain B-cell malignancies and other cancers
    "NOTCH3",     # Notch receptor; overexpression/mutations are associated with ovarian and lung cancers
    "NOTCH4",     # Less frequently altered Notch receptor; dysregulation can contribute to tumorigenesis in select cancers
    "NF1",        # Tumor suppressor gene; mutations cause neurofibromatosis type 1 and predispose to malignant peripheral nerve sheath tumors
    "NF2",        # Tumor suppressor; loss-of-function mutations lead to neurofibromatosis type 2, often causing vestibular schwannomas and meningiomas
    "CDK4",       # Cyclin-dependent kinase; alterations drive cell cycle progression and are implicated in melanoma and other cancers
    "CDK6",       # Works with CDK4 in cell cycle regulation; dysregulation observed in lymphoma, breast cancer, etc.
    "IDH1",       # Metabolic enzyme; mutations (affecting cellular metabolism and epigenetics) are common in gliomas and acute myeloid leukemia (AML)
    "IDH2",       # Similar to IDH1; mutations are implicated in AML and certain gliomas
    "ARID1A",     # Component of the SWI/SNF chromatin remodeling complex; frequently mutated in ovarian clear cell carcinoma and endometrial cancers
    "CTCF"        # Transcriptional regulator involved in chromatin organization; alterations can disrupt genome architecture and gene expression in various cancers
]

# Input folder for full DepMap files (modify as needed)
input_dir = r"C:\Users\param\OneDrive\Documents\PyCharmProjects\UCSD-CSE\MED 263\final-project\CRISPR-multi-omics\DepMap"
# Output folder for subset files
output_dir = r"C:\Users\param\OneDrive\Documents\PyCharmProjects\UCSD-CSE\MED 263\final-project\CRISPR-multi-omics\DepMap-preprocessed"
os.makedirs(output_dir, exist_ok=True)

def fix_header(header):
    """If the first column name is blank, set it to 'Unnamed: 0'."""
    if header and header[0].strip() == "":
        header[0] = "Unnamed: 0"
    return header

# -------------------------------
# 1. Subset CRISPRGeneEffect.csv
# -------------------------------
effect_input = os.path.join(input_dir, "CRISPRGeneEffect.csv")
effect_output = os.path.join(output_dir, "CRISPRGeneEffect.csv")

with open(effect_input, 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    header_effect = fix_header(next(reader))

selected_effect_cols = [col for col in header_effect if col == 'Unnamed: 0' or col.split(' ')[0] in genes_of_interest]
print("CRISPRGeneEffect: Selected columns:", selected_effect_cols)

effect_chunks = []
for chunk in pd.read_csv(effect_input, usecols=selected_effect_cols, chunksize=5000):
    effect_chunks.append(chunk)
gene_effect_subset = pd.concat(effect_chunks, axis=0)
# Rename the identifier column (whether it is "Unnamed: 0" or blank) to "ModelID"
gene_effect_subset = gene_effect_subset.rename(columns={"Unnamed: 0": "ModelID", "": "ModelID"})
gene_effect_subset.to_csv(effect_output, index=False)
print("Saved subset of CRISPRGeneEffect to", effect_output)

# -------------------------------
# 2. Subset CRISPRGeneDependency.csv
# -------------------------------
dep_input = os.path.join(input_dir, "CRISPRGeneDependency.csv")
dep_output = os.path.join(output_dir, "CRISPRGeneDependency.csv")

with open(dep_input, 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    header_dep = fix_header(next(reader))
selected_dep_cols = [col for col in header_dep if col == 'Unnamed: 0' or col.split(' ')[0] in genes_of_interest]
print("CRISPRGeneDependency: Selected columns:", selected_dep_cols)

dep_chunks = []
for chunk in pd.read_csv(dep_input, usecols=selected_dep_cols, chunksize=5000):
    dep_chunks.append(chunk)
gene_dep_subset = pd.concat(dep_chunks, axis=0)
gene_dep_subset = gene_dep_subset.rename(columns={"Unnamed: 0": "ModelID", "": "ModelID"})
gene_dep_subset.to_csv(dep_output, index=False)
print("Saved subset of CRISPRGeneDependency to", dep_output)

# -------------------------------
# 3. Subset OmicsCNGene.csv
# -------------------------------
cn_input = os.path.join(input_dir, "OmicsCNGene.csv")
cn_output = os.path.join(output_dir, "OmicsCNGene.csv")

with open(cn_input, 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    header_cn = fix_header(next(reader))
selected_cn_cols = [col for col in header_cn if col == 'Unnamed: 0' or col.split(' ')[0] in genes_of_interest]
if "ModelID" not in selected_cn_cols and "Unnamed: 0" in header_cn:
    selected_cn_cols.insert(0, "Unnamed: 0")
print("OmicsCNGene: Selected columns:", selected_cn_cols)

cn_chunks = []
for chunk in pd.read_csv(cn_input, usecols=selected_cn_cols, chunksize=5000):
    cn_chunks.append(chunk)
cn_subset = pd.concat(cn_chunks, axis=0)
cn_subset = cn_subset.rename(columns={"Unnamed: 0": "ModelID", "": "ModelID"})
cn_subset.to_csv(cn_output, index=False)
print("Saved subset of OmicsCNGene to", cn_output)

# -------------------------------
# 4. Subset OmicsSomaticMutationsMatrixDamaging.csv
# -------------------------------
mut_input = os.path.join(input_dir, "OmicsSomaticMutationsMatrixDamaging.csv")
mut_output = os.path.join(output_dir, "OmicsSomaticMutationsMatrixDamaging.csv")

with open(mut_input, 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    header_mut = fix_header(next(reader))
selected_mut_cols = [col for col in header_mut if col == 'Unnamed: 0' or col.split(' ')[0] in genes_of_interest]
if "Unnamed: 0" not in selected_mut_cols:
    selected_mut_cols.insert(0, "Unnamed: 0")
print("OmicsSomaticMutationsMatrixDamaging: Selected columns:", selected_mut_cols)

mut_chunks = []
for chunk in pd.read_csv(mut_input, usecols=selected_mut_cols, chunksize=5000):
    mut_chunks.append(chunk)
mut_subset = pd.concat(mut_chunks, axis=0)
mut_subset = mut_subset.rename(columns={"Unnamed: 0": "ModelID", "": "ModelID"})
mut_subset.to_csv(mut_output, index=False)
print("Saved subset of OmicsSomaticMutationsMatrixDamaging to", mut_output)

# -------------------------------
# 5. Subset OmicsExpressionProteinCodingGenesTPMLogp1.csv
# -------------------------------
expr_input = os.path.join(input_dir, "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
expr_output = os.path.join(output_dir, "OmicsExpressionProteinCodingGenesTPMLogp1.csv")

with open(expr_input, 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    header_expr = fix_header(next(reader))
# Select columns from the expanded list if present in the expression file header.
selected_expr_cols = [col for col in header_expr if col == 'Unnamed: 0' or col.split(' ')[0] in genes_of_interest]
if "Unnamed: 0" not in selected_expr_cols and "Unnamed: 0" in header_expr:
    selected_expr_cols.insert(0, "Unnamed: 0")
print("OmicsExpressionProteinCodingGenesTPMLogp1: Selected columns:", selected_expr_cols)

expr_chunks = []
for chunk in pd.read_csv(expr_input, usecols=selected_expr_cols, chunksize=5000):
    expr_chunks.append(chunk)
expr_subset = pd.concat(expr_chunks, axis=0)
expr_subset = expr_subset.rename(columns={"Unnamed: 0": "ModelID", "": "ModelID"})
expr_subset.to_csv(expr_output, index=False)
print("Saved subset of OmicsExpressionProteinCodingGenesTPMLogp1 to", expr_output)

# -------------------------------
# 6. Subset OmicsExpressionProteinCodingGenesTPMLogp1Stranded.csv
# -------------------------------
expr_stranded_input = os.path.join(input_dir, "OmicsExpressionProteinCodingGenesTPMLogp1Stranded.csv")
expr_stranded_output = os.path.join(output_dir, "OmicsExpressionProteinCodingGenesTPMLogp1Stranded.csv")

with open(expr_stranded_input, 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    header_expr_stranded = fix_header(next(reader))
selected_expr_stranded_cols = [col for col in header_expr_stranded if col == 'Unnamed: 0' or col.split(' ')[0] in genes_of_interest]
if "Unnamed: 0" not in selected_expr_stranded_cols and "Unnamed: 0" in header_expr_stranded:
    selected_expr_stranded_cols.insert(0, "Unnamed: 0")
print("OmicsExpressionProteinCodingGenesTPMLogp1Stranded: Selected columns:", selected_expr_stranded_cols)

expr_stranded_chunks = []
for chunk in pd.read_csv(expr_stranded_input, usecols=selected_expr_stranded_cols, chunksize=5000):
    expr_stranded_chunks.append(chunk)
expr_stranded_subset = pd.concat(expr_stranded_chunks, axis=0)
expr_stranded_subset = expr_stranded_subset.rename(columns={"Unnamed: 0": "ModelID", "": "ModelID"})
expr_stranded_subset.to_csv(expr_stranded_output, index=False)
print("Saved subset of OmicsExpressionProteinCodingGenesTPMLogp1Stranded to", expr_stranded_output)

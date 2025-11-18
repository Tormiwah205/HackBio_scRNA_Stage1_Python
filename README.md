# Adeyeye Daniel Oludayo
# HackBio Internship
## Project: Stage 1 of scRNA Internship: Coding for Bio (Python)

### Tasks
1. **DNA to Protein Translation**  
   A Python function that translates a DNA sequence into the corresponding protein sequence.

2. **Hamming Distance Calculation**  
   A Python function to calculate the Hamming distance between my Slack username and my Twitter/X handle.

3. **Literature Review**  
   Collaborate with my team to review one of the assigned papers and produce a concise essay of 300–400 words summarizing key findings and insights.

4. **Reproduce Figures**  
   Recreate figures from the selected study using Python data libraries (**pandas**, **seaborn**, and optionally **matplotlib**). This includes scatter plots, heatmaps, volcano plots, and density plots to visualize gene expression and breast cancer dataset features.

## 1. DNA to Protein Translation

```python

# TASK 1: DNA to Protein Translation

# RNA codon table for translating RNA codons into amino acids
codon_table = {
    'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L',
    'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S',
    'UAU':'Y', 'UAC':'Y', 'UAA':'_', 'UAG':'_',
    'UGU':'C', 'UGC':'C', 'UGA':'_', 'UGG':'W',
    'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',
    'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    'CGU':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'AUU':'I', 'AUC':'I', 'AUA':'I', 'AUG':'M',
    'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    'AGU':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V',
    'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
}

def dna_to_protein(dna):
    """Convert a DNA sequence into a protein string."""
    dna = dna.upper()  # standardize to uppercase

    # Validate DNA sequence contains only A, T, G, C
    for base in dna:
        if base not in "ATGC":
            return "INVALID"

    # Convert DNA to RNA by replacing T with U
    rna = dna.replace("T", "U")

    protein = ""
    # Translate RNA in codons (3 nucleotides each)
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        protein += codon_table.get(codon, "X")  # X for unknown codons

    return protein

# Example DNA sequence
dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
protein_seq = dna_to_protein(dna_seq)

# Print protein or error if invalid
if protein_seq == "INVALID":
    print("Error: DNA sequence contains invalid characters.")
else:
    print("DNA:", dna_seq)
    print("Protein:", protein_seq)

```

## 2. Hamming Distance Calculation

```python

# TASK 2: Hamming Distance Calculation
def hamming_distance(s1, s2):
    """Compute Hamming distance by padding shorter string."""
    max_len = max(len(s1), len(s2))  # make lengths equal
    s1 = s1.ljust(max_len, "_")
    s2 = s2.ljust(max_len, "_")

    distance = 0
    # Count differing characters
    for a, b in zip(s1, s2):
        if a != b:
            distance += 1

    return distance

# Example strings to compare
slack = "daniel"
twitter = "danolud"

print("Slack:", slack)
print("Twitter/X:", twitter)
print("Hamming distance:", hamming_distance(slack, twitter))

```

## 3. Literature Review: Single-Cell RNA-seq Analysis

This is a review of the assigned essay on "Best Practices in Single-Cell RNA-seq Analysis: A Tutorial".

## Best Practices in Single-Cell RNA-seq Analysis: A Tutorial

Single-cell RNA sequencing (scRNA-seq) has fundamentally advanced the study of biological systems, enabling the analysis of gene expression at an unprecedented resolution and leading to the discovery of previously obscured cellular populations. However, the relative immaturity of the field and this proliferation of methods present significant challenges to establishing standardized scRNA-seq analysis workflows. To address this, current best practices outline a comprehensive workflow that progresses from raw count matrices through downstream analysis. 

The initial stage, pre-processing, involves several steps: quality control (QC), normalization, data correction, feature selection, and dimensionality reduction. Cell QC is crucial, requiring the joint examination of key covariates, to identify and filter out outliers like dying cells or doublets. For normalization, which aims to correct for variability in captured mRNA molecules due to sampling effects, the scran pooling-based size factor estimation method is highly recommended for non-full-length datasets, consistently demonstrating top performance in comparison studies. Subsequently, data correction addresses unwanted variability; for straightforward batch effects within a single experiment, the linear method ComBat is recommended as a strong performer. A critical recommendation is to perform statistical testing on "measured data" rather than on "corrected data," which is better reserved for visual comparison. For efficiency in downstream analysis, dimensionality reduction is performed. For visualization, non-linear methods are preferred; UMAP (Uniform Approximation and Projection) is regarded as best practice for exploratory data visualization due to its speed and ability to scale to large datasets.

Downstream analysis focuses on biological inference through cell-level and gene-level approaches. Cell clustering, often achieved using Louvain community detection on a K-Nearest Neighbour (KNN) graph, has been shown to be robust for scRNA-seq data. Trajectory inference (TI) models continuous biological processes; the Slingshot method is highly effective for reconstructing simpler linear or bifurcating trajectories. Finally, differential expression (DE) testing is performed, typically within cell-identity clusters. Recommended tools include the single-cell specific hurdle model MAST or bulk methods adapted for single-cell data.

### Reference

Luecken, M. D., & Theis, F. J. (2019). Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology, 15*(e8746). https://doi.org/10.15252/msb.20188746

## 4. Reproduce Figures
#### Figures 1 through 6 display the results of the code, allowing you to review the outputs of each analysis step.

```python

# TASK 2
# PART A: Gene Expression Analysis (HBR vs UHR)
# This section performs two analyses: 1. Clustered heatmap of the top differentially expressed genes; 2. Volcano plot of differential expression results.

# Load normalized gene expression matrix of top DEGs.
file_path = "hbr_uhr_top_deg_normalized_counts.csv"
df = pd.read_csv(file_path, index_col=0)

print("Loaded normalized expression data.")
print("Shape:", df.shape)

# Check for missing values, as they can interfere with clustering.
if df.isnull().values.any():
    print("Warning: Missing values detected.")

# A(a): Clustered Heatmap of Top DEGs
# Create a hierarchical clustered heatmap to reveal similarity patterns among samples and genes based on their expression levels.
g = sns.clustermap(
    df,
    cmap="Blues",            # Visual scale for expression intensity
    linewidths=0.4,
    linecolor="black",
    figsize=(12, 10),
    xticklabels=True,
    yticklabels=True
)

g.fig.suptitle(
    "Clustered Heatmap of Top Differentially Expressed Genes (HBR vs UHR)",
    fontsize=14,
    fontweight="bold",
    y=1
)
plt.show()

# A(b): Volcano Plot of Differential Expression Results
# This plot visualizes the relationship between fold change and statistical significance.

diff_path = "hbr_uhr_deg_chr22_with_significance.csv"
deg_df = pd.read_csv(diff_path)
print("\nLoaded DEG data. Shape:", deg_df.shape)

# Ensure required columns are present.
if not all(col in deg_df.columns for col in ["log2FoldChange", "PAdj"]):
    raise ValueError("Missing one or more required columns: log2FoldChange, PAdj")

# Remove rows with missing or invalid values.
deg_df = deg_df.dropna(subset=["log2FoldChange", "PAdj"])
deg_df = deg_df[deg_df["PAdj"] > 0]

# Convert adjusted p-values to -log10 scale to better visualize significance.
deg_df["neg_log10_padj"] = -np.log10(deg_df["PAdj"])

# Classify genes based on fold change and significance thresholds.
# Thresholds used:- |log2FC| > 1 : meaningful biological change; PAdj < 0.05 : statistically significant change

def label_gene(row):
    if row["PAdj"] < 0.05 and row["log2FoldChange"] > 1:
        return "Upregulated"
    elif row["PAdj"] < 0.05 and row["log2FoldChange"] < -1:
        return "Downregulated"
    return "Not significant"

deg_df["Significance"] = deg_df.apply(label_gene, axis=1)

# Colors used to distinguish categories in the volcano plot.
colors = {
    "Upregulated": "green",
    "Downregulated": "orange",
    "Not significant": "grey"
}

plt.figure(figsize=(10, 7))

# Plot each gene category separately for clarity.
for category, color in colors.items():
    subset = deg_df[deg_df["Significance"] == category]
    plt.scatter(
        subset["log2FoldChange"],
        subset["neg_log10_padj"],
        color=color,
        label=category,
        alpha=0.7,
        edgecolor="black",
        linewidth=0.3
    )

# Add threshold reference lines.
plt.axvline(x=1, color="black", linestyle="--")    # Upregulation threshold
plt.axvline(x=-1, color="black", linestyle="--")   # Downregulation threshold
plt.axhline(y=0, color="grey", linestyle="--")     # Baseline significance

plt.title("Volcano Plot of Differential Expression Results (HBR vs UHR)", fontsize=14, fontweight="bold", pad=10)
plt.xlabel("log2(Fold Change)")
plt.ylabel("-log10(Adjusted p-value)")
plt.legend(title="Significance")
plt.tight_layout()
plt.show()

# PART B: Breast Cancer Data Exploration (Tasks C–F)
# This section explores a breast cancer dataset to visualize relationships among tumor features and compare malignant vs benign cases.

breast_cancer_file = "data-3.csv"
df = pd.read_csv(breast_cancer_file)
print("Dataset loaded. Shape:", df.shape)

# B(c): Scatter plot (texture_mean vs radius_mean)
# Purpose: Examine how these two features separate malignant from benign tumors.

required_cols_c = ["diagnosis", "radius_mean", "texture_mean"]
missing_c = [c for c in required_cols_c if c not in df.columns]

if missing_c:
    print("Missing columns for scatter plot:", missing_c)
else:
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x="radius_mean",
        y="texture_mean",
        hue="diagnosis",     # Distinguish malignant vs benign
        palette={"M": "red", "B": "blue"},
        alpha=0.7,
        edgecolor="black"
    )
    plt.title("Scatter Plot of Texture vs Radius (Breast Cancer Data)", fontsize=14, fontweight="bold")
    plt.xlabel("Radius Mean")
    plt.ylabel("Texture Mean")
    plt.tight_layout()
    plt.show()

# B (d): Correlation heatmap of selected features
# Purpose: Identify how key morphological features relate to each other.

features_d = [
    "radius_mean",
    "texture_mean",
    "perimeter_mean",
    "area_mean",
    "smoothness_mean",
    "compactness_mean"
]

missing_d = [c for c in features_d if c not in df.columns]

if missing_d:
    print("Missing columns for correlation heatmap:", missing_d)
else:
    df_corr = df[features_d].corr()
    plt.figure(figsize=(8, 6))
    sns.heatmap(df_corr, annot=True, cmap="Blues", fmt=".2f", linewidths=0.5)
    plt.title("Correlation Heatmap of Breast Cancer Features", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.show()


# B (e) Scatter Plot: smoothness_mean vs compactness_mean
# Purpose: Observe whether these shape descriptors help separate tumor classes.

required_cols_e = ["diagnosis", "smoothness_mean", "compactness_mean"]
missing_e = [c for c in required_cols_e if c not in df.columns]

if missing_e:
    print("Missing columns for scatter plot:", missing_e)
else:
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x="smoothness_mean",
        y="compactness_mean",
        hue="diagnosis",
        palette={"M": "red", "B": "blue"},
        alpha=0.7,
        edgecolor="black"
    )
    plt.title("Scatter Plot of Compactness vs Smoothness", fontsize=14, fontweight="bold")
    plt.xlabel("Smoothness Mean")
    plt.ylabel("Compactness Mean")
    plt.tight_layout()
    plt.show()

# B(f) Density Plot: area_mean distribution in malignant vs benign tumors
# Purpose: Compare how tumor area differs between the two classes.

required_cols_f = ["diagnosis", "area_mean"]
missing_f = [c for c in required_cols_f if c not in df.columns]

if missing_f:
    print("Missing columns for density plot:", missing_f)
else:
    plt.figure(figsize=(8, 6))

    # Density curve for malignant tumors
    sns.kdeplot(
        data=df[df["diagnosis"] == "M"],
        x="area_mean",
        label="Malignant",
        color="red",
        fill=True
    )

    # Density curve for benign tumors
    sns.kdeplot(
        data=df[df["diagnosis"] == "B"],
        x="area_mean",
        label="Benign",
        color="blue",
        fill=True
    )

    plt.title("Density Plot of Area Mean (Malignant vs Benign)", fontsize=14, fontweight="bold")
    plt.xlabel("Area Mean")
    plt.ylabel("Density")
    plt.legend(title="Diagnosis")
    plt.tight_layout()
    plt.show()

```

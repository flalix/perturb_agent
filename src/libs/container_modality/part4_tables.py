"""Part 4: presentation views.
  4.1 only central (unimodal, k=1) distributions
  4.2 ordered by maximum total up modality
  4.3 ordered by maximum total down modality
4.2/4.3 are given at both levels: genes (samples falling in up/down modes) and samples (genes in up/down classes).
"""
import pandas as pd

def central_only(gene_summary):                                   # 4.1
    return gene_summary[gene_summary["k"] == 1]

def genes_by_up(gene_summary):                                    # 4.2 gene level
    return gene_summary.sort_values(["n_up_samples", "n_up_modes"], ascending=False)

def genes_by_down(gene_summary):                                  # 4.3 gene level
    return gene_summary.sort_values(["n_down_samples", "n_down_modes"], ascending=False)

def samples_by_up(sample_counts):                                 # 4.2 sample level
    return sample_counts.sort_values(["total_up", "+2 up"], ascending=False)

def samples_by_down(sample_counts):                               # 4.3 sample level
    return sample_counts.sort_values(["total_down", "-2 down"], ascending=False)

if __name__ == "__main__":
    pd.set_option("display.width", 200)
    gs = pd.read_csv("gene_summary.csv", index_col=0)
    sc = pd.read_csv("sample_class_counts.csv", index_col=0)
    print("4.1 central-only genes\n", central_only(gs).head(), "\n")
    print("4.2 genes by total up\n", genes_by_up(gs).head(), "\n")
    print("4.2 samples by total up\n", samples_by_up(sc).head(), "\n")
    print("4.3 genes by total down\n", genes_by_down(gs).head(), "\n")
    print("4.3 samples by total down\n", samples_by_down(sc).head())
    central_only(gs).to_csv("t41_central_only.csv")
    genes_by_up(gs).to_csv("t42_genes_by_up.csv"); samples_by_up(sc).to_csv("t42_samples_by_up.csv")
    genes_by_down(gs).to_csv("t43_genes_by_down.csv"); samples_by_down(sc).to_csv("t43_samples_by_down.csv")

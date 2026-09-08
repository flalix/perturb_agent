"""Part 3: sample table split into -2 down / -1 down / central / +1 up / +2 up, with genes per class.

Offsets beyond +-2 are collapsed into the +-2 classes.
Outputs: sample_classes.csv (gene lists), sample_class_counts.csv (gene counts)
"""
import numpy as np
import pandas as pd

CLASSES = {-2: "-2 down", -1: "-1 down", 0: "central", 1: "+1 up", 2: "+2 up"}

def build_sample_table(offset_matrix):
    off = offset_matrix.clip(-2, 2)
    lists = pd.DataFrame(index=off.index, columns=list(CLASSES.values()), dtype=object)
    counts = pd.DataFrame(0, index=off.index, columns=list(CLASSES.values()), dtype=int)
    for code, name in CLASSES.items():
        mask = off == code
        lists[name] = [";".join(off.columns[row]) for row in mask.to_numpy()]
        counts[name] = mask.sum(axis=1).to_numpy()
    counts["total_down"] = counts["-2 down"] + counts["-1 down"]
    counts["total_up"] = counts["+1 up"] + counts["+2 up"]
    return lists, counts

if __name__ == "__main__":
    offset_matrix = pd.read_csv("offset_matrix.csv", index_col=0)
    lists, counts = build_sample_table(offset_matrix)
    lists.to_csv("sample_classes.csv")
    counts.to_csv("sample_class_counts.csv")
    print(counts.head(10))

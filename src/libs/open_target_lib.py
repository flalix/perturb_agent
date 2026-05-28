#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/05/28
# Udated  on 2026/05/28
# @self.author: Flavio Lichtenstein

import os
from pathlib import Path
import duckdb
import pandas as pd
from typing import Optional, Iterable, Set, Tuple, Any, List

from IPython.display import display, HTML
# display(HTML("<style>.container { width:100% !important; }</style>"))
display(HTML("<style>:root_ot { --jp-notebook-max-width: 100% !important; }</style>"))

"""
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/evidence_cancer_gene_census .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/drug_mechanism_of_action .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/pharmacogenomics .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/study .
rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/26.03/output/evidence_reactome .

rsync -rpltvz --delete rsync.ebi.ac.uk::pub/databases/opentargets/platform/25.03/output/target .

colab/
    open_targets/
    ├── target/
    ├── disease/
    ├── drug_mechanism_of_action/
    ├── evidence_cancer_gene_census/
    ├── evidence_reactome/
    ├── interactions/
    ├── pharmacogenomics/
    └── study/

https://www.opentargets.org/
"""

from libs.Basic import create_dir, pdwritecsv, pdreadcsv

class OpenTarget(object):
    def __init__(self, root_colab:Path):

        self.root_colab = root_colab
        self.root_ot = root_colab / 'open_targets'

        self.con = duckdb.connect()

        self.tables = {
            "target": self.root_ot / "target",
            "disease": self.root_ot / "disease",
            "drug_moa": self.root_ot / "drug_mechanism_of_action",
            "cgc": self.root_ot / "evidence_cancer_gene_census",
            "reactome": self.root_ot / "evidence_reactome",
            "interactions": self.root_ot / "interactions",
            "pharmacogenomics": self.root_ot / "pharmacogenomics",
            "study": self.root_ot / "study",
            "association_by_datasource_direct": self.root_ot / "association_by_datasource_direct",
        }

        self._create_views()


    def is_installed(self):
        print("root exists:", self.root_ot.exists())
        print("root is dir:", self.root_ot.is_dir())

        print("\nTop-level content:")
        for p in sorted(self.root_ot.iterdir()):
            print(
                "DIR " if p.is_dir() else "FILE",
                p.name
            )


    def _glob(self, path: Path) -> str:
        return str(path / "**" / "*.parquet")

    def _create_views(self):
        for name, path in self.tables.items():
            if path.exists():
                self.con.execute(
                    f"""
                    CREATE OR REPLACE VIEW {name} AS
                    SELECT * FROM read_parquet('{self._glob(path)}')
                    """
                )
            else:
                print(f"Warning: missing table directory: {path}")

    def schema(self, table: str) -> pd.DataFrame:
        return self.con.execute(f"DESCRIBE {table}").df()
    
    # Resolve gene symbol or Ensembl ID
    def resolve_target(self, gene_or_ensembl: str) -> pd.DataFrame:
        symbol = gene_or_ensembl.strip()

        sql = """
        SELECT
            id AS target_id,
            approvedSymbol AS symbol,
            approvedName AS name,
            biotype
        FROM target
        WHERE
            id = ?
            OR upper(approvedSymbol) = upper(?)
        """

        return self.con.execute(sql, [symbol, symbol]).df()
    

    def show_schema_summary(self):

        for table in self.tables.keys():
            print("\n###", table)
            display(self.schema(table))


    def is_target_related_to_cancer(self, target_id: str) -> pd.DataFrame:
        df = self.con.execute(
            """
            SELECT
                c.*,
                d.name AS disease_name, d.synonyms, d.ontology
            FROM cgc c
            LEFT JOIN disease d
                ON c.diseaseId = d.id
            WHERE c.targetId = ?
            ORDER BY c.score DESC NULLS LAST
            LIMIT 50
            """, 
            [target_id]
        ).df()

        return df
    
    def is_target_related_to_disease(self, target_id: str, disease: str) -> pd.DataFrame:

        disease_like = f"%{disease.strip().lower()}%"

        df = self.con.execute(
            """
            SELECT
                a.targetId,
                t.approvedSymbol,
                t.approvedName,
                a.diseaseId,
                d.name AS disease_name,
                a.score,
                d.synonyms,
                d.ontology
            FROM association_by_datasource_direct a
            LEFT JOIN disease d
                ON a.diseaseId = d.id
            LEFT JOIN target t
                ON a.targetId = t.id
            WHERE a.targetId = ?
            AND (
                lower(coalesce(d.name, '')) LIKE ?
                OR lower(coalesce(d.description, '')) LIKE ?
                OR lower(coalesce(array_to_string(d.synonyms, ' '), '')) LIKE ?
            )
            ORDER BY a.score DESC NULLS LAST
            """,
            [target_id, disease_like, disease_like, disease_like]
        ).df()

        return df


    def get_drugs_for_disease(self, disease: str, limit: int = 200) -> pd.DataFrame:
        
        disease_like = f"%{disease.lower()}%"

        query = """
            SELECT
                d.id AS diseaseId,
                d.name AS disease_name,

                kd.drugId,
                kd.drugName,
                kd.targetId,
                t.approvedSymbol,
                t.approvedName,

                kd.phase,
                kd.status,
                kd.urls
            FROM known_drug kd
            LEFT JOIN disease d
                ON kd.diseaseId = d.id
            LEFT JOIN target t
                ON kd.targetId = t.id
            WHERE
                lower(coalesce(d.name, '')) LIKE ?
                OR lower(coalesce(d.description, '')) LIKE ?
                OR lower(coalesce(array_to_string(d.synonyms, ' '), '')) LIKE ?
            ORDER BY
                kd.phase DESC NULLS LAST,
                kd.drugName
            LIMIT ?
            """

        if isinstance(limit, int) and limit > 0:
            query += f" LIMIT {limit}"

        df = self.con.execute(
            query,
            [disease_like, disease_like, disease_like]
        ).df()

        return df


    def get_drugs_for_disease_with_moa(self, disease: str, limit: int = 300) -> pd.DataFrame:

        disease_like = f"%{disease.lower()}%"

        query =  """
            SELECT
                d.id AS diseaseId,
                d.name AS disease_name,

                kd.drugId,
                kd.drugName,
                kd.targetId,
                t.approvedSymbol,
                t.approvedName,

                kd.phase,
                kd.status,

                moa.mechanismOfAction,
                moa.actionType,
                moa.description AS moa_description

            FROM known_drug kd

            LEFT JOIN disease d 
                ON kd.diseaseId = d.id

            LEFT JOIN target t
                ON kd.targetId = t.id

            LEFT JOIN drug_moa moa
                ON kd.drugId = moa.drugId

            WHERE
                lower(coalesce(d.name, '')) LIKE ?
                OR lower(coalesce(d.description, '')) LIKE ?

            ORDER BY
                kd.phase DESC NULLS LAST,
                kd.drugName,
                t.approvedSymbol
            """
 
        if isinstance(limit, int) and limit > 0:
            query += f" LIMIT {limit}"
        
        df = self.con.execute(
            query,
            [disease_like, disease_like]
        ).df()

        return df


    def has_target_Reactome_evidence(self, target_id: str) -> pd.DataFrame:
        df = self.con.execute(
            """
            SELECT
                r.*,
                d.name AS disease_name
            FROM reactome r
            LEFT JOIN disease d
                ON r.diseaseId = d.id
            WHERE r.targetId = ?
            ORDER BY r.score DESC NULLS LAST
            LIMIT 50
            """,
            [target_id]
        ).df()

        return df
    

    def has_target_moa(self, target_id: str) -> pd.DataFrame:
        df = self.con.execute(
            """
            SELECT *
            FROM drug_moa
            WHERE CAST(targets AS VARCHAR) ILIKE ?
            """,
            [f"%{target_id}%"]
        ).df()

        return df
    
    def has_target_interactions(self, target_id: str, limit: int = 100) -> pd.DataFrame:

        query = \
            """
            SELECT *
            FROM interactions
            WHERE targetA = ?
            OR targetB = ?
            """
        
        if isinstance(limit, int) and limit > 0:
            query += f" LIMIT {limit}"

        df = self.con.execute(
            query,
            [target_id, target_id]
        ).df()

        return df
    

     
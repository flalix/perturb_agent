#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/26
# Udated  on 2026/03/26
# @author: Flavio Lichtenstein
# @local: Home sweet home

import os, sys
import pandas as pd

from kivy.app import App
from kivy.clock import Clock
from kivy.lang import Builder
from kivy.properties import ListProperty, ObjectProperty, StringProperty
from kivy.uix.boxlayout import BoxLayout

from pathlib import Path

ROOT = Path().resolve().parent.parent
SRC = os.path.join(ROOT, "src")

if str(SRC) not in sys.path:
    sys.path.append(str(SRC))

print("ROOT:", ROOT)
print("SRC added:", SRC)

from libs.tcga_gdc_lib import *
from libs.Basic import *


ROOT = Path().resolve().parent
root_data = os.path.join(ROOT, "data/tcga")

gdc = GDC(root_data=root_data)

verbose = True


KV = """
#:import dp kivy.metrics.dp

<PrimarySiteOption@SpinnerOption>:
    size_hint_y: None
    height: dp(36)
    background_normal: ""
    background_color: 0.18, 0.18, 0.22, 1
    color: 1, 1, 1, 1

<PrimarySiteDropDown@DropDown>:
    max_height: dp(36) * 3

<RowCell@Label>:
    size_hint_x: 1
    color: 0.95, 0.95, 0.95, 1
    text_size: self.size
    halign: "left"
    valign: "middle"
    padding: dp(6), 0

<HeaderCell@Label>:
    bold: True
    size_hint_x: 1
    color: 1, 1, 1, 1
    text_size: self.size
    halign: "left"
    valign: "middle"
    padding: dp(6), 0

<BoxHeader@BoxLayout>:
    size_hint_y: None
    height: dp(30)
    canvas.before:
        Color:
            rgba: 0.20, 0.24, 0.32, 1
        Rectangle:
            pos: self.pos
            size: self.size

<DataRow6>:
    orientation: "horizontal"
    size_hint_y: None
    height: dp(28)

    canvas.before:
        Color:
            rgba: root.bg_color
        Rectangle:
            pos: self.pos
            size: self.size

    RowCell:
        text: root.c1
    RowCell:
        text: root.c2
    RowCell:
        text: root.c3
    RowCell:
        text: root.c4
    RowCell:
        text: root.c5
    RowCell:
        text: root.c6

<RootWidget>:
    orientation: "vertical"
    padding: dp(8)
    spacing: dp(8)

    canvas.before:
        Color:
            rgba: 0.10, 0.10, 0.12, 1
        Rectangle:
            pos: self.pos
            size: self.size

    BoxLayout:
        size_hint_y: None
        height: dp(42)
        spacing: dp(8)

        Label:
            text: "Program"
            color: 1, 1, 1, 1
            size_hint_x: 0.18
            halign: "left"
            valign: "middle"
            text_size: self.size

        Spinner:
            id: sp_program
            text: root.program_selected if root.program_selected else "Select program"
            values: root.program_values
            background_normal: ""
            background_color: 0.18, 0.18, 0.22, 1
            color: 1, 1, 1, 1

    BoxLayout:
        size_hint_y: None
        height: dp(42)
        spacing: dp(8)

        Label:
            text: "Primary Site"
            color: 1, 1, 1, 1
            size_hint_x: 0.18
            halign: "left"
            valign: "middle"
            text_size: self.size

        Spinner:
            id: sp_primary_site
            text: root.primary_site_selected if root.primary_site_selected else "Select primary site"
            values: root.primary_site_values
            option_cls: "PrimarySiteOption"
            dropdown_cls: "PrimarySiteDropDown"
            background_normal: ""
            background_color: 0.18, 0.18, 0.22, 1
            color: 1, 1, 1, 1

    BoxLayout:
        size_hint_y: None
        height: dp(42)
        spacing: dp(8)

        Button:
            text: "Find cases, subtypes, tumor class and stages"
            background_normal: ""
            background_color: 0.24, 0.45, 0.78, 1
            color: 1, 1, 1, 1
            on_release: root.fetch_cases()

    Label:
        id: lbl_status
        text: root.status_text
        color: 0.85, 0.90, 1, 1
        size_hint_y: None
        height: dp(28)
        halign: "left"
        valign: "middle"
        text_size: self.size

    Label:
        text: "df_subt - subtype summary"
        bold: True
        color: 1, 1, 1, 1
        size_hint_y: None
        height: dp(28)
        halign: "left"
        valign: "middle"
        text_size: self.size

    BoxHeader:
        HeaderCell:
            text: "pid"
        HeaderCell:
            text: "subtype_global"
        HeaderCell:
            text: "tumor_class"
        HeaderCell:
            text: "subtype_tissue"
        HeaderCell:
            text: "sstage"
        HeaderCell:
            text: "n"

    RecycleView:
        id: rv_subt
        size_hint_y: 0.42
        viewclass: "DataRow6"
        data: root.subt_table_data
        scroll_type: ["bars", "content"]
        bar_width: dp(10)

        RecycleBoxLayout:
            default_size: None, dp(28)
            default_size_hint: 1, None
            size_hint_y: None
            height: self.minimum_height
            orientation: "vertical"

    Label:
        text: "df_cases2 - cases"
        bold: True
        color: 1, 1, 1, 1
        size_hint_y: None
        height: dp(28)
        halign: "left"
        valign: "middle"
        text_size: self.size

    BoxHeader:
        HeaderCell:
            text: "pid"
        HeaderCell:
            text: "case_id"
        HeaderCell:
            text: "subtype_global"
        HeaderCell:
            text: "tumor_class"
        HeaderCell:
            text: "subtype_tissue"
        HeaderCell:
            text: "stage"

    RecycleView:
        id: rv_cases
        size_hint_y: 0.42
        viewclass: "DataRow6"
        data: root.cases_table_data
        scroll_type: ["bars", "content"]
        bar_width: dp(10)

        RecycleBoxLayout:
            default_size: None, dp(28)
            default_size_hint: 1, None
            size_hint_y: None
            height: self.minimum_height
            orientation: "vertical"
"""


class DataRow6(BoxLayout):
    c1 = StringProperty("")
    c2 = StringProperty("")
    c3 = StringProperty("")
    c4 = StringProperty("")
    c5 = StringProperty("")
    c6 = StringProperty("")
    bg_color = ListProperty([0.16, 0.16, 0.18, 1])


class RootWidget(BoxLayout):
    verbose = True
    force = False

    program_values = ListProperty([])
    primary_site_values = ListProperty([])

    program_selected = StringProperty("TCGA")
    primary_site_selected = StringProperty("")

    status_text = StringProperty("Ready.")

    df_primary_sites = ObjectProperty(allownone=True)
    df_subt = ObjectProperty(allownone=True)
    df_cases2 = ObjectProperty(allownone=True)

    subt_table_data = ListProperty([])
    cases_table_data = ListProperty([])

    pid = StringProperty("")

    def __init__(self, **kwargs):
        self.df_primary_sites = pd.DataFrame()
        self.df_subt = pd.DataFrame()
        self.df_cases2 = pd.DataFrame()
        super().__init__(**kwargs)

    def on_kv_post(self, base_widget) -> None:
        self.ids.sp_program.bind(text=self._program_changed)
        self.ids.sp_primary_site.bind(text=self._primary_site_changed)
        Clock.schedule_once(lambda dt: self.load_programs(), 0)

    def load_programs(self) -> None:
        try:
            self.status_text = "Loading programs..."
            prog_list = gdc.get_gdc_progams(force=self.force, verbose=self.verbose)

            if not isinstance(prog_list, list):
                raise TypeError("gdc.get_gdc_progams() did not return a list")

            self.program_values = [str(x) for x in prog_list]

            if "TCGA" in self.program_values:
                self.program_selected = "TCGA"
            elif self.program_values:
                self.program_selected = self.program_values[0]
            else:
                self.program_selected = ""

            self.ids.sp_program.text = self.program_selected

            self.primary_site_values = []
            self.primary_site_selected = ""
            self.ids.sp_primary_site.text = "Select primary site"

            self.subt_table_data = []
            self.cases_table_data = []

            if self.program_selected:
                self.load_primary_sites(self.program_selected)

        except Exception as e:
            self.status_text = f"Error loading programs: {e}"

    def _program_changed(self, instance, value) -> None:
        self.program_selected = value

        self.primary_site_values = []
        self.primary_site_selected = ""
        self.ids.sp_primary_site.text = "Select primary site"

        self.subt_table_data = []
        self.cases_table_data = []

        self.df_subt = pd.DataFrame()
        self.df_cases2 = pd.DataFrame()
        self.pid = ""

        if value:
            self.load_primary_sites(value)

    def load_primary_sites(self, program: str) -> None:
        try:
            self.status_text = f"Loading primary sites for {program}..."

            dfc = gdc.get_primary_sites(
                program=program,
                force=False,
                verbose=self.verbose
            )

            if not isinstance(dfc, pd.DataFrame):
                raise TypeError("gdc.get_primary_sites() did not return a DataFrame")

            expected = {"pid", "primary_site", "project_id", "disease_type", "name"}
            missing = expected - set(dfc.columns)
            if missing:
                raise ValueError(f"Missing columns in dfc: {sorted(missing)}")

            self.df_primary_sites = dfc.copy()

            primary_sites = (
                self.df_primary_sites["primary_site"]
                .dropna()
                .astype(str)
                .sort_values()
                .unique()
                .tolist()
            )

            self.primary_site_values = primary_sites
            self.primary_site_selected = ""
            self.ids.sp_primary_site.text = "Select primary site"

            self.status_text = f"Loaded {len(primary_sites)} primary sites."

        except Exception as e:
            self.status_text = f"Error loading primary sites: {e}"
            self.primary_site_values = []
            self.primary_site_selected = ""
            self.ids.sp_primary_site.text = "Select primary site"

    def _primary_site_changed(self, instance, value) -> None:
        if value == "Select primary site":
            self.primary_site_selected = ""
        else:
            self.primary_site_selected = value

    def fetch_cases(self) -> None:
        try:
            if not self.program_selected:
                self.status_text = "Please select a program."
                return

            if not self.primary_site_selected:
                self.status_text = "Please select a primary site."
                return

            if self.df_primary_sites is None or self.df_primary_sites.empty:
                self.status_text = "Primary site table is empty."
                return

            df_match = self.df_primary_sites[
                self.df_primary_sites["primary_site"].astype(str) == str(self.primary_site_selected)
            ]

            if df_match.empty:
                self.status_text = "Could not find pid for the selected primary site."
                return

            self.pid = str(df_match.iloc[0]["pid"])
            self.status_text = f"Loading cases and subtypes for {self.pid}..."

            df_cases, df_subt, _ = gdc.get_cases_and_subtypes(
                pid=self.pid,
                batch_size=200,
                do_filter=False,
                force=False,
                verbose=self.verbose
            )

            if not isinstance(df_cases, pd.DataFrame):
                raise TypeError("df_cases is not a DataFrame")
            if not isinstance(df_subt, pd.DataFrame):
                raise TypeError("df_subt is not a DataFrame")

            subt_cols = ["pid", "subtype_global", "tumor_class", "subtype_tissue", "sstage", "n"]
            missing_subt = [c for c in subt_cols if c not in df_subt.columns]
            if missing_subt:
                raise ValueError(f"df_subt missing columns: {missing_subt}")

            case_cols = ["pid", "case_id", "subtype_global", "tumor_class", "subtype_tissue", "stage"]
            missing_cases = [c for c in case_cols if c not in df_cases.columns]
            if missing_cases:
                raise ValueError(f"df_cases missing columns: {missing_cases}")

            self.df_subt = df_subt[subt_cols].copy()
            self.df_cases2 = df_cases[case_cols].copy()

            self.subt_table_data = self._df_to_table_data(self.df_subt, subt_cols)
            self.cases_table_data = self._df_to_table_data(self.df_cases2, case_cols)

            self.status_text = (
                f"Loaded {len(self.df_subt)} subtype rows and "
                f"{len(self.df_cases2)} case rows for {self.pid}."
            )

        except Exception as e:
            self.status_text = f"Error loading cases/subtypes: {e}"

    @staticmethod
    def _safe_str(x) -> str:
        if pd.isna(x):
            return ""
        return str(x)

    def _df_to_table_data(self, df: pd.DataFrame, cols: list[str]) -> list[dict]:
        rows = []

        for i, record in enumerate(df[cols].to_dict("records")):
            bg = [0.16, 0.16, 0.18, 1] if i % 2 == 0 else [0.12, 0.12, 0.14, 1]

            rows.append(
                {
                    "c1": self._safe_str(record.get(cols[0], "")),
                    "c2": self._safe_str(record.get(cols[1], "")),
                    "c3": self._safe_str(record.get(cols[2], "")),
                    "c4": self._safe_str(record.get(cols[3], "")),
                    "c5": self._safe_str(record.get(cols[4], "")),
                    "c6": self._safe_str(record.get(cols[5], "")),
                    "bg_color": bg,
                }
            )

        return rows


class GDCExplorerApp(App):
    def build(self):
        Builder.load_string(KV)
        return RootWidget()


if __name__ == "__main__":
    GDCExplorerApp().run()
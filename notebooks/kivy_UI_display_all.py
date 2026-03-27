#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/03/26
# Udated  on 2026/03/26
# @author: Flavio Lichtenstein
# @local: Home sweet home

import pandas as pd

from kivy.app import App
from kivy.lang import Builder
from kivy.properties import ListProperty, ObjectProperty, StringProperty
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.recycleview.views import RecycleDataViewBehavior
from kivy.uix.widget import Widget

KV = """
#:import dp kivy.metrics.dp

<RowCell@Label>:
    size_hint_y: None
    height: dp(30)
    text_size: self.size
    halign: "left"
    valign: "middle"
    padding: dp(6), 0

<DataRow>:
    orientation: "horizontal"
    size_hint_y: None
    height: dp(30)

    RowCell:
        text: root.col1
    RowCell:
        text: root.col2
    RowCell:
        text: root.col3
    RowCell:
        text: root.col4
    RowCell:
        text: root.col5

<RootWidget>:
    orientation: "vertical"
    padding: dp(8)
    spacing: dp(8)

    BoxLayout:
        size_hint_y: None
        height: dp(42)
        spacing: dp(8)

        Spinner:
            id: sp_country
            text: root.country_selected
            values: root.country_values

        Spinner:
            id: sp_state
            text: root.state_selected
            values: root.state_values

        Spinner:
            id: sp_city
            text: root.city_selected
            values: root.city_values

    BoxLayout:
        orientation: "horizontal"
        size_hint_y: None
        height: dp(32)

        Label:
            text: "country"
            bold: True
        Label:
            text: "state"
            bold: True
        Label:
            text: "city"
            bold: True
        Label:
            text: "gene"
            bold: True
        Label:
            text: "value"
            bold: True

    RecycleView:
        id: rv
        viewclass: "DataRow"
        data: root.table_data

        RecycleBoxLayout:
            default_size: None, dp(30)
            default_size_hint: 1, None
            size_hint_y: None
            height: self.minimum_height
            orientation: "vertical"

    BoxLayout:
        size_hint_y: None
        height: dp(42)
        spacing: dp(8)

        TextInput:
            id: gene_filter
            hint_text: "Filter gene contains..."
            multiline: False

        TextInput:
            id: min_value
            hint_text: "Min value"
            multiline: False
            input_filter: "float"

        TextInput:
            id: max_value
            hint_text: "Max value"
            multiline: False
            input_filter: "float"
"""


class DataRow(RecycleDataViewBehavior, BoxLayout):
    col1 = StringProperty("")
    col2 = StringProperty("")
    col3 = StringProperty("")
    col4 = StringProperty("")
    col5 = StringProperty("")


class RootWidget(BoxLayout):
    df = ObjectProperty(None)

    country_values = ListProperty(["All"])
    state_values = ListProperty(["All"])
    city_values = ListProperty(["All"])

    country_selected = StringProperty("All")
    state_selected = StringProperty("All")
    city_selected = StringProperty("All")

    table_data = ListProperty([])

    def __init__(self, **kwargs):
        self.df = pd.DataFrame(
            {
                "country": ["Brazil", "Brazil", "Brazil", "USA", "USA", "Germany", "Germany"],
                "state": ["SP", "SP", "RJ", "CA", "NY", "BW", "BY"],
                "city": ["Santo Andre", "Sao Paulo", "Rio", "Los Angeles", "New York", "Heidelberg", "Munich"],
                "gene": ["TP53", "EGFR", "BRCA1", "KRAS", "PIK3CA", "MYC", "PTEN"],
                "value": [10, 20, 15, 30, 12, 18, 25],
            }
        )
        super().__init__(**kwargs)

    def on_kv_post(self, base_widget: Widget) -> None:
        self.country_values = ["All"] + sorted(self.df["country"].dropna().unique().tolist())
        self.state_values = ["All"] + sorted(self.df["state"].dropna().unique().tolist())
        self.city_values = ["All"] + sorted(self.df["city"].dropna().unique().tolist())

        self.ids.sp_country.bind(text=self._country_changed)
        self.ids.sp_state.bind(text=self._state_changed)
        self.ids.sp_city.bind(text=self._city_changed)

        self.ids.gene_filter.bind(text=lambda *_: self.apply_filters())
        self.ids.min_value.bind(text=lambda *_: self.apply_filters())
        self.ids.max_value.bind(text=lambda *_: self.apply_filters())

        self.apply_filters()

    def _country_changed(self, instance, value):
        self.country_selected = value
        self.state_selected = "All"
        self.city_selected = "All"

        df = self.df.copy()
        if value != "All":
            df = df[df["country"] == value]

        self.state_values = ["All"] + sorted(df["state"].dropna().unique().tolist())
        self.city_values = ["All"] + sorted(df["city"].dropna().unique().tolist())

        self.ids.sp_state.text = "All"
        self.ids.sp_city.text = "All"

        self.apply_filters()

    def _state_changed(self, instance, value):
        self.state_selected = value
        self.city_selected = "All"

        df = self.df.copy()
        if self.country_selected != "All":
            df = df[df["country"] == self.country_selected]
        if value != "All":
            df = df[df["state"] == value]

        self.city_values = ["All"] + sorted(df["city"].dropna().unique().tolist())
        self.ids.sp_city.text = "All"

        self.apply_filters()

    def _city_changed(self, instance, value):
        self.city_selected = value
        self.apply_filters()

    def apply_filters(self):
        df = self.df.copy()

        if self.country_selected != "All":
            df = df[df["country"] == self.country_selected]
        if self.state_selected != "All":
            df = df[df["state"] == self.state_selected]
        if self.city_selected != "All":
            df = df[df["city"] == self.city_selected]

        gene_txt = self.ids.gene_filter.text.strip()
        if gene_txt:
            df = df[df["gene"].str.contains(gene_txt, case=False, na=False)]

        min_txt = self.ids.min_value.text.strip()
        if min_txt:
            try:
                df = df[df["value"] >= float(min_txt)]
            except ValueError:
                pass

        max_txt = self.ids.max_value.text.strip()
        if max_txt:
            try:
                df = df[df["value"] <= float(max_txt)]
            except ValueError:
                pass

        self.table_data = [
            {
                "col1": str(row["country"]),
                "col2": str(row["state"]),
                "col3": str(row["city"]),
                "col4": str(row["gene"]),
                "col5": str(row["value"]),
            }
            for _, row in df.iterrows()
        ]

class ChainedDropdownTableApp(App):
    def build(self):
        Builder.load_string(KV)
        return RootWidget()


if __name__ == "__main__":
    ChainedDropdownTableApp().run()
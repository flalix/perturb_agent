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
from kivy.metrics import dp
from kivy.properties import ListProperty, StringProperty
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.button import Button
from kivy.uix.dropdown import DropDown
from kivy.uix.recycleview import RecycleView
from kivy.uix.recycleview.views import RecycleDataViewBehavior


KV = """
<RowWidget>:
    orientation: "horizontal"
    size_hint_y: None
    height: dp(32)
    canvas.before:
        Color:
            rgba: (0.95, 0.95, 0.95, 1) if self.index % 2 == 0 else (1, 1, 1, 1)
        Rectangle:
            pos: self.pos
            size: self.size

    Label:
        text: root.country
        size_hint_x: 0.18
    Label:
        text: root.state
        size_hint_x: 0.14
    Label:
        text: root.city
        size_hint_x: 0.20
    Label:
        text: root.gene
        size_hint_x: 0.20
    Label:
        text: root.value
        size_hint_x: 0.12
    Label:
        text: root.category
        size_hint_x: 0.16

<ResultsTable>:
    viewclass: "RowWidget"
    RecycleBoxLayout:
        default_size: None, dp(32)
        default_size_hint: 1, None
        size_hint_y: None
        height: self.minimum_height
        orientation: "vertical"

<MainWidget>:
    orientation: "vertical"
    padding: dp(10)
    spacing: dp(10)

    BoxLayout:
        size_hint_y: None
        height: dp(42)
        spacing: dp(8)

        Button:
            id: btn_country
            text: root.selected_country
            on_release: root.country_dropdown.open(self)

        Button:
            id: btn_state
            text: root.selected_state
            on_release: root.state_dropdown.open(self)

        Button:
            id: btn_city
            text: root.selected_city
            on_release: root.city_dropdown.open(self)

    BoxLayout:
        orientation: "horizontal"
        size_hint_y: None
        height: dp(32)
        spacing: dp(2)
        canvas.before:
            Color:
                rgba: 0.75, 0.75, 0.75, 1
            Rectangle:
                pos: self.pos
                size: self.size

        Label:
            text: "Country"
            bold: True
            size_hint_x: 0.18
        Label:
            text: "State"
            bold: True
            size_hint_x: 0.14
        Label:
            text: "City"
            bold: True
            size_hint_x: 0.20
        Label:
            text: "Gene"
            bold: True
            size_hint_x: 0.20
        Label:
            text: "Value"
            bold: True
            size_hint_x: 0.12
        Label:
            text: "Category"
            bold: True
            size_hint_x: 0.16

    ResultsTable:
        id: rv

    BoxLayout:
        orientation: "horizontal"
        size_hint_y: None
        height: dp(42)
        spacing: dp(8)

        TextInput:
            id: gene_filter
            hint_text: "Filter gene contains..."
            multiline: False
            on_text: root.apply_filters()

        TextInput:
            id: min_value
            hint_text: "Min value"
            multiline: False
            input_filter: "float"
            on_text: root.apply_filters()

        TextInput:
            id: max_value
            hint_text: "Max value"
            multiline: False
            input_filter: "float"
            on_text: root.apply_filters()

        Button:
            text: "Reset"
            size_hint_x: None
            width: dp(100)
            on_release: root.reset_filters()
"""
class RowWidget(RecycleDataViewBehavior, BoxLayout):
    index = 0
    country = StringProperty("")
    state = StringProperty("")
    city = StringProperty("")
    gene = StringProperty("")
    value = StringProperty("")
    category = StringProperty("")


class ResultsTable(RecycleView):
    pass


class MainWidget(BoxLayout):
    selected_country = StringProperty("Country: All")
    selected_state = StringProperty("State: All")
    selected_city = StringProperty("City: All")

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.df = pd.DataFrame({
            "country": ["Brazil", "Brazil", "Brazil", "USA", "USA", "Germany", "Germany", "Brazil"],
            "state":   ["SP", "SP", "RJ", "CA", "NY", "BW", "BY", "SP"],
            "city":    ["Santo Andre", "Sao Paulo", "Rio", "Los Angeles", "New York", "Heidelberg", "Munich", "Campinas"],
            "gene":    ["TP53", "EGFR", "BRCA1", "KRAS", "PIK3CA", "MYC", "PTEN", "EGFR"],
            "value":   [10, 20, 15, 30, 12, 18, 25, 22],
            "category": ["A", "B", "A", "B", "C", "A", "B", "C"]
        })

        self.country_value = "All"
        self.state_value = "All"
        self.city_value = "All"

        self.country_dropdown = DropDown()
        self.state_dropdown = DropDown()
        self.city_dropdown = DropDown()

        self.build_country_dropdown()
        self.build_state_dropdown()
        self.build_city_dropdown()
        self.apply_filters()

    def make_button(self, text, callback):
        btn = Button(text=text, size_hint_y=None, height=dp(36))
        btn.bind(on_release=callback)
        return btn

    def build_country_dropdown(self):
        self.country_dropdown.clear_widgets()
        values = ["All"] + sorted(self.df["country"].dropna().unique().tolist())
        for value in values:
            btn = self.make_button(value, lambda btn, v=value: self.set_country(v))
            self.country_dropdown.add_widget(btn)

    def build_state_dropdown(self):
        self.state_dropdown.clear_widgets()
        df = self.df.copy()
        if self.country_value != "All":
            df = df[df["country"] == self.country_value]

        values = ["All"] + sorted(df["state"].dropna().unique().tolist())
        for value in values:
            btn = self.make_button(value, lambda btn, v=value: self.set_state(v))
            self.state_dropdown.add_widget(btn)

    def build_city_dropdown(self):
        self.city_dropdown.clear_widgets()
        df = self.df.copy()

        if self.country_value != "All":
            df = df[df["country"] == self.country_value]
        if self.state_value != "All":
            df = df[df["state"] == self.state_value]

        values = ["All"] + sorted(df["city"].dropna().unique().tolist())
        for value in values:
            btn = self.make_button(value, lambda btn, v=value: self.set_city(v))
            self.city_dropdown.add_widget(btn)

    def set_country(self, value):
        self.country_value = value
        self.state_value = "All"
        self.city_value = "All"

        self.selected_country = f"Country: {value}"
        self.selected_state = "State: All"
        self.selected_city = "City: All"

        self.country_dropdown.dismiss()
        self.build_state_dropdown()
        self.build_city_dropdown()
        self.apply_filters()

    def set_state(self, value):
        self.state_value = value
        self.city_value = "All"

        self.selected_state = f"State: {value}"
        self.selected_city = "City: All"

        self.state_dropdown.dismiss()
        self.build_city_dropdown()
        self.apply_filters()

    def set_city(self, value):
        self.city_value = value
        self.selected_city = f"City: {value}"
        self.city_dropdown.dismiss()
        self.apply_filters()

    def get_filtered_df(self):
        df = self.df.copy()

        if self.country_value != "All":
            df = df[df["country"] == self.country_value]
        if self.state_value != "All":
            df = df[df["state"] == self.state_value]
        if self.city_value != "All":
            df = df[df["city"] == self.city_value]

        gene_text = self.ids.gene_filter.text.strip().lower()
        if gene_text:
            df = df[df["gene"].str.lower().str.contains(gene_text, na=False)]

        min_text = self.ids.min_value.text.strip()
        max_text = self.ids.max_value.text.strip()

        if min_text:
            try:
                df = df[df["value"] >= float(min_text)]
            except ValueError:
                pass

        if max_text:
            try:
                df = df[df["value"] <= float(max_text)]
            except ValueError:
                pass

        return df

    def apply_filters(self):
        df = self.get_filtered_df()

        self.ids.rv.data = [
            {
                "index": i,
                "country": str(row["country"]),
                "state": str(row["state"]),
                "city": str(row["city"]),
                "gene": str(row["gene"]),
                "value": str(row["value"]),
                "category": str(row["category"]),
            }
            for i, (_, row) in enumerate(df.iterrows())
        ]

    def reset_filters(self):
        self.country_value = "All"
        self.state_value = "All"
        self.city_value = "All"

        self.selected_country = "Country: All"
        self.selected_state = "State: All"
        self.selected_city = "City: All"

        self.ids.gene_filter.text = ""
        self.ids.min_value.text = ""
        self.ids.max_value.text = ""

        self.build_state_dropdown()
        self.build_city_dropdown()
        self.apply_filters()


class KivyFilterApp(App):
    def build(self):
        Builder.load_string(KV)
        return MainWidget()


if __name__ == "__main__":
    KivyFilterApp().run()
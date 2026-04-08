#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/04/07
# Updated  on 2026/04/07
# @author: Flavio Lichtenstein
# @local: Home sweet home

from pathlib import Path
import runpy

target = Path(__file__).parent / "notebooks" / "streamlit_UI_GDC_test03.py"
runpy.run_path(str(target), run_name="__main__")
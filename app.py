#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/04/07
# Updated  on 2026/04/07
# @author: Flavio Lichtenstein
# @local: Home sweet home

from pathlib import Path
import runpy

target = Path(__file__).resolve().parent / "src" / "app_main.py"
print(f"Running {target}")
runpy.run_path(str(target), run_name="__main__")
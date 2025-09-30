 #!/usr/bin/python3
# astrolabe.py
# -*- coding: utf-8 -*-
#
# The python script in this file makes the various parts of a model astrolabe.
#
# Copyright (C) 2010-2024 Dominic Ford <https://dcford.org.uk/>
#
# This code is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# You should have received a copy of the GNU General Public License along with
# this file; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA  02110-1301, USA

# ----------------------------------------------------------------------------

"""
This is the top level script for drawing all the parts needed to build astrolabes which work at a range of
different latitudes. They are rendered in PDF, SVG and PNG image formats.

Additionally, we use LaTeX to build a summary document for each latitude, which includes all of the parts needed
to build an astrolabe for that latitude, and instructions as to how to put them together.
"""

import os
import subprocess
import time
from typing import Dict, Union

import text
from climate import Climate
from graphics_context import GraphicsPage, CompositeComponent
from mother_back import MotherBack
from mother_front import MotherFront
from rete import Rete
from rule import Rule
from settings import fetch_command_line_arguments


#arguments: Dict[str, Union[int, str]] = fetch_command_line_arguments()
#theme: str = arguments['theme']

# Render astrolabe in all available languages
language: str = "fr"
format="svg"
theme="default"

places = ((43.584533, 7.115050, "Antibes"),
         (41.823150, 8.787917, "Acellasca"),
         )

# Render simplified and full astrolabes
astrolabe_type: str
for astrolabe_type in ["full"]: #, "simplified"

    for latitude, longitude, placename in places[0:1]:
        
        # Boolean flag for which hemisphere we're in
        southern: bool = latitude < 0

        # A dictionary of common substitutions
        subs: Dict[str, Union[str, float]] = {
            "dir_parts": ".",
            "abs_lat": int(abs(latitude)),
            "ns": "S" if southern else "N",
            "astrolabe_type": astrolabe_type,
            "lang": language,
            "lang_short": "" if language == "en" else "_{}".format(language)
        }

        settings: Dict[str, Union[str, float]] = {
            'language': language,
            "astrolabe_type": astrolabe_type,
            'latitude': latitude,
            'placename': placename,
            'theme': theme
        }

        # Render the parts of the astrolabe that do not change with geographic location
        if False:
            MotherFront(settings=settings).render_to_file(
                filename="{dir_parts}/mother_front_{abs_lat:02d}{ns}_{lang}_{astrolabe_type}".format(**subs),
                img_format=format
            )
        
        if True:
            MotherBack(settings=settings).render_to_file(
                filename="{dir_parts}/mother_back_{abs_lat:02d}{ns}_{lang}_{astrolabe_type}".format(**subs),
                img_format=format,
            )
        #
        if False:
            Rete(settings=settings).render_to_file(
                filename="{dir_parts}/rete_{abs_lat:02d}{ns}_{lang}_{astrolabe_type}".format(**subs),
                img_format=format
            )
        
        if False:
            Rule(settings=settings).render_to_file(
                filename="{dir_parts}/rule_{abs_lat:02d}{ns}_{lang}_{astrolabe_type}".format(**subs),
                img_format=format
            )

        # Render the climate of the astrolabe
        if False:
            Climate(settings=settings).render_to_file(
                filename="{dir_parts}/climate_{abs_lat:02d}{ns}_{lang}_{astrolabe_type}".format(**subs),
                img_format=format,
            )

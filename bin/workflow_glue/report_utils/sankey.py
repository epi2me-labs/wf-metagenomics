#!/usr/bin/env python
"""Interactive Sankey plot."""
import json
import os

import dominate.tags as html_tags
from dominate.util import raw as html_raw


sankey_stuff = os.path.dirname(os.path.realpath(__file__))


def sankey_plot(all_json, report):
    """Create a sankey plot from lineage json files.

    :param all_json (dict): Dictionary with all the info from lineages json files.
    :param report (html): Report template.
    :return (html): report template with the sankey.
    """
    scripts = [
        "https://d3js.org/d3.v6.min.js",
        "https://cdn.jsdelivr.net/npm/d3-color@3",
        "https://cdn.jsdelivr.net/npm/d3-dispatch@3",
        "https://cdn.jsdelivr.net/npm/d3-ease@3",
        "https://cdn.jsdelivr.net/npm/d3-interpolate@3",
        "https://cdn.jsdelivr.net/npm/d3-selection@3",
        "https://cdn.jsdelivr.net/npm/d3-timer@3",
        "https://cdn.jsdelivr.net/npm/d3-transition@3",
        "https://cdn.jsdelivr.net/npm/d3-drag@3",
        "https://cdn.jsdelivr.net/npm/d3-zoom@3",
        "https://unpkg.com/d3-array@1",
        "https://unpkg.com/d3-collection@1",
        "https://unpkg.com/d3-path@1",
        "https://unpkg.com/d3-shape@1",
        "https://unpkg.com/d3-sankey@0",
    ]
    with report.header:
        for script in scripts:
            html_tags.script(src=script)

    # Load sankey code and replace data
    replace_me = json.dumps(all_json).replace('"', '\"')
    with open(f'{sankey_stuff}/sankey/sankey.js') as sankey_code_js:
        sankey_js = sankey_code_js.read()
    sankey_data = sankey_js.replace(
            "replace_me",
            replace_me.replace('"', '\\"'))

    with report.add_section('Lineage plot', 'Lineage'):
        # Sankey div
        with html_tags.div(className="container"):
            # Load sankey style
            with open(f'{sankey_stuff}/sankey/sankey_style.css') as sankey_style_css:
                sankey_css = sankey_style_css.read()
            html_tags.style(html_raw(sankey_css))
            with html_tags.main(className="wf-metagenomics-sankey"):
                html_tags.div(html_raw(
                    """
                    <ul>
                        <li>
                            <label for="sample-select">Select sample</label>
                            <select id="sample-select"></select>
                        </li>
                        <li>
                            <label for="rank-select">Select rank</label>
                            <select id="rank-select"></select>
                        </li>
                        <li>
                            <label for="cutoff-select">Select cutoff</label>
                            <select id="cutoff-select"></select>
                        </li>
                        <li>
                            <button onclick="zoomIn()">[+] Zoom in</button>
                        </li>
                        <li>
                            <button onclick="zoomOut()">[-] Zoom out</button>
                        </li>
                        <li>
                            <button onclick="resetZoom()">[ ] Reset zoom</button>
                        </li>
                    </ul>
                    """), id="controls")
            with html_tags.div(id="visualisation"):
                html_tags.div(id="sankey-plot")
                html_tags.div(id="tooltip")
                html_tags.script(html_raw(sankey_data))

    return report

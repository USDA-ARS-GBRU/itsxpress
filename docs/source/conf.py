# -*- coding: utf-8 -*-
import os
import sys
import subprocess

# Go up two levels from docs/source/ to find the package root folder
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
project = 'ITSxpress'
copyright = 'CC0 Public Domain Attribution 2018-2026'
author = 'Adam R. Rivers'

# Dynamic version parsing via Git tag history tracking
try:
    git_version = subprocess.check_output(
        ["git", "describe", "--tags", "--always"], stderr=subprocess.STDOUT
    ).strip().decode("utf-8")
    if git_version.startswith("v"):
        git_version = git_version[1:]
except Exception:
    git_version = "2.1.4"  # Modern baseline fallback matching current release context

release = git_version
version = ".".join(git_version.split(".")[:2]) if "." in git_version else git_version

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
source_suffix = '.rst'
master_doc = 'index'
language = 'en'

# -- API Reference Extension Control -----------------------------------------
autodoc_default_options = {
    'special-members': '__init__',
    'private-members': None,
}

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Automatic API Documentation Generation Hook ----------------------------
def run_apidoc(app):
    """Automatically tracks layout modifications on documentation initialization."""
    from sphinx.ext import apidoc
    
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    output_path = os.path.join(cur_dir, "api")
    module_path = os.path.join(cur_dir, "../..", "itsxpress")
    exclude_paths = [os.path.join(cur_dir, "../..", "recipes")]
    
    apidoc.main(["-M", "-f", "-e", "-o", output_path, module_path] + exclude_paths)

def setup(app):
    """Hooks dynamic compilation tracking pipelines into the Sphinx runner."""
    app.connect("builder-inited", run_apidoc)

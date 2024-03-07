import sys
import os

sys.path.insert(0, os.path.abspath("../.."))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Picked Group FDR"
copyright = "2024, Matthew The"
author = "Matthew The"
release = "0.6.6"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
]

autodoc_mock_imports = [
    "pandas",
    "numpy",
    "scipy",
    "networkx",
    "matplotlib",
    "triqler",
    "cython",
    "llvmlite",
    "mokapot",
    "Bottleneck",
    "toml",
    "pyqt5",
    "job-pool",
]
autosummary_generate = True
autosummary_ignore_module_all = False
autosummary_imported_members = True

autoclass_content = "both"

napoleon_custom_sections = [("Required columns", "returns_style")]

templates_path = ["_templates"]
exclude_patterns = []

toc_object_entries = False


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# -- Options for HTML output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_css_files = [
    "custom.css",
]

pygments_style = "sphinx"

# -- Options for EPUB output
epub_show_urls = "footnote"

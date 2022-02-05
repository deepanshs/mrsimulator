# -*- coding: utf-8 -*-
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config
# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import datetime
import os
import subprocess
import sys
import warnings

from sphinx_gallery.sorting import ExplicitOrder
from sphinx_gallery.sorting import FileNameSortKey

# import plotly.io as pio
# from plotly.io import _sg_scraper

sys.path.insert(0, os.path.abspath("../.."))

# scraper = _sg_scraper.plotly_sg_scraper
# -- Project information -----------------------------------------------------
now = datetime.datetime.now()
year = now.year

project = "mrsimulator"
copyright = f"2019-{year}, The Mrsimulator Developers"
author = "The Mrsimulator Developers"


# get version number from the file
with open("../src/mrsimulator/__init__.py", "r") as f:
    for line in f.readlines():
        if "__version__" in line:
            before_keyword, keyword, after_keyword = line.partition("=")
            __version__ = after_keyword.strip()[1:-1]
            # Must break since want only first use of __version__ within file
            break

# The short X.Y version
version = ".".join(__version__.split(".")[:2])
# The full version, including alpha/beta/rc tags
release = __version__

# -- General configuration ---------------------------------------------------
show_authors = True
# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = "2.0"

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "matplotlib.sphinxext.plot_directive",
    "sphinx.ext.mathjax",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    # "sphinxcontrib.bibtex",
    "breathe",
    "sphinxjp.themes.basicstrap",
    "sphinx_gallery.gen_gallery",
    "sphinx.ext.intersphinx",
    "sphinx_tabs.tabs",
    "sphinx.ext.todo",
    "recommonmark",
    "versionwarning.extension",
]

# generate autosummary even if no references
autosummary_generate = True

# ---------------------------------------------------------------------------- #
#                              Sphinx Version warning                          #
# ---------------------------------------------------------------------------- #
# sphinx-version-warning config
versionwarning_messages = {
    "latest": (
        "This document is for the development version. "
        'For the stable version documentation, see <a href="/en/stable/">here</a>.'
    )
}
# Show warning at top of page
versionwarning_body_selector = "div.document"
# versionwarning_banner_title = ""
# For debugging locally
# versionwarning_project_version = "latest"

# ---------------------------------------------------------------------------- #
#                                  Plotly config                               #
# ---------------------------------------------------------------------------- #
# pio.renderers.default = "sphinx_gallery_png"

# ---------------------------------------------------------------------------- #
#                               Plot directive config                          #
# ---------------------------------------------------------------------------- #
plot_html_show_source_link = False
plot_rcparams = {
    "figure.figsize": [6, 3],
}

# ---------------------------------------------------------------------------- #
#                               Sphinx Gallery config                          #
# ---------------------------------------------------------------------------- #
# filter sphinx matplotlib warning
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message="Matplotlib is currently using agg, which is a"
    " non-GUI backend, so cannot show the figure.",
)
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=(
        "The physical quantity name, 'plane angle', is not "
        "defined in the astropy.units package. Continuing "
        "with 'plane angle' as the physical quantity name "
        "for unit deg."
    ),
)
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    message=(
        "The physical quantity name, 'electric field strength', is not "
        "defined in the astropy.units package. Continuing "
        "with 'electric field strength' as the physical quantity name "
        "for unit N / C."
    ),
)


# sphinx gallery config
sphinx_gallery_conf = {
    "examples_dirs": [
        "../signal_processing_source",
        "../examples_source",
        "../fitting_source",
    ],
    "remove_config_comments": True,
    "gallery_dirs": [
        "signal_processing",
        "examples",
        "fitting",
    ],  # path to where to save gallery generated output
    "within_subsection_order": FileNameSortKey,
    # "show_memory": True,
    # "line_numbers": True,
    "subsection_order": ExplicitOrder(
        [
            "../signal_processing_source",
            "../examples_source/1D_simulation(crystalline)",
            "../examples_source/1D_simulation(macro_amorphous)",
            "../examples_source/2D_simulation(crystalline)",
            "../examples_source/2D_simulation(macro_amorphous)",
            "../fitting_source/1D_fitting",
            "../fitting_source/2D_fitting",
        ]
    ),
    "reference_url": {
        # The module you locally document uses None
        "mrsimulator": None
    },
    "backreferences_dir": "examples",
    "doc_module": ("mrsimulator"),
    "image_scrapers": ["matplotlib"],  # , scraper],
    "compress_images": ("images", "thumbnails"),
    "capture_repr": ("_repr_html_", "__repr__"),
    # --- uncomment to include binder links for notebooks ---
    # "first_notebook_cell": (
    #     "# This cell is added by sphinx-gallery\n"
    #     "!pip install mrsimulator --quiet\n\n\n"
    #     "%matplotlib inline\n\n"
    #     "import mrsimulator\n"
    #     "print(f'You are using mrsimulator v{mrsimulator.__version__}')"
    # ),
    # "binder": {
    #     # Required keys
    #     "org": "deepanshs",
    #     "repo": "mrsimulator",
    #     "branch": "master",
    #     "binderhub_url": "https://mybinder.org",
    #     "dependencies": "../requirements.txt",
    #     # Optional keys
    #     "filepath_prefix": "docs/_build/html",
    #     "notebooks_dir": "../../notebooks",
    #     "use_jupyter_lab": True,
    # },
}

intersphinx_mapping = {
    "matplotlib": ("https://matplotlib.org", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference/", None),
    "csdmpy": ("https://csdmpy.readthedocs.io/en/stable/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "lmfit": ("https://lmfit-py.readthedocs.io/en/stable/", None),
}

# ---------------------------------------------------------------------------- #
#                              Sphinx copybutton                               #
# ---------------------------------------------------------------------------- #
copybutton_prompt_text = ">>> |\\$ |\\[\\d*\\]: |\\.\\.\\.: |[.][.][.] "
copybutton_prompt_is_regexp = True

# ---------------------------------------------------------------------------- #
#                            Doxygen C docs config                             #
# ---------------------------------------------------------------------------- #
try:
    subprocess.run("doxygen", shell=False)
    doxy_output = os.path.abspath("./xml")

    # Setup the breathe extension
    breathe_projects = {"My Project": doxy_output}
    breathe_default_project = "My Project"
    breathe_domain_by_extension = {"h": "c", "py": "py"}
    breathe_use_project_refids = True
    breathe_doxygen_config_options = {
        "PREDEFINED": "DOXYGEN_SHOULD_SKIP_THIS",
        # "GENERATE_XML": True,
        # "XML_PROGRAMLISTING": True,
        # "INPUT": "../../src/c_lib/include",
    }
except Exception:
    print("Skipping C-docs (doxygen)....")
# ---------------------------------------------------------------------------- #

# numfig config
numfig = True
numfig_secnum_depth = 2
numfig_format = {"figure": "Figure %s", "table": "Table %s", "code-block": "Listing %s"}

# math
math_number_all = True

# Tell sphinx what the primary language being documented is.
primary_domain = "py"

# Tell sphinx what the pygments highlight language should be.
highlight_language = "c"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}

# The master toctree document.
master_doc = "index"

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = None


# ---------------------------------------------------------------------------- #
#                                  HTML theme                                  #
# ---------------------------------------------------------------------------- #

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# Some html_theme options are 'alabaster', 'bootstrap', 'sphinx_rtd_theme',
# 'classic', 'basicstrap'
html_theme = "basicstrap"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    # Set the lang attribute of the html tag. Defaults to 'en'
    "lang": "en",
    # Disable showing the sidebar. Defaults to 'false'
    "nosidebar": False,
    # Show header searchbox. Defaults to false. works only "nosidebar=True",
    "header_searchbox": False,
    # Put the sidebar on the right side. Defaults to false.
    "rightsidebar": False,
    # Set the width of the sidebar. Defaults to 3
    "sidebar_span": 3,
    # Fix navbar to top of screen. Defaults to true
    "nav_fixed_top": True,
    # Fix the width of the sidebar. Defaults to false
    "nav_fixed": True,
    # Set the width of the sidebar. Defaults to '900px'
    "nav_width": "300px",
    # Fix the width of the content area. Defaults to false
    "content_fixed": False,
    # Set the width of the content area. Defaults to '900px'
    "content_width": "900px",
    # Fix the width of the row. Defaults to false
    "row_fixed": False,
    # Disable the responsive design. Defaults to false
    "noresponsive": False,
    # Disable the responsive footer relbar. Defaults to false
    "noresponsiverelbar": False,
    # Disable flat design. Defaults to false.
    # Works only "bootstrap_version = 3"
    "noflatdesign": False,
    # Enable Google Web Font. Defaults to false
    # "googlewebfont": True,
    # Set the URL of Google Web Font's CSS.
    # Defaults to 'http://fonts.googleapis.com/css?family=Text+Me+One'
    # "googlewebfont_url": "http://fonts.googleapis.com/css?family=Roboto+Script+One",  # NOQA
    # "googlewebfont_url": "http://fonts.googleapis.com/css2?family=Inter",
    # Set the Style of Google Web Font's CSS.
    # Defaults to "font-family: 'Text Me One', sans-serif;"
    # "googlewebfont_style": "font-family: Helvetica",
    # "googlewebfont_style": "font-family: 'Inter', sans-serif;",
    # Set 'navbar-inverse' attribute to header navbar. Defaults to false.
    "header_inverse": True,
    # Set 'navbar-inverse' attribute to relbar navbar. Defaults to false.
    "relbar_inverse": False,
    # Enable inner theme by Bootswatch. Defaults to false
    "inner_theme": False,
    # Set the name of inner theme. Defaults to 'bootswatch-simplex'
    # "inner_theme_name": "bootswatch-Yeti",
    # Select Twitter bootstrap version 2 or 3. Defaults to '3'
    "bootstrap_version": "3",
    # Show "theme preview" button in header navbar. Defaults to false.
    "theme_preview": False,
}

html_style = "style.css"
html_title = f"mrsimulator:docs v{__version__}"
html_logo = "_static/mrsimulator_logo.png"
html_favicon = "_static/favicon.ico"
html_last_updated_fmt = ""

html_sidebars = {
    "**": ["searchbox.html", "globaltoc.html"],
    "using/windows": ["searchbox.html", "windowssidebar.html"],
}
# ---------------------------------------------------------------------------- #

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "MRSimulatordoc"


# ---------------------------------------------------------------------------- #
#                              LaTeX setup                                     #
# ---------------------------------------------------------------------------- #
latex_engine = "pdflatex"
latex_logo = "_static/mrsimulator_logo.pdf"
latex_show_pagerefs = True
latex_toplevel_sectioning = "part"
latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    "papersize": "letterpaper",
    # The font size ('10pt', '11pt' or '12pt').
    #
    "pointsize": "10pt",
    "fontenc": "\\usepackage[utf8]{inputenc}",
    "fontpkg": "\\usepackage{amsmath,amsfonts,amssymb,amsthm}",
    # "geometry": "\\usepackage[vmargin=2.5cm, hmargin=1.5cm]{geometry}",
    # "fncychap": "\usepackage[Rejne]{fncychap}",
    # Additional stuff for the LaTeX preamble.
    # \usepackage[T1]{fontenc}
    "preamble": r"""
        %%%add number to subsubsection 2=subsection, 3=subsubsection
        %%% below subsubsection is not good idea.
        \setcounter{secnumdepth}{2}

        %%%% Table of content upto 2=subsection, 3=subsubsection
        \setcounter{tocdepth}{2}

        \usepackage[utf8]{inputenc}
        \usepackage[T1]{fontenc}
        \usepackage{helvet}
        \usepackage{amsfonts, amsmath, amssymb, mathbbol}
        \usepackage{graphicx}
        \usepackage{caption}
        \usepackage{xcolor}

        \definecolor{ocre}{RGB}{64,64,64}
        \usepackage[font={color=ocre}]{caption}

        %% unicode characters
        \usepackage{newunicodechar}
        \newunicodechar{⁹}{$^9$}
        \newunicodechar{⁸}{$^8$}
        \newunicodechar{⁷}{$^7$}
        \newunicodechar{⁶}{$^6$}
        \newunicodechar{⁵}{$^5$}
        \newunicodechar{⁴}{$^4$}
        \newunicodechar{³}{$^3$}
        \newunicodechar{²}{$^2$}
        \newunicodechar{¹}{$^1$}
        \newunicodechar{⁰}{$^0$}

        \newunicodechar{₉}{$_9$}
        \newunicodechar{₈}{$_8$}
        \newunicodechar{₇}{$_7$}
        \newunicodechar{₆}{$_6$}
        \newunicodechar{₅}{$_5$}
        \newunicodechar{₄}{$_4$}
        \newunicodechar{₃}{$_3$}
        \newunicodechar{₂}{$_2$}
        \newunicodechar{₁}{$_1$}
        \newunicodechar{₀}{$_0$}

        \newunicodechar{⟨}{$\langle$}
        \newunicodechar{⟩}{$\rangle$}
        \newunicodechar{−}{$-$}
        \newunicodechar{⟶}{$\longrightarrow$}
        \newunicodechar{Δ}{$\Delta$}
        \newunicodechar{⭐}{$\star$}


        %%% reduce spaces for Table of contents, figures and tables
        %%% it is used "\addtocontents{toc}{\vskip -1.2cm}" etc. in the document
        \usepackage[notlot,nottoc,notlof]{}

        \usepackage{setspace}
        \singlespacing

        %% table setting
        \renewcommand{\arraystretch}{1.75}

        %%%%%%%%%%% datetime
        \usepackage{datetime}

        \newdateformat{MonthYearFormat}{%
            \monthname[\THEMONTH], \THEYEAR}

        %% RO, LE will not work for 'oneside' layout.
        %% Change oneside to twoside in document class
        \usepackage{fancyhdr}
        \pagestyle{fancy}
        \fancyhf{}

        %%% Alternating Header for oneside
        \fancyhead[L]{
            \ifthenelse{\isodd{\value{page}}}{ \small \nouppercase{\leftmark} }{}
        }
        \fancyhead[R]{
            \ifthenelse{\isodd{\value{page}}}{}{ \small \nouppercase{\rightmark} }
        }

        %%% page number
        \fancyfoot[CO, CE]{\thepage}

        %%reduce spacing for itemize
        \usepackage{enumitem} \setlist{nosep}

        \makeatletter
            \renewcommand{\sphinxtableofcontents}{%
            %
            % before resetting page counter, let's do the right thing.
            \if@openright\cleardoublepage\else\clearpage\fi
            \addcontentsline{toc}{chapter}{Table of Contents}%
            \pagenumbering{roman}%
            \begingroup
                \parskip \z@skip
                \tableofcontents
            \endgroup
            %
            %% additional lists
            \if@openright\cleardoublepage\else\clearpage\fi
            \addcontentsline{toc}{chapter}{List of Figures}%
            \listoffigures
            %
            \if@openright\cleardoublepage\else\clearpage\fi
            \addcontentsline{toc}{chapter}{List of Tables}%
            \listoftables
            %
            \if@openright\cleardoublepage\else\clearpage\fi
            \addcontentsline{toc}{chapter}{List of Code Blocks}%
            \listof{literalblock}{List of Code Blocks}%
            %
            \if@openright\cleardoublepage\else\clearpage\fi
            \pagenumbering{arabic}%
            }
        \makeatother
    """,
    # "maketitle": r"""
    #     \pagenumbering{Roman} %%% to avoid page 1 conflict with actual page 1
    #     \begin{titlepage}
    #         \begin{figure}[!h]
    #             \centering
    #             \includegraphics[scale=0.3]{mrsimulator.png}
    #         \end{figure}
    #         %% \vfill adds at the bottom
    #         \vfill
    #     \end{titlepage}
    #     \clearpage
    #     \pagenumbering{roman}
    #     \tableofcontents
    #     \listoffigures
    #     \listoftables
    #     \clearpage
    #     \pagenumbering{arabic}
    # """,
    # Latex figure (float) alignment
    #
    # "figure_align": "htbp",
    "sphinxsetup": "\
        hmargin={0.7in,0.7in}, vmargin={1in,1in}, \
        verbatimwithframe=true, \
        HeaderFamily=\\rmfamily\\bfseries, \
        InnerLinkColor={rgb}{0,0.3,0.7}, \
        OuterLinkColor={rgb}{0.4,0,0.6}, \
        VerbatimBorderColor={rgb}{0.75,0.75,0.75}, \
        VerbatimColor={rgb}{0.98,0.98,0.985}, \
    ",
}  # TitleColor={rgb}{0.25,0.1,0.1}, \

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "mrsimulator.tex", "Mrsimulator Documentation", author, "manual")
]

# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "mrsimulator",
        "Mrsimulator Documentation",
        author,
        "mrsimulator",
        "Toolbox for simulating NMR spectrum.",
        "Miscellaneous",
    )
]
# ---------------------------------------------------------------------------- #


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (
        master_doc,
        "mrsimulator",
        "Mrsimulator Documentation",
        ["Deepansh J. Srivastava"],
        1,
    )
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_basename = project
epub_theme = "epub"
# epub_theme_options = {"relbar1": False, "footer": False}
epub_cover = ("_static/launch_2048x2732.png", "")
# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ["search.html", "_static/style.css"]


def setup(app):
    app.add_css_file("style.css")

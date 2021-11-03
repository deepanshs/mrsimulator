.. _installing_vscode:

Setting up VS Code
------------------

VS Code is a free code editor popular among the programming community. We recommend using VS Code
because of its flexibility and community extensions. The following is a guide to get the basics
set up and start running ``mrsimulator`` in VS Code. We also provide some suggested extension which
can make your programming much easier.

Instillation
""""""""""""

Microsoft offers detailed installation guides for all major platforms.
Visit the `VS Code setup page <https://code.visualstudio.com/docs/setup/setup-overview>`__ and
follow the instructions for your operating system.

We strongly recommend going through the `user interface guide
<https://code.visualstudio.com/docs/getstarted/userinterface>`__ for VS Code. You can also
find helpful tutorial and overview videos on this site as well as a more extensive list of
VS Code's features.

Suggested Extensions for running Python
"""""""""""""""""""""""""""""""""""""""

**Bare Minimums**

At bare minimum to run Python in VS Code, the `python extension
<https://marketplace.visualstudio.com/items?itemName=ms-python.python>`__ needs to be downloaded.
You can either follow the above hyperlink or search for "Python" in the VS Code app by clicking
the 5th icon from the top on the left-hand menu bar. The extension should have
*IntelliSense (Pylance)* in the description.

**Jupyter Notebooks**

If you want to run Jupyter Notebooks within VS Code, install the `Jupyter extension
<https://marketplace.visualstudio.com/items?itemName=ms-toolsai.jupyter>`__ from the marketplace.
To open a new notebook, press ``command + shift + p`` (``ctrl + shift + p`` on windows) or ``f1``
to bring up the command pallette. Type ``Create New Jupyter Notebook`` and select the first option.

.. note::

    VS Code may need to install other extensions for Notebooks to run. Install any additional
    extensions VS Code says are required.

VS Code may ask which Python environment to use. Make sure you select an environment where
``mrsimulator`` is installed. For more info on Python environments, see our `setting up
Anaconda <_virtual_environment_troubleshooting>`__ page.

Now you're ready to run Jupyter Notebooks in VS Code! Check out our
`examples <_example_gallery>`__ to download an example's notebook and open it with VS Code.

**Code Formatter**

The `Prettier - Code formatter <https://marketplace.visualstudio.com/items?itemName=esbenp.prettier-vscode>`__
lets you quickly format Python code to PEP standards. In any ``.py`` file, simply right-click
and select format.

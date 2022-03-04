===========================
Common Python Syntax Errors
===========================

Python syntax is slightly different than other langues which cause some confusion.
The following are some common syntax errors encountered and how to solve them.

However, the easiest way to find and prevent syntax errors is to use a dedicated code editor.
We recommend using VS Code on your local machine
or Google Colab which runs everything online. These programs check for syntax errors as you
write code.

IndentationError
""""""""""""""""

If you are encountering an ``IndentationError``, this means you have an extra/missing
whitespace somewhere in your code. Code editors make finding troublesome whitespace easier,
but the error should also show the code snippet which threw the error.

``IndentationError: expected an indented block`` means some code after a class/method/loop
deceleration is missing an indent.

``IndentationError: unindent does not match any outer indentation level`` means the code didn't
return to a previous indentation level.

``IndentationError: unexpected indent`` means Python encountered unexpected whitespace.

Code blocks in Python rely on indentation levels (1 level = 4 spaces) so whitespace can't
be placed randomly. Code blocks are preceded by a ``:`` and all code in one
block has the same indentation. To get out of a code block, simply remove an indentation level.

As an example of indentation, here is some code which adds the numbers 0 to 9:

.. plot::
    :context: close-figs

    # Add numbers 0 through 9
    total = 0
    for i in range(10):
        # New code block (4 spaces)
        total += i
    # Exit loop code block (0 spaces)
    print(total)

Mismatched Brackets and Square Brackets
"""""""""""""""""""""""""""""""""""""""

Nesting many lists and dictionaries inside each other quickly becomes hard to read. If you have
mismatched or missing brackets, Python will throw ``SyntaxError: invalid syntax``.
Code editors can automatically format large nestings and highlight which openings and closings go
together making the code easer to understand.

Make sure all brackets are balanced and that opening and closing brackets match. Python uses three
types of brackets:

* ``()`` is used when creating a `tuple <https://docs.python.org/3/library/stdtypes.html#tuple>`__ or when creating/calling method signatures.
* ``[]`` is used when creating a `list <https://docs.python.org/3/library/stdtypes.html#list>`__ or when indexing an item in a list or tuple.
* ``{}`` is used when creating a `dict <https://docs.python.org/3/library/stdtypes.html#mapping-types-dict>`__ or `set <https://docs.python.org/3/library/stdtypes.html#set>`__.

TypeError: object is not callable
"""""""""""""""""""""""""""""""""

The most common reason ``TypeError: object is not callable`` is when ``()`` is used instead of
``[]``. Parentheses are used to call functions. For example

.. plot::
    :context: close-figs

    def foo(n):
        print("I received", n)

    foo(1)
    # I received 1

But parentheses aren't valid for indexing a subscriptable object (list, tuple, etc). For
example, the following code will throw an TypeError

.. skip: start

.. plot::
    :context: close-figs

    bar = [1, 2, 3, 4]
    bar(1)

.. plot::

    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: 'list' object is not callable

.. skip:: end

but the following code is valid

.. plot::
    :context: close-figs

    bar = [1, 2, 3, 4]
    print(bar[1])
    # 2

The same applies for dictionaries, but instead of indexing with an integer you would index with
a keyword. For example

.. skip: start

.. plot::
    :context: close-figs

    spam = {"ham": "Hello World!", "eggs": 54.73}
    print(spam["ham"])  # prints Hello World!
    print(spam("ham"))

.. plot::

    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: 'dict' object is not callable

TypeError: object is not subscriptable
""""""""""""""""""""""""""""""""""""""

``TypeError: object is not subscriptable`` is thrown when indexing a non-subscriptable object.
For example

.. plot::
    :context: close-figs

    some_num = 42
    some_num[3]

.. plot::

    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: 'int' object is not subscriptable

.. skip: end

Also subscriptable objects can only be indexed so many times. A 1D list can only be indexed once,
2D twice, and so on. If you are using nested lists/dicts, make sure you aren't exceeding the number
of indexes possible.

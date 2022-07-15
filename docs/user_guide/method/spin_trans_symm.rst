.. note::
    
    In the `symmetry pathway approach
    <https://doi.org/10.1016/j.pnmrs.2010.11.003>`_,  the idea of coherence order is extended to form
    a complete set of spin transition symmetry functions, :math:`\xi_\ell
    (i,j)`, given by

    .. math::

        \xi_\ell(i,j) = \bra{j}  \hat{T}_{\ell,0} \ket{j} - \bra{i}  \hat{T}_{\ell,0} \ket{i},

    where the :math:`\hat{T}_{l,0}` are irreducible tensor operators.  The function
    symbol :math:`\xi_\ell(i,j)` is replaced with the lower-case symbols
    :math:`\mathbb{p}(i,j)`, :math:`\mathbb{d}(i,j)`, :math:`\mathbb{f}
    (i,j)`, :math:`\ldots`, i.e., we follow the spectroscopic sub-shell letter
    designations:

    .. math::

        \begin{array}{cccccccccccccccl}
        \ell = & 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10  &11  &12  &13  & \leftarrow \text{numerical value} \\
        \xi_\ell \equiv	& \mathbb{s} &  \mathbb{p} &  \mathbb{d} &  \mathbb{f} &  \mathbb{g} &  \mathbb{h} &  \mathbb{i} & \mathbb{k} &\mathbb{l} & \mathbb{m} & \mathbb{o} & \mathbb{q} & \mathbb{r} &\mathbb{t} & \leftarrow \text{symbol}\\
        \end{array}

    To simplify usage in figures and discussions, we scale the transition symmetry
    functions to integers values according to

    .. math::

        \text{p}(i,j) = \mathbb{p}(i,j), ~~~~~
        \text{d}(i,j) = \sqrt{\frac{2}{3}} \, \mathbb{d}(i,j), ~~~~~
        \text{f}(i,j) = \sqrt{\frac{10}{9}} \, \mathbb{f}(i,j),
        ~~~~~
        \cdots

    The :math:`\ell=0` function is dropped as it always evaluates to zero. For a
    single spin, :math:`I`, a complete set of functions is defined up to
    :math:`\ell = 2I`.

    For weakly coupled nuclei, we define the transition symmetry functions

    .. math::

        \xi_{\ell_1,\ell_2, \ldots, \ell_n} (i,j) =
        \left \langle j \right|\hat{T}_{\ell_1,0}({\bf I}_1)\hat{T}_{\ell_2,0}({\bf I}_2)\ldots\hat{T}_{\ell_n,0}({\bf I}_n) \left|j \right \rangle
        -
        \left \langle i \right|\hat{T}_{\ell_1,0}({\bf I}_1)\hat{T}_{\ell_2,0}({\bf I}_2)\ldots\hat{T}_{\ell_n,0}({\bf I}_n) \left|i \right \rangle

    Replacing the symmetry function symbol using sub-shell letter designations becomes
    more cumbersome in this case.  When the :math:`\ell` are zero on all nuclei except one,
    we identify these functions as

    .. math::

        \begin{array}{cccc}
        \mathbb{p}_1 = \xi_{1,0, \ldots, 0} (i,j), &
        \mathbb{p}_2 = \xi_{0,1, \ldots, 0} (i,j), &
        \ldots, &
        \mathbb{p}_n = \xi_{0,0, \ldots, 1} (i,j),\\
        \\
        \mathbb{d}_1 = \xi_{2, 0, \ldots, 0} (i,j), &
        \mathbb{d}_2 = \xi_{0,2, \ldots, 0} (i,j), &
        \ldots, &
        \mathbb{d}_n = \xi_{0,0, \ldots, 2} (i,j), \\
        \\
        \mathbb{f}_1 = \xi_{3, 0, \ldots, 0} (i,j), &
        \mathbb{f}_2 = \xi_{0,3, \ldots, 0} (i,j), &
        \ldots, &
        \mathbb{f}_n = \xi_{0,0, \ldots, 3} (i,j), \\
        \vdots & \vdots &  & \vdots
        \end{array}

    For weakly coupled homonuclear spins it is also convenient to define

    .. math::

        \begin{array}{c}
        \mathbb{p}_{1,2,\ldots,n} =  \mathbb{p}_{1}
        + \mathbb{p}_{2} + \cdots \mathbb{p}_{n}, \\
        \\
        \mathbb{d}_{1,2,\ldots,n} =  \mathbb{d}_{1}
        + \mathbb{d}_{2} + \cdots \mathbb{d}_{n}, \\
        \\
        \mathbb{f}_{1,2,\ldots,n} =  \mathbb{f}_{1}
        + \mathbb{f}_{2} + \cdots \mathbb{f}_{n}, \\
        \vdots
        \end{array}

    When the :math:`\ell` are zero on all nuclei except two, then we identify
    these functions using a concatenation of sub-shell letter designations, e.g.,

    .. math::

        \begin{array}{cccc}
        (\mathbb{pp})_{1,2} = \xi_{1,1,0, \ldots, 0} (i,j), &
        (\mathbb{pp})_{1,3} = \xi_{1,0,1, \ldots, 0} (i,j), &
        \ldots, &
        (\mathbb{pp})_{1,n} = \xi_{1,0,0, \ldots, 1} (i,j),\\
        \\
        (\mathbb{pd})_{1,2} = \xi_{1, 2, 0, \ldots, 0} (i,j), &
        (\mathbb{pd})_{1,3} = \xi_{1,0,2 \ldots, 0} (i,j), &
        \ldots, &
        (\mathbb{pd})_{1,n} = \xi_{1,0, \ldots, 2} (i,j), \\
        \\
        (\mathbb{dp})_{1,2} = \xi_{2, 1, 0, \ldots, 0} (i,j), &
        (\mathbb{dp})_{1,3} = \xi_{2 ,0, 1 \ldots, 0} (i,j), &
        \ldots, &
        (\mathbb{dp})_{1,n} = \xi_{2, 0, \ldots, 1} (i,j), \\
        \vdots & \vdots &  & \vdots
        \end{array}

    Below is an energy level diagram of two coupled spin :math:`I=1/2` nuclei with
    transition labeled according to their transition symmetry function values.  Note
    that each transition has a unique set of transition symmetry function values.

    .. figure:: ../../_static/CoupledOneHalf.*
        :width: 650
        :alt: figure
        :align: center

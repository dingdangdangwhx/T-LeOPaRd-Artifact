"""Parameter selection for LeOPaRd and its threshold version T-LeOPaRd.

The modules follow the order in which the paper chooses the parameters:

* ``widths``      Gaussian widths from Theorem ``thm:MLWEP-to-HMLWEP``,
* ``modulus``     the modulus q from the uniqueness bound p_u,
* ``hardness``    MSIS and MLWE hardness, i.e. the ranks n_s, n_c, m and ell,
* ``parameters``  the fixed-point search that ties these together,
* ``sizes``       communication and public-key sizes,
* ``threshold``   the same for the threshold scheme T-LeOPaRd.

Labels such as ``thm:MLWEP-to-HMLWEP`` refer to the LaTeX sources of the
LeOPaRd paper, except in ``threshold``, which refers to the threshold paper.
"""

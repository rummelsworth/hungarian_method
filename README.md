<!-- This is a Markdown-formatted document. -->

The Hungarian method: an implementation in C
============================================

###### by William Rummler (w.a.rummler@gmail.com), 21 November 2010

Package Contents
----------------

+   `hungarian_method.c`

+   `hungarian_method.h`

+   `brute_force_assignment.c`

+   `brute_force_assignment.h`

+   `hm_test.c`

+   `supplement.pdf`

+   `README`

+   `COPYING`

Summary
-------

This implementation of the Hungarian method is derived almost entirely from
Chapter 11 of _Combinatorial Optimization: Algorithms and Complexity_ by
Christos Papadimitriou and Kenneth Steiglitz. Specifically, the code within
`hungarian_method.c` should be able to be read with an approximately one-to-one
correspondence to the pseudo-code presented in Figure 11-2. This package also
contains an implementation of a brute-force solution to the _assignment
problem_, the problem that the Hungarian method solves so much more
efficiently. The brute-force implementation is included for the sake of
comparison and testing. Finally, there is a very basic testing program
contained in `hm_test.c`.

Just to be totally clear, here is an example command line that compiles and
runs the test program:

    $ gcc hm_test.c brute_force_assignment.c hungarian_method.c
    $ ./a.out

Motivation
----------

This source package was prepared while I was debugging an initial attempt at an
implementation of the Hungarian method according to P&S Figure 11-2. After much
testing of the code and examination of the text, especially Example 11.1, I
found what appeared to be several errors in the pseudo-code. The purpose of
this source package is to document these errata with an appropriately
annotated, working implementation.

To offer somewhat more concise documentation of the errata, this package also
includes a supplementary errata file (`supplement.pdf`) whose format is
essentially the same as that of the [8 October 2000 errata file][latest]
located at Prof. Steiglitz's Princeton homepage as of 21 November 2010. For
details-in-context on the errata, please view the comments embedded within the
implementation file `hungarian_method.c`.

[latest]: http://www.cs.princeton.edu/~ken/latest.pdf

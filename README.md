# The Hungarian method: an implementation in C/C++

- 5 March 2017
- 21 November 2010

## Summary

This implementation of the Hungarian method is derived almost entirely from
Chapter 11 of _Combinatorial Optimization: Algorithms and Complexity_ by
Christos Papadimitriou and Kenneth Steiglitz. Specifically, the code within
`hungarian_method.cc` should be able to be read with an approximately one-to-one
correspondence to the pseudo-code presented in Figure 11-2. This package also
contains an implementation of a brute-force solution to the _assignment
problem_, the problem that the Hungarian method solves so much more efficiently.
The brute-force implementation is included for the sake of comparison and
testing. Finally, there is a very basic testing program contained in
`hm_test.cc`.

This implementation was _not_ designed as a reusable library, with qualities
like API user-friendliness and performance in mind. The purpose of development
was didactic, to produce as close a correct analog as possible of the Figure
11-2 pseudocode.

The source was originally written in ANSI C to be compiled with a plain `gcc`
command line, but has since been updated (with very few changes) to compile as a
Visual C++ project in Visual Studio 2017.

## Motivation

This source package was prepared while I was debugging an initial attempt at an
implementation of the Hungarian method according to P&S Figure 11-2. After much
testing of the code and examination of the text, especially Example 11.1, I
found what appeared to be several errors in the pseudo-code. The purpose of this
source package is to document these errata with an appropriately annotated,
working implementation.

To offer somewhat more concise documentation of the errata, this package also
includes a [supplementary errata file](supplement.pdf) whose format is
essentially the same as that of the [8 October 2000 errata file][latest] located
at Prof. Steiglitz's Princeton homepage as of 21 November 2010. For
details-in-context on the errata, please view the comments embedded within the
implementation file `hungarian_method.cc`.

[latest]: http://www.cs.princeton.edu/~ken/latest.pdf

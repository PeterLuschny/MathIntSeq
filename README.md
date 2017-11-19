# MathIntSeq

Julia implementations of mathematical integer sequences.

For Julia 0.6.1 on Linux, Mac OS X and Win64.

[![Build status](https://travis-ci.org/PeterLuschny/MathIntSeq.svg?branch=master)](https://travis-ci.org/PeterLuschny/MathIntSeq) [![Build status](https://ci.appveyor.com/api/projects/status/d60u86n44fc43awi?svg=true)](https://ci.appveyor.com/project/PeterLuschny/mathintseq)


First steps and usage can be seen in this

- [Example notebook](http://olms.onl/julia/MathIntSeq/notebook/SeqNotebook.html).

For a quick first view install MathIntSeq and run the

- [demo](https://github.com/OpenLibMathSeq/MathIntSeq/blob/master/src/demo.jl).

The API documentation is located at the

- [Online documentation](http://olms.onl/julia/MathIntSeq).

MathIntSeq serves also as an Julia interface to the On-Line Encyclopedia of
Integer Sequences.

- [OEIS](http://oeis.org/)

In the OEIS you will find comments and references which explain the sequences.

MathIntSeq is based on the multiprecision arithmetic of the computer algebra package

- [Nemo](http://nemocas.github.io/Nemo.jl/latest/)

We also use some of Nemo's libraries. [Nemo](http://nemocas.org/) was designed and written by William Hart and is maintained by him and others.

If you want to contribute to the 'Prog' section in the OEIS, this would be a
typical entry (here on the page A000041):

    (Julia)
    using MathIntSeq
    println(L000041(50)) # _your-OEIS-signature_

The example returns the list of the first 50 Partition numbers (number of
partitions of n).    

If you want to contribute to MathIntSeq see

- [Contributing](https://github.com/OpenLibMathSeq/MathIntSeq/blob/master/CONTRIBUTING.md).                

Any feedback is welcome. Issues and suggestions should be posted to

-  [MathIntSeq's issue tracker](https://github.com/OpenLibMathSeq/MathIntSeq/issues).

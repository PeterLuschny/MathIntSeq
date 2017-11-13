# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module OLMS
using MathIntSeq

doc"""
Generating the OLMS logo.

⍊(n, k) means k is strong prime to n and ι() are the Iverson brackets.
Insert the line sum(L) == n - 5 && println(n) if you want to see why this
triangle is called the Save Prime Triangle.
"""
function SavePrimeTriangle()
    for n in 5:23
        T = [k for k in 1:n if ⍊(n, k)]
        L = [ι(k ∈ T) for k in 3:(n - 2)]
        println(L)
    end
end

SavePrimeTriangle()

end # module

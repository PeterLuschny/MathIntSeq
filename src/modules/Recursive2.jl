# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Recursive2
using SeqBase, Nemo

export Recurrence2, L053602, L002467, L005442, L058797, L001040, L058279
export L001046, L036246, L024167, L056953, L286032, L000246, L166474

# In this module we look at recursive sequences of the form
# a(n) = f(n) × a(n-1) + g(n) × a(n-2) with given a(0) and a(1) (default 1),
# both f(n) and g(n) not zero and at least one of the functions f(n), g(n)
# _not_ constant.
#
# This recurence pattern includes the factorial and the Fibonacci numbers.
# The references to the OEIS are approximate. Offset, signs and the first few
# values might differ (often a '0' is prepended and sometimes the leading '1'
# is deleted in the OEIS).

doc"""
Computes recursive sequences of the form
``a(n) = f(n) a(n-1) + g(n) a(n-2)`` with given initial values ``a(0)`` and
``a(1)`` (with defaults 1), ``f(n)`` and ``g(n)`` not zero and at least one
of the functions ``f(n)``, ``g(n)`` is _not_ constant.
"""
function Recurrence2(len::Int, f::Function, g::Function, a1=1, a2=1)
    len ≤ 0 && return fmpz[]
    A = Array{fmpz}(len)
    for n in 0:len - 1
        if n == 0
            A[1] = a1
        elseif n == 1
            A[2] = a2
        else
            A[n + 1] = f(n) * A[n] + g(n) * A[n - 1]
        end
    end
    SeqArray(A)
end

doc"""
Return an array of length len of the Hankel transform of the Bessel numbers
starting at n = 1. (David Callan)
"""
L058797(len::Int) = Recurrence2(len, n -> n, n -> -1)

doc"""
Return an array of length len of the continuant transform of ``1, 2, 3, 4, 5, …``.
"""
L001040(len::Int) = Recurrence2(len, n -> n, n -> 1)

doc"""
Return the number of palindromic compositions of ``n`` into odd parts. (Emeric Deutsch)
"""
L053602(len::Int) = Recurrence2(len, n -> 1, n -> (-1)^n)

doc"""
Return the number of permutations in the symmetric group ``S_n`` that have
a fixed point.
"""
L002467(len::Int) = Recurrence2(len, n -> n, n -> n)

doc"""
Return an array of length len of the continuant transform of squares ``1, 4, 9, …``.
"""
L036246(len::Int) = Recurrence2(len, n -> n^2, n -> 1)

doc"""
Return an array of length len of the sequence with generating function
``2 ((√(2) + x)^2/(2 - x^2))^{\frac{1}{√(2)}} /(2 - x^2)``. (Vaclav Kotesovec)
"""
L166474(len::Int) = Recurrence2(len, n -> 1, n -> div(n^2 - n, 2))

doc"""
Return an array of length len of the sequence with exponential generating function
``2 π(I`` Bessel``Y(3, 2I)`` Bessel``I(2, 2 √(1-x)) +``
Bessel``I(3, 2)`` Bessel``Y(2, 2 I √(1-x)))/(1-x)``. (Wolfdieter Lang)
"""
L058279(len::Int) = Recurrence2(len, n -> n + 1, n -> 1)

doc"""
``a(n) = n(n-1)a(n-1)/2 + a(n-2)``.
"""
L001046(len::Int) = Recurrence2(len, n -> div(n^2 - n, 2), n -> 1)

doc"""
Return the first len terms defined by
``a(n) = n! [x^n] (1 - √(Pi/2) \exp(-((x - 1)^2) / 2) (x - 1)
(erfi((x - 1) / √(2)) + erfi(1 / √(2))))``. For 'erfi' see MathWorld.
"""
L286032(len::Int) = Recurrence2(len, n -> 1, n -> -n)

doc"""
Return the number of n-permutations that have a cycle with length greater
than ``n/2``. (Geoffrey Critzer)
"""
L024167(len::Int) = Recurrence2(len, n -> 1, n -> n^2)

doc"""
Return the number of involutions such that every 2-cycle contains one odd and
one even element. (Alois P. Heinz)
"""
L056953(len::Int) = Recurrence2(len, n -> 1, n -> div(n, 2))

doc"""
Return the number of permutations in the symmetric group ``S_n`` that have odd order.
"""
L000246(len::Int) = Recurrence2(len, n -> 1, n -> (n - 1) * (n - 2))

doc"""
Return ``n!`` Fibonacci``(n+1)``.
 """
L005442(len::Int) = Recurrence2(len, n -> n, n -> n^2 - n)

# Fibonacci numbers
# L000045(len::Int) = Recurrence2(len, n -> -(-1)^n, n -> -1)

# Factorial numbers
# L000142(len::Int) = Recurrence2(len, n -> n-1, n -> n-1)

# Similar sequences:
# L058798(len::Int) = Recurrence2(len, n -> n, n -> -1, 0)
# L001053(len::Int) = Recurrence2(len, n -> n, n -> 1, 0)
# L051792(len::Int) = Recurrence2(len, n -> 1, n -> (-1)^n, 0)
# L002467(len::Int) = Recurrence2(len, n -> n, n -> n, 0)

end # module

module Recursive2Test
using Base.Test, SeqTests, SeqBase, Recursive2, Nemo, OEISUtils

function test()

# In the OEIS sometimes a '0' is prepended or a '1' is missing.
# N.B. Refrences to the OEIS are always approximative!
Data = Dict{Int, Array{fmpz}}(
058797 => [1, 1, 1, 2, 7, 33, 191, 1304, 10241],
001040 => [1, 1, 3, 10, 43, 225, 1393, 9976, 81201],
058279 => [1, 1, 4, 17, 89, 551, 3946, 32119, 293017],
001046 => [1, 1, 2, 7, 44, 447, 6749, 142176, 3987677],
036246 => [1, 1, 5, 46, 741, 18571, 669297, 32814124, 2100773233],
024167 => [1, 1, 5, 14, 94, 444, 3828, 25584, 270576],
056953 => [1, 1, 2, 3, 7, 13, 34, 73, 209],
286032 => [1, 1, -1, -4, 0, 20, 20, -120, -280],
000246 => [1, 1, 1, 3, 9, 45, 225, 1575, 11025],
166474 => [1, 1, 2, 5, 17, 67, 322, 1729, 10745],
053602 => [1, 1, 2, 1, 3, 2, 5, 3, 8],
002467 => [1, 1, 4, 15, 76, 455, 3186, 25487, 229384],
005442 => [1, 1, 4, 18, 120, 960, 9360, 105840, 1370880]
)

    @testset "Recursive2" begin

        if oeis_isinstalled()

            Seq = [ L058797, L001040, L058279, L001046, L036246, L024167, L056953,
                L286032, L000246, L166474, L053602, L002467, L005442 ]

            for seq in Seq
                S = seq(9)
                anum = SeqNum(seq)
                data = SeqArray(Data[anum])
                # println(anum); println(S); println(data)
                @test all(S[0:8] .== data[0:8])
            end
        end
    end
end

function demo()
    println("A058797 ", L058797(12))
end

doc"""
L286032(10000) :: 0.095866 seconds (30.00 k allocations: 546.891 KB)
"""
function perf()
    gc()
    @time L286032(10000)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

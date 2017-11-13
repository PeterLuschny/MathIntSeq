# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Deleham
using Nemo, SeqBase, NumberTheory

export DeléhamΔ, T084938, T060693, T106566, T094665, T090238, T225478
export T055883, T184962, T088969, T090981, T011117

# Applying Deléham's T-operation often gives an additional first column
# or an additional main diagonal in the resulting triangle compared to what is
# listed in the OEIS.

doc"""
Return the product of two integer sequences introduced by Philippe Deléham in
[A084938](http://oeis.org/A084938).
"""
function DeléhamΔ(n::Int, S::Function, T::Function)
    n ≤ 0 && return fmpz[]
    R, x = PolynomialRing(ZZ, "x")
    A = [R(S(k) + x * T(k)) for k in 0:n - 2]
    C = [R(1) for i in 0:n]; C[1] = R(0)
    len = div((n + 1) * n, 2)
    M = SeqArray(len)
    m = 0

    for k in 0:n - 1
        for j in k + 1:-1:2
            C[j] = C[j - 1] + C[j + 1] * A[j - 1]
        end
        for j in 0:k
            M[m] = coeff(C[2], j)
            m += 1
        end
    end
    M
end

doc"""
Return the number of permutations of ``{1,2,...,n}`` having ``k`` cycles such that the
elements of each cycle of the permutation form an interval. (Ran Pan)
"""
T084938(n::Int) = DeléhamΔ(n, i -> div(i + 1, 2), i -> 0^i)

doc"""
Return the number of lattice paths from ``(0,0)`` to ``(x,y)`` that never pass below ``y = x``
and use step set ``{(0,1), (1,0), (2,0), (3,0), ...}``.
"""
T011117(n::Int) = DeléhamΔ(n, i -> 0^i, i -> IsOdd(i) ? 1 : (i > 0 ? 2 : 0))

doc"""
Return the number of Schroeder paths (i.e., consisting of steps ``U=(1,1), D=(1,-1), H=(2,0)``
and never going below the x-axis) from ``(0,0)`` to ``(2n,0)``, having ``k`` peaks.
(Emeric Deutsch)
"""
T060693(n::Int) = DeléhamΔ(n, i -> 1, i -> IsOdd(i) ? 0 : 1)

doc"""
Return the the Catalan convolution triangle.
"""
T106566(n::Int) = DeléhamΔ(n, i -> i == 0 ? 0 : 1, i -> i == 0 ? 1 : 0)

doc"""
Return the number of increasing 0-2 trees (A002105) on 2n edges in which the minimal path
from the root has length k. (David Callan)
"""
T094665(n::Int) = DeléhamΔ(n, i -> div(i * (i + 1), 2), i -> i + 1)

doc"""
Return the number of lists of k unlabeled permutations whose total length is n. (David Callan)
"""
T090238(n::Int) = DeléhamΔ(n, i -> div(i, 2) + (IsOdd(i) ? 2 : 0), i -> i == 0 ? 1 : 0)

doc"""
Return the triangle ``4^k S_4(n, k)`` where ``S_m(n, k)`` are the Stirling-Frobenius cycle
numbers of order m.
"""
T225478(n::Int) = DeléhamΔ(n, i -> 2(i + 1) + (i + 1) % 2, i -> IsOdd(i) ? 0 : 4)

doc"""
Return the exponential transform of Pascal's triangle.
"""
T055883(n::Int) = DeléhamΔ(n, i -> IsOdd(i) ? div(i + 1, 2) : 1, i -> IsOdd(i) ? div(i + 1, 2) : 1)

doc"""
Return the number of Schroeder paths of length ``2n`` and having ``k`` ascents.
"""
T090981(n::Int) = DeléhamΔ(n, i -> i == 0 ? 1 : (IsOdd(i) ? 0 : 2), i -> IsOdd(i) ? 1 : 0)

doc"""
Return a triangle related to the median Euler numbers.
"""
T088969(n::Int) = DeléhamΔ(n, i -> i^2, i -> IsOdd(i) ? 3div(i, 2) + 2 : 5div(i, 2) + 1)

doc"""
Return the Bell transform of the Fubini numbers.
"""
T184962(n::Int) = DeléhamΔ(n, i -> div((i + 1) - (i + 1) % 2, 2 - (i + 1) % 2), i -> IsOdd(i) ? 0 : 1)

end # module

module DelehamTest
using Base.Test, SeqBase, Nemo, Deleham, OEISUtils

# References to the OEIS A-numbers are always approximately only!

function test()

Data = Dict{Int, Array{fmpz}}(
084938 => [1, 0, 1, 0, 1, 1, 0, 2, 2, 1],
060693 => [1, 1, 1, 2, 3, 1, 5, 10, 6, 1],
106566 => [1, 0, 1, 0, 1, 1, 0, 2, 2, 1],
094665 => [1, 0, 1, 0, 1, 3, 0, 4, 15, 15],
090238 => [1, 0, 1, 0, 2, 1, 0, 6, 4, 1],
225478 => [1, 3, 4, 21, 40, 16, 231, 524, 336, 64],
055883 => [1, 1, 1, 2, 4, 2, 5, 15, 15, 5],
184962 => [1, 0, 1, 0, 1, 1, 0, 3, 3, 1],
088969 => [1, 0, 1, 0, 1, 3, 0, 5, 20, 21],
090981 => [1, 1, 0, 1, 1, 0, 1, 4, 1, 0],
011117 => [1, 1, 0, 1, 1, 0, 1, 2, 3, 0]
)

#[0] 1
#[1] 1		1
#[2] 2		4		2
#[3] 5		15		15		5
#[4] 15		60		90		60		15
#[5] 52		260		520		520		260		52
#[6] 203	1218	3045	4060	3045	1218	203

    @testset "Deléham" begin

        n = 7
        B = SeqArray([bell(n) * binomial(n, j) for j in 0:n])
        R = Row(T055883(n+1), n)
        @test all(B .== R)

        a = SeqArray([0, 5040, 2208, 828, 272, 70, 12, 1])
        b = Row(T090238(8), 7)
        @test all(a .== b)

        Seq = [ T084938, T060693, T106566, T094665, T090238, T225478, T055883,
              T184962, T088969 ]

        for seq in Seq
            S = seq(10)
            anum = SeqNum(seq)
            data = SeqArray(Data[anum])
            # println(anum); println(S); println(data)
            @test all(S[0:9] .== data[0:9])
        end
    end
end

function demo()
    dim = 5
    for n in 0:dim
        println(T055883(n))
    end

    Show(T055883(dim))
end

doc"""
T225478(100) :: 0.071950 seconds (15.27 k allocations: 638.438 KiB)
T055883(100) :: 0.076918 seconds (15.27 k allocations: 638.438 KiB)
"""
function perf()
    gc()
    @time T225478(100)
    @time T055883(100)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

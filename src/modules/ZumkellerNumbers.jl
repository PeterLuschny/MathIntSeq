# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module ZumkellerNumbers
using Nemo, SeqBase, NumberTheory, CombinationsIterator

export IsZumkeller, ZumkellerNumberList, L083207

doc"""
Is ``n`` a Zumkeller number?

A Zumkeller number ``n`` is an integer whose divisors can be partitioned
into two disjoint sets whose sums are both ``σ(n)/2``.
"""
function IsZumkeller(n::Int)
    n == 0 && return false
    T = Divisors(n)
    s = sum(T)
    ((s % 2 ≠ 0) || (s < 2n)) && return false
    S = s >> 1 - n
    D = [d for d in T if d ≤ S]
    D == [] && return true
    for c in combinations(D)
        S == sum(c) && return true
    end
    return false 
end

doc"""
Return a list of length len of the first Zumkeller numbers.
"""
ZumkellerNumberList(len) = SeqArray(len, IsZumkeller)

doc"""
Return a list of length len of the first Zumkeller numbers.
"""
L083207(len) = SeqArray(len, IsZumkeller)

end # module

module ZumkellerNumbersTest
using Base.Test, Nemo, NumberTheory, ZumkellerNumbers

function test()
    @testset "Zumkeller" begin
        @test IsZumkeller(17000) == true
        @test IsZumkeller(27472) == true
        @test IsZumkeller(29062) == false
        @test IsZumkeller(43464) == true
    end
end

function demo()
    println(ZumkellerNumberList(10))

    for n in 20:30
        println(n, " ↦ ",  IsZumkeller(n))
    end
    for n in 0:6
        println(n, " ↦ ", ZumkellerNumberList(n))
    end
end

doc"""
for n in 1:2000 IsZumkeller(n) end :: 0.476108 seconds (2.89 M allocations: 119.340 MB, 25.32% gc time)
"""
function perf()
    IsZumkeller(10)
    gc()
    @time (for n in 1:2000 IsZumkeller(n) end)
end

function main()
    demo()
    test()
    perf()
end

main()

end # module

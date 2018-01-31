# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Derangements
using Nemo
export Derangement, Subfactorial, A000166, A000255, A105927

doc"""
Return the number of fixpoit free permutations of an n-set.
"""
function Derangement(n)
    prec = 250
    while prec < 10000
        CC = AcbField(prec)
        c = exp(CC(-1)) * Nemo.gamma(CC(n+1), CC(-1))
        b, i = unique_integer(c)
        b && return i
        prec *= 2
    end
    println("n = $n gives an InexactError!")
    # throw(InexactError())
end

doc"""
Return the number of fixpoit free permutations of an n-set.
"""
Subfactorial(n) = Derangement(n)
doc"""
Return the number of fixpoit free permutations of an n-set.
"""
A000166(n) = Derangement(n)

doc"""
Return the number of fixpoit free permutations of an n-set beginning with 2.
"""
function A000255(n)
    prec = 250
    while prec < 10000
        CC = AcbField(prec)
        # Do not use Subfactorial(n+2)/(n+1) !
        c = exp(CC(-1)) * Nemo.gamma(CC(n+3), CC(-1)) / CC(n+1)
        b, i = unique_integer(c)
        b && return i
        prec *= 2
    end
    println("n = $n gives an InexactError!")
    # throw(InexactError())
end

doc"""
Return the number of derangements of a set of n + 2 elements such that the first
two elements belong to distinct cycles. (I. Mezo)
"""
function A105927(n)
    prec = 250
    while prec < 10000
        CC = AcbField(prec)
        c = (CC(n^2+n-1)*exp(CC(-1)) * Nemo.gamma(CC(n+1),CC(-1)) - (-1)^n*CC(n-1))/CC(2)
        b, i = unique_integer(c)
        b && return i
        prec *= 2
    end
    println("n = $n gives an InexactError!")
    # throw(InexactError())
end

# see also A105928

end # module

module DerangementsTest
using Base.Test, SeqBase, SeqTests, Nemo, OEISUtils, Derangements

function test()
    @testset "Derangements" begin
        @test Derangement(1) == 0
        @test A000255(8) == 148329
        @test A105927(6) == 5430

        if oeis_isinstalled()
            A = [A000166, A000255, A105927]
            SeqTest(A, 'A')
        end
    end
end

function demo()
    for n in 0:30 println(Derangement(n)) end
    println(Derangement(450))
    for n in 0:30 println(A000255(n)) end
    println(A000255(100))
    for n in 0:30 println(A105927(n)) end
    println(A105927(300))
end

doc"""
(for n in 0:450 Subfactorial(n) end) ::  0.342082 seconds (55.69 k allocations: 1.775 MiB)
"""
function perf()
    gc()
    @time (for n in 0:450 Subfactorial(n) end)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

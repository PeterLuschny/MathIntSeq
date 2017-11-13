# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Abundant
using Nemo, SeqBase, NumberTheory

export IsAbundant, C005101, L005101, I005101, B005101, C002093, L002093

doc"""
Is ``n`` an abundant number, i.e. is ``σ(n) > 2n``?
"""
IsAbundant(n) = σ(n) - n > n

doc"""
Generate the abundant numbers which are the numbers such that ``σ(n) > 2n``.
"""
C005101() = Channel(csize=4) do c
    n = 0
    while true
        IsAbundant(n) && put!(c, ZZ(n))
        n += 1
    end
end

doc"""
Generate the highly abundant numbers which are the numbers where record values
of ``σ(n)`` occur.
"""
C002093() = Channel(csize=2) do c
    for n in 1:4 put!(c, ZZ(n)) end
    n = high = ZZ(4)
    while true
        n += 2
        s = σ(n)
        if s > high
            high = s
            put!(c, n)
        end
    end
end

doc"""
Return a list of length len of highly abundant numbers.

julia> L002093(8)
[0, 1, 2, 3, 4, 6, 8, 10, 12]
"""
L002093(len) = SeqArray(len, C002093())

doc"""
Return a list of length len of abundant numbers.

julia> L005101(9)
[12, 18, 20, 24, 30, 36, 40, 42, 48]
"""
L005101(len) = SeqArray(len, IsAbundant)

# We include an alternative implementation to compare the running times.
doc"""
Return a list of length len of abundant numbers. (L005101 is faster.)
"""
I005101(len) = SeqArray(len, C005101())

doc"""
Return a list of abundant numbers which do not exceed the given bound.

julia> B005101(50)
[12, 18, 20, 24, 30, 36, 40, 42, 48]
"""
B005101(bound) = FindUpTo(bound, IsAbundant)

end # module

module AbundantTest
using Abundant, Base.Test, SeqTests, SeqBase, OEISUtils, Nemo

function test()
    @testset "Abundant" begin
        @test IsAbundant(100800) == true
        @test IsAbundant(2402400) == true
        @test IsAbundant(49008960) == true

        channel = C002093()
        for _ in 1:9 take!(channel) end
        @test take!(channel) == 18
        close(channel)

        @test SeqSize(B005101(100)) == 22

        if oeis_isinstalled()
            SeqTest([L005101, L002093], 'L')
            SeqTest([B005101], 'B')
        end
    end
end

function demo()

    println(L005101(5))
    println(B005101(30))
    println()

    c = C002093()
    foreach(i -> println(i), Iterator(20, c))
    close(c)

    c = C002093()
    for i in Iterator(20, c) println(i) end
    close(c)

    c = C002093()
    println(collect(Iterator(20, c)))
    close(c)

    c = C002093()
    println(SeqArray(20, c))
    close(c)

    for n in 40:50 println(n, " ↦ ", IsAbundant(n)) end
    for n in 0:6   println(n, " ↦ ", L005101(n)) end

    println(L002093(11))
    for n in 0:6   println(n, " ↦ ", L002093(n)) end
end

# Unsurprisingly the channel implementation is slower.
doc"""
L005101(80000) :: 0.558171 seconds (1.05 M allocations: 16.634 MiB)
I005101(80000) :: 0.798868 seconds (1.17 M allocations: 18.468 MiB)
L002093(200)   :: 0.773658 seconds (720.73 k allocations: 11.006 MB, 6.30% gc time)
"""
function perf()
    L005101(10)
    @time L005101(80000)
    gc()
    I005101(10)
    @time I005101(80000)
    gc()
    L002093(10)
    @time L002093(200)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

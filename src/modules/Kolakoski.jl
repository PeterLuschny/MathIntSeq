# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Kolakoski

export KolakoskiList, C000002, L000002

doc"""
Generates the Kolakoski sequence which is the unique sequence
over the alphabet ``{1, 2}`` starting with ``1`` and having the
sequence of run lengths identical with itself.
"""
C000002() = Channel(csize=10) do c
    x = y = Int(-1)

    while true
        put!(c, [2, 1][(x & 1) + 1])
        f = y & ~(y + 1)
        x = xor(x, f)
        y = (y + 1) | (f & (x >> 1))
    end
end

doc"""
Return the list of the first ``n`` terms of the Kolakoski sequence.
"""
function KolakoskiList(len::Int)
    len ≤ 0 && return []
    generator = C000002()
    L = [take!(generator) for _ in 1:len]
    close(generator)
    L
end

doc"""
Return the list of the first ``n`` terms of the Kolakoski sequence.
"""
L000002(n::Int) = KolakoskiList(n)

end # module

module KolakoskiTest
using Base.Test, Kolakoski

function test()

    @testset "Kolakoski" begin
        K = KolakoskiList(100)
        @test K[1]  == 1
        @test K[33] == 2
        @test K[72] == 2

        generator = C000002()
        for n in [1, 33, 72]
            k = take!(generator)
            @test K[n] == k
        end
        close(generator)
    end
end

function demo()
    println(KolakoskiList(20))

    generator = C000002()
    o = e = 0
    for n in 1:80
        take!(generator) == 1 ? o += 1 : e += 1
        print(o - e, " ")
    end
    println()
    close(generator)
end

doc"""
KolakoskiList(10000) :: 0.013653 seconds (12.04 k allocations: 1.028 MiB)
"""
function perf()
    @time KolakoskiList(10000)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

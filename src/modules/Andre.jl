# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Andre
using Memoize, Nemo

export André, C000111, A000111, A178963, A178964, A181936, A250283

# n\k [0][1][2][3][4] [5] [6]  [7]   [8]   [9]  [10]    [11]
# [1]  1, 1, 1, 1, 1,  1,  1,   1,    1,    1,    1,       1  [A000012]
# [2]  1, 1, 1, 2, 5, 16, 61, 272, 1385, 7936, 50521, 353792  [A000111]
# [3]  1, 1, 1, 1, 3,  9, 19,  99,  477, 1513, 11259,  74601  [A178963]
# [4]  1, 1, 1, 1, 1,  4, 14,  34,   69,  496,  2896,  11056  [A178964]
# [5]  1, 1, 1, 1, 1,  1,  5,  20,   55,  125,   251,   2300  [A181936]
# [6]  1, 1, 1, 1, 1,  1,  1,   6,   27,   83,   209,    461  [A250283]

doc"""
Return the generalized André numbers which are the ``m``-alternating
permutations of length ``n``. Cf. A181937.

julia> [André(2, n) for n in 0:10]
[1, 1, -1, -2, 5, 16, -61, -272, 1385, 7936, -50521]
"""
André(m, n) = A(m, n)
# Currently you cannot document a memoized function.
# This is the only reason for this workaround.
@memoize function A(m, n)
    n ≤ 0 && return ZZ(1)
    # r = [m*k for k in 0:div(n-1,m)]
    r = range(0, m, div(n + m - 1, m))
    S = sum(binom(n, k) * A(m, k) for k in r)
    n % m == 0 ? -S : S
end

#@memoize function André(m::Int, N::Int)::Int
#    n = 5
#    n ≤ 0 && return ZZ(1)
#    # r = [m*k for k in 0:div(n-1,m)]
#    r = range(0, m, div(n + m - 1, m))
#    S = sum(binom(n, k) * A(m, k) for k in r)
#    n % m == 0 ? -S : S
#end

doc"""
Return the generalized André numbers which are the ``m``-alternating permutations
of length ``n``.

julia> [A181937(3, n) for n in 0:10]
[1, 1, 1, 1, 3, 9, 19, 99, 477, 1513, 11259]
"""
A181937(m::Int, n::Int) = abs(A(m, n))

doc"""
Return the up-down numbers (2-alternating permutations).

julia> A000111(10)
50521
"""
A000111(n::Int) = abs(André(2, n))

doc"""
Return the number of 3-alternating permutations.

julia> [A178963(n) for n in 0:10]
[1, 1, 1, 1, 3, 9, 19, 99, 477, 1513, 11259]
"""
A178963(n::Int) = abs(André(3, n))

doc"""
Return the number of 4-alternating permutations.

julia> [A178964(n) for n in 0:10]
[1, 1, 1, 1, 1, 4, 14, 34, 69, 496, 2896]
"""
A178964(n::Int) = abs(André(4, n))

doc"""
Return the number of 5-alternating permutations.

julia> [A181936(n) for n in 0:10]
[1, 1, 1, 1, 1, 1, 5, 20, 55, 125, 251]
"""
A181936(n::Int) = abs(André(5, n))

doc"""
Return the number of 6-alternating permutations.

julia> [A250283(n) for n in 0:10]
[1, 1, 1, 1, 1, 1, 1, 6, 27, 83, 209]
"""
A250283(n::Int) = abs(André(6, n))

doc"""
Generate the André numbers (a.k.a. Euler-up-down numbers A000111).
Don't confuse with the Euler numbers A122045.
"""
C000111() = Channel(csize=2) do c
    D = Dict{Int,fmpz}(0 => 1, -1 => 0)
    i = k = 0
    s = 1

    while true
        A = 0; D[k + s] = 0; s = -s
        for j in 0:i
            A += D[k]; D[k] = A; k += s
        end
        put!(c, A)
        i += 1
    end
end

end # module

module AndreTest
using Andre, SeqBase, Base.Test, SeqTests, OEISUtils, Nemo

function test()
    @testset "André" begin

        @test isa(André(2, 10), fmpz)

        @test André(2, 10) == -50521
        @test André(2, 50) == -6053285248188621896314383785111649088103498225146815121
        @test A178963(30) == 2716778010767155313771539
        @test A178964(40) == 11289082167259099068433198467575829

        if oeis_isinstalled()
            A = [A000111, A178963, A178964, A181936, A250283]
            SeqTest(A, 'A')
        end

        V = [1, 1, 1, 2, 5, 16, 61, 272, 1385, 7936, 50521]
        generator = C000111()
        for n in 1:10
            v = take!(generator)
            @test V[n]  == v
        end

        v = take!(generator)
        close(generator)

        @test isa(v, fmpz)
    end
end

function demo()
    for m in 1:8
        println([André(m, n) for n in 0:9])
    end

    generator = C000111()
    for n in 0:10
        v = take!(generator)
        println(v)
    end
    close(generator)
end

doc"""
for m in 1:20, n in 0:100 André(m,n) end :: 0.025871 seconds (98.44 k allocations: 2.320 MB)
"""
function perf()
    gc()
    @time (for m in 1:20, n in 0:100 André(m, n) end)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module GaussFactorial
using Nemo, NumberTheory, Products, SeqBase

export GaußFactorial, GaußLcm, A216919, A001783, L001783
export A124441, L124441, A124442, L124442, A066570, A160377, L160377
export A103131, L103131, A055634, A038610, L038610, A128247, L128247
export B033948, B033949, HasPrimitiveRoot, HasNoPrimitiveRoot

# Reminder: We write m ⊥ n or ⊥(m, n) for 'm prime to n'.
# Knuth, Graham and Patashnik write in "Concrete Mathematics":
# "Hear us, O mathematicians of the world! Let us not wait any longer! We can
# make many formulas clearer by defining a new notation now! Let us agree to
# write m ⊥ n, and to say "m is prime to n", if m and n are relatively prime."

doc"""
Return ``∏_{1 ≤ j ≤ N, j ⊥ n} j``, the product of the positive integers which
are ``≤ N`` and are prime to ``n``.
"""
GaußFactorial(N::Int, n::Int) = n == 0 ? 0 : ∏([j for j in 1:N if ⊥(j, n)])

doc"""
Return ``∏_{1 ≤ j ≤ N, j ⊥ n} j``, the product of the positive integers which
are ``≤ N`` and are prime to ``n``.
"""
A216919(N::Int, n::Int) = GaußFactorial(N, n)

doc"""
Return ``∏_{1 ≤ j ≤ n, j ⊥ n} j``, the product of the positive integers which
are ``≤ n`` and are prime to ``n``.
"""
A001783(n::Int) = GaußFactorial(n, n)

doc"""
Return a list of the first len terms of A001783.
"""
L001783(len::Int) = SeqArray(len, A001783)

doc"""
Return the product of the positive integers which are ``≤ n/2`` and are
prime to ``n``.
"""
A124441(n::Int) = GaußFactorial(div(n, 2), n)

doc"""
Return a list of the first len terms of A124441.
"""
L124441(len::Int) = SeqArray(len, A124441)

doc"""
Return the product of the positive integers which are ``≥ n/2`` and are
prime to ``n``.
"""
A124442(n::Int) = ∏([j for j in cld(n, 2):n if ⊥(j, n)])

doc"""
Return a list of the first len terms of A124442.
"""
L124442(len::Int) = SeqArray(len, A124442)

doc"""
Return the product of numbers ``≤ n`` that have a prime factor in common with ``n``.
"""
A066570(n::Int) = div(GaußFactorial(n, 1), GaußFactorial(n, n))

doc"""
The product of the residues in ``[1, n]`` relatively prime to n taken modulo n,
where the absolute representation of the integers modulo n is used.
A160377(n) = mod(GaußFactorial(n, n), n).
"""
function A160377(n::Int)
    r = 1
    for j in 1:n
        ⊥(j, n) && (r = mod(r * j, n))
    end
    r
end

doc"""
Return a list of the first len terms of A160377.
"""
L160377(len::Int) = SeqArray(len, A160377)

doc"""
Return the product of the residues in ``[1,n]`` which are prime to n taken
modulo n, where the symmetric (or minimal) representation of the integers
modulo n is used.
A103131(n) = mods(GaußFactorial(n, n), n).
"""
function A103131(n::Int)
    r = 1
    for j in 1:n
        ⊥(j, n) && (r = mods(r * j, n))
    end
    r
end

doc"""
Return a list of the first len terms of A103131.
"""
L103131(len::Int) = SeqArray(len, A103131)

doc"""
Return the 2-adic factorial function.
"""
A055634(n::Int) = (-1)^n * GaußFactorial(n, 2)

doc"""
Return lcm``_{1 ≤ j ≤ n, j ⊥ n} j``, the least common multiple of the positive
integers which are ``≤ n`` and are prime to ``n``.
"""
GaußLcm(N::Int, n::Int) = n == 0 ? 0 : lcm([j for j in 1:N if ⊥(j, n)])

doc"""
Return the least common multiple of integers less than and prime to n.
"""
A038610(n::Int) = GaußLcm(n, n)

doc"""
Return a list of the first len terms of A038610.
"""
L038610(len::Int) = SeqArray(len, A038610)

#A128247(n::Int) = div(GaußFactorial(n, n), GaußLcm(n, n))
doc"""
Return ∏(COP(n)) / lcm(COP(n)) where COP(n) denotes the positive integers which
are ``≤ n`` and are prime to ``n``.
"""
function A128247(n::Int)
    COP = [j for j in 1:n if ⊥(j, n)]
    COP == [] ? 0 : div(∏(COP), lcm(COP))
end

doc"""
Return a list of the first len terms of A128247.
"""
L128247(len::Int) = SeqArray(len, A128247)

# From Gauß's generalization of Wilson's theorem it follows that, for ``n>1``,
# mods(GaußFactorial``(n, n), n) = -1`` if and only if there exists a primitive
# root modulo ``n`` (cf. Hardy and Wright, Theorem 129).
doc"""
Query if there exists a primitive root modulo ``n``.
"""
HasPrimitiveRoot(n::Int) = 0 < n ≤ 2 ? true : mods(GaußFactorial(n, n), n) == -1

doc"""
Return numbers that have a primitive root (the multiplicative group
modulo ``n`` is cyclic).
"""
B033948(bound::Int) = FindUpTo(bound, HasPrimitiveRoot)

doc"""
Query if the discriminant of the n-th cyclotomic polynomial is a square.
"""
HasNoPrimitiveRoot(n::Int) = n ≤ 2 ? false : mod(GaußFactorial(n, n), n) == 1

doc"""
Return numbers ``n`` such that the discriminant of the n-th cyclotomic polynomial
is a square.
"""
B033949(bound::Int) = FindUpTo(bound, HasNoPrimitiveRoot)

end # module

module GaussFactorialTest
using Base.Test, Nemo, SeqBase, NumberTheory, Products, OEISUtils, GaussFactorial

function test()
    @testset "GaußFactorial" begin

        if oeis_isinstalled()

            B = [B033949, B033948]

            for seq in B
                name = SeqName(seq)
                # the parameter is not 'length' but 'search bound'.
                O = oeis_local(name, 12)
                S = seq(100)
                # println(name); println(O); println(S)
                @test all(S[0:9] .== O[0:9])
            end

            A = [A001783, A124441, A124442, A066570, A160377, A128247]

            for seq in A
                name = SeqName(seq)
                O = oeis_local(name, 10)
                # Our sequences start at offset 0 wheras in the OEIS at offset 1.
                S = SeqArray([seq(n) for n in 1:10])
                # println(name); println(O); println(S)
                @test all(S .== O)
            end

            L = [L124442, L124441, L128247, L160377, L001783, L038610, L103131]

            for seq in L
                name = SeqName(seq)
                O = oeis_local(name, 10)
                S = seq(11)
                # println(name); println(O); println(S)
                # Our sequences start at offset 0 wheras in the OEIS at offset 1.
                @test all(S[1:10] .== O[0:9])
            end
        end
    end
end

function demo()
    for n in 0:6   println(n, " ↦ ", L124442(n)) end
end

doc"""
 L001783(1000) :: 0.625464 seconds (2.44 M allocations: 44.045 MB, 35.76% gc time)
 [A055634(n) for n in 1:1000] :: 0.446586 seconds (2.36 M allocations: 41.835 MB, 26.88% gc time)
 [A103131(n) for n in 1:1000] :: 0.368846 seconds (2.00 M allocations: 30.556 MB, 35.05% gc time)
"""
function perf()
    gc()
    A001783(10)
    @time L001783(1000)
    gc()
    A055634(10)
    @time [A055634(n) for n in 1:1000]
    gc()
    A103131(10)
    @time [A103131(n) for n in 1:1000]

    # Remark on the implementation of A103131.
    # A(n::Int) = mods(GaußFactorial(n, n), n)
    # @time [A(n) for n in 1:2000]
    # 2.804732 seconds (10.19 M allocations: 182.562 MB, 38.02% gc time)
    # function B(n::Int)
    #    r = 1
    #    for j in 1:n
    #        ⊥(j, n) && (r = mods(r*j, n))
    #    end
    #    r
    # end
    # @time [B(n) for n in 1:2000]
    # 1.712230 seconds (8.00 M allocations: 122.147 MB, 40.57% gc time)

end

function main()
    test()
    demo()
    perf()
end

main()

end # module

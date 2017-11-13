# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Products
using Nemo, SeqBase

export ∏, Product, F!, RisingFactorial, ↑, FallingFactorial, ↓
export MultiFactorial, A000407, A124320, A265609, A000142
export A000165, A032031, A007559, A008544, A007696, A001813, A008545, A047053
export A081125, A001147

doc"""
If ``a ≤ b`` then return the product of ``i`` in ``a:b`` else return ``1``.
"""
function ∏(a::Int, b::Int)
    n = b - a
    if n < 24
        p = ZZ(1)
        for k in a:b
            p *= k
        end
        return p
    end
    m = div(a + b, 2)
    ∏(a, m) * ∏(m + 1, b)
end

doc"""
If ``a ≤ b`` then return the product of ``i`` in ``a:b`` else return ``1``.
"""
Product(a::Int, b::Int) = ∏(a, b)

doc"""
Return the accumulated product of an array.
"""
function ∏(A)
    function prod(a::Int, b::Int)
        n = b - a
        if n < 24
            p = ZZ(1)
            for k in a:b
                p *= A[k]
            end
            return p
        end
        m = div(a + b, 2)
        prod(a, m) * prod(m + 1, b)
    end
    A == [] && return 1
    if IsSeqArray(A)
        prod(0, length(linearindices(A))-1)
    else
        prod(1, length(A))
    end
end

doc"""
Return the accumulated product of an array.
"""
Product(A) = ∏(A)

doc"""
Return ``\frac{n!} {⌊n/2⌋!}``.
"""
A081125(n::Int) = ∏(div(n, 2) + 1, n)

doc"""
Return the rising factorial which is the product of ``i`` in ``n:(n + k - 1)``.
"""
RisingFactorial(n::Int, k::Int) = ∏(n, n + k - 1)

doc"""
Return the rising factorial which is the product of ``i`` in ``n:(n + k - 1)``.
A convenient infix syntax for the rising factorial is n ↑ k.
"""
↑(n, k) = RisingFactorial(n, k)

doc"""
Return the rising factorial i.e. the product of ``i`` in ``n:(n + k - 1)``.
"""
A265609(n::Int, k::Int) = RisingFactorial(n, k)

# *** deprecated, use n ↑ k instead ***
#doc"""
#Return 'Pochhammer(n, k)', which is ambiguous in the literature, as the
#RisingFactorial(n,k).
#"""
#Pochhammer(n::Int, k::Int) = ∏(n, n + k - 1)
# *** deprecated ***

doc"""
Return the falling factorial which is the product of ``i`` in ``(n - k + 1):n``.
"""
FallingFactorial(n::Int, k::Int) = ∏(n - k + 1, n)

doc"""
Return the falling factorial which is the product of ``i`` in ``(n - k + 1):n``.
A convenient infix syntax for the falling factorial is n ↓ k.
"""
↓(n, k) = FallingFactorial(n, k)

doc"""
Return the number of permutations of n letters, ``n! = ∏(1, n)``.
'' is a shortcut for 'Factorial'.
"""
F!(n::Int) = Nemo.fac(n)

# function !(n) F!(n) end

doc"""
Return the factorial numbers.
"""
A000142(n::Int) = Nemo.fac(n)

doc"""
Return the central rising factorial ``(n+1) ↑ (n+1) = (2n+1)! / n!``.
"""
A000407(n::Int) = (n + 1) ↑ (n + 1)
# A000407(n::Int) = ∏(n + 1, 2n + 1)

doc"""
Return the restricted rising factorial which is zero for ``n < 0`` or ``k > n``.
"""
A124320(n::Int, k::Int) = (n < 0 || k > n) ? 0 : ∏(n, n + k - 1)

doc"""
Return the multi-factorial which is the function ``n → ∏(a + b, a(n-1) + b)``
"""
MultiFactorial(a::Int, b::Int) = n -> ∏([a * k + b for k in 0:(n - 1)])

doc"""
Return the double factorial of odd numbers, ``1×3×5×...×(2n-1) = (2n-1)!!``.
"""
A001147(n::Int) = MultiFactorial(2, 1)(n)

doc"""
Return the double factorial of even numbers: ``2^n n! = (2n)!!``.
"""
A000165(n::Int) = MultiFactorial(2, 2)(n)

doc"""
Return the triple factorial numbers with shift 1, ``3^n n! = (3n)!!!``.
"""
A007559(n::Int) = MultiFactorial(3, 1)(n)

doc"""
Return the triple factorial numbers with shift 2.
"""
A008544(n::Int) = MultiFactorial(3, 2)(n)

doc"""
Return the triple factorial numbers with shift 3.
"""
A032031(n::Int) = MultiFactorial(3, 3)(n)

doc"""
Return the quadruple factorial numbers with shift 1.
"""
A007696(n::Int) = MultiFactorial(4, 1)(n)

doc"""
Return the quadruple factorial numbers with shift 2, ``(2n)!/n!``.
"""
A001813(n::Int) = MultiFactorial(4, 2)(n) # = ∏(n + 1, 2n)

doc"""
Return the quadruple factorial numbers with shift 3.
"""
A008545(n::Int) = MultiFactorial(4, 3)(n)

doc"""
Return the quadruple factorial numbers ``4^n n!``.
"""
A047053(n::Int) = MultiFactorial(4, 4)(n)

end # module

# RisingFactorial(n, k)
# n\k [0  1   2    3     4      5        6         7          8]
# --------------------------------------------------------------
# [0] [1, 0,  0,   0,    0,     0,       0,        0,         0]
# [1] [1, 1,  2,   6,   24,   120,     720,     5040,     40320]
# [2] [1, 2,  6,  24,  120,   720,    5040,    40320,    362880]
# [3] [1, 3, 12,  60,  360,  2520,   20160,   181440,   1814400]
# [4] [1, 4, 20, 120,  840,  6720,   60480,   604800,   6652800]
# [5] [1, 5, 30, 210, 1680, 15120,  151200,  1663200,  19958400]
# [6] [1, 6, 42, 336, 3024, 30240,  332640,  3991680,  51891840]
# [7] [1, 7, 56, 504, 5040, 55440,  665280,  8648640, 121080960]
# [8] [1, 8, 72, 720, 7920, 95040, 1235520, 17297280, 259459200]

module ProductsTest
using Base.Test, OEISUtils, SeqBase, SeqTests, Products, NumberTheory

function test()
    @testset "FallingFact" begin
        @test FallingFactorial(100, 100) == factorial(BigInt(100))
        @test (100 ↓ 100) == factorial(BigInt(100))
        @test FallingFactorial(333, 333) == factorial(BigInt(333))
        @test FallingFactorial(111, 0) == 1
    end
    @testset "RisingFact" begin
        @test RisingFactorial(11, 11) == 14079294028800
        @test (11 ↑ 11) == 14079294028800
        @test RisingFactorial(33, 33) == 31344295059422473624824839739793024055460338073600000000
        @test RisingFactorial(111, 0) == 1
    end
    @testset "MultiFact" begin
        a = SeqArray([MultiFactorial(2, 1)(n) for n in 0:6])
        b = SeqArray([1, 1, 3, 15, 105, 945, 10395])
        @test all(a .== b)
    end
    if oeis_isinstalled()
        A = [A000142, A000165, A007696, A001813, A047053, A001147, A008545,
            A081125, A000407, A032031, A007559, A008544]

            @testset "Products" begin
            SeqTest(A, 'A')
        end
    end
end

function demo()
    for n in 0:6
        println([RisingFactorial(n, k) for k in 0:7])
    end

    for n in 0:6
        println([FallingFactorial(n, k) for k in 0:7])
    end

    println(∏([FallingFactorial(n, n) for n in 0:12]))
    println(Product([FallingFactorial(n, n) for n in 0:12]))

    for n in 0:9, k in 0:9 println(n ↑ k) end
    for n in 0:9, k in 0:9 println(n ↓ k) end

end

doc"""
for n in 1:10000 F!(n) end :: 1.561170 seconds (10.00 k allocations: 156.250 KB)
for n in 1:1000 A000407(n) end :: 0.173578 seconds (995.90 k allocations: 15.196 MB)
for n in 1:200, k in 1:200 Pochhammer(n, k) end :: 1.855906 seconds (4.43 M allocations: 67.596 MB, 19.49% gc time)
"""
function perf()
    gc()
    @time (for n in 1:10000 F!(n) end)
    gc()
    @time (for n in 1:1000 A000407(n) end)
    gc()
    @time (for n in 1:200, k in 1:200 RisingFactorial(n, k) end)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module SwingingFactorial
using Nemo, SeqBase, Products, PrimeSieve, NumberTheory

export Swing, CatalanNumber, CentralBinomial, Factorial, Apéry
export A163590, A001790, A001803, A056040, A000984, A002457, A057977
export A281594, A080397, A000108, A163641, A049606, A005430
export L163085

SwingOddpart = fmpz[1,1,1,3,3,15,5,35,35, 315, 63, 693, 231,
   3003, 429, 6435, 6435, 109395,12155,230945,46189,969969,
   88179,2028117, 676039,16900975,1300075,35102025,5014575,
   145422675,9694845,300540195,300540195]

doc"""
Return the odd part of the swinging factorial ``n≀``. Cf. A163590.
"""
function swing_oddpart(n)
    n < 33 && return SwingOddpart[n+1]

    sqrtn = isqrt(n)
    factors = Primes(div(n, 2) + 1, n)
    r = Primes(sqrtn + 1, div(n, 3))
    s = filter(x -> IsOdd(div(n, x)), r)
    append!(factors, s)

    for prime in Primes(3, sqrtn)
        p, q = 1, n
        while true
            q = div(q, prime)
            q == 0 && break
            IsOdd(q) && (p *= prime)
        end
        p > 1 && push!(factors, p)
    end

    ∏(factors)
end

doc"""
Return the swinging factorial (a.k.a. Swing numbers) ``n≀``. Cf. A056040.
"""
Swing(n) = swing_oddpart(n) << count_ones(div(n, 2))

doc"""
Return the odd part of the swinging factorial ``n≀``
"""
A163590(n) = swing_oddpart(n)

doc"""
Return the odd part of the central binomial coefficient.
"""
A001790(n) = swing_oddpart(2n)

doc"""
Return the odd part of the Apéry numbers.
"""
A001803(n) = swing_oddpart(2n + 1)

doc"""
Return the swinging factorial (a.k.a. Swing numbers) ``n≀``.
"""
A056040(n) = Swing(n)

doc"""
Return the central binomial coefficient.
"""
A000984(n) = Swing(2n)

doc"""
Return the central binomial coefficient.
"""
CentralBinomial(n) = Swing(2n)

doc"""
Return the n-th Apéry number, ``n \binom{2n}{n}``.
"""
Apéry(n) = n * Swing(2n)

doc"""
Return the n-th Apéry numbers, ``n \binom{2n}{n}``.
"""
A005430(n) = Apéry(n)

doc"""
Return ``(2n+1)! / n!^2``.
"""
A002457(n) = Swing(2n + 1)

doc"""
Return the squarefree kernel of the swinging factorial.
"""
A163641(n) = Radical(Swing(n))

doc"""
Return the squarefree kernel of the central binomial coefficient.
"""
A080397(n) = Radical(Swing(2n))

doc"""
Return the cumulative product of the swinging factorial.
"""
function L163085(len)
    L = SeqArray(len)
    L[0] = 1
    for n in 1:len-1
        L[n] = L[n-1] * Swing(n)
    end
    L
end

doc"""
Return the n-th extended Catalan number. Cf. the exposition
[Lost Catalan Numbers](http://oeis.org/wiki/User:Peter_Luschny/TheLostCatalanNumbers).
"""
A057977(n) = div(Swing(n), div(n, 2) + 1)

doc"""
Return the n-th Catalan number.
"""
CatalanNumber(n) = div(Swing(2n), n + 1)

doc"""
Return the n-th Catalan number.
"""
A000108(n) = CatalanNumber(n)

doc"""
Return the squarefree kernel of the Catalan number.
"""
A281594(n) = Radical(CatalanNumber(n))

# Fast factorial http://luschny.de/math/factorial/FastFactorialFunctions.htm

FactorialOddPart = fmpz[1, 1, 1, 3, 3, 15, 45, 315, 315, 2835, 14175, 155925,
    467775, 6081075, 42567525, 638512875, 638512875, 10854718875, 97692469875,
    1856156927625, 9280784638125, 194896477400625, 2143861251406875,
    49308808782358125, 147926426347074375, 3698160658676859375]

doc"""
Return the largest odd divisor of ``n!``. Cf. A049606.
"""
function factorial_oddpart(n)
    n < length(FactorialOddPart) && return ZZ(FactorialOddPart[n + 1])
    swing_oddpart(n)*(factorial_oddpart(div(n, 2))^2)
end

doc"""
Return the factorial ``n! = 1×2×...×n``. Alternatively you can call F!(n)
which uses the Nemo implementation. Cf. A000142.
"""
function Factorial(n)
    n < 0 && ArgumentError("Argument must be ≥ 0")
    sh = n - count_ones(n)
    factorial_oddpart(n) << sh
end

doc"""
Return the odd part of the factorial.
"""
A049606(n) = factorial_oddpart(n)

end # module

module SwingingFactorialTest
using Base.Test, SeqTests, SeqBase, Nemo, Products, OEISUtils
using NumberTheory, SwingingFactorial

function test()
    @testset "Swing" begin

        @test isa(Swing(30), fmpz)

        @test Swing(0)  == 1
        @test Swing(27) == 280816200
        @test Swing(77) == 530731789949381124304200
        @test Swing(98) == 25477612258980856902730428600
        @test Swing(1000) == binom(1000, 500)

        N = Int(2)^10
        @test Swing(N) == binom(N, N>>1)

        a = L163085(9)
        b = SeqArray([1, 1, 2, 12, 72, 2160, 43200, 6048000, 423360000])
        @test all(a .== b)

        if oeis_isinstalled()

            A = [ A163590, A001790, A001803, A056040, A000984, A002457,
            A281594, A080397, A000108, A163641, A049606, A005430, A057977 ]
            SeqTest(A, 'A')
        end
    end
end

function demo()
    println([A163590(n) for n in 0:10])
    println([A001790(n) for n in 0:10])
    println([A001803(n) for n in 0:10])
    println([A056040(n) for n in 0:10])
    println([A000984(n) for n in 0:10])
    println([A002457(n) for n in 0:10])
    println([A080397(n) for n in 0:10])
    println([A000108(n) for n in 0:10])
    println([A281594(n) for n in 0:10])
end

doc"""
for n in 0:999 Swing(n) end :: 0.034996 seconds (106.70 k allocations: 4.163 MiB)
for n in 0:999 A000984(n) end :: 0.071966 seconds (208.77 k allocations: 7.221 MiB)
"""
function perf()
    gc()
    @time (for n in 0:999 Swing(n) end)
    @time (for n in 0:999 A000984(n) end)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

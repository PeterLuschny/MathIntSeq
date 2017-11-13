# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module NumberTheory
using Memoize, Nemo, SeqBase, Products
import Nemo.zero

export τ, σ, ϕ, ω, Ω, ⊥, ⍊, Divisors, PrimeDivisors, Factors, Radical, mods
export A000005, A000010, A000203, A001222, A001221, A008683, A181830, A034444
export L003277, A061142, A034386, A002110, Divides, IsPrime, IsCyclic, IsOdd
export IsPrimeTo, IsStrongPrimeTo, IsNonnegative, IsPositive, IsEven, IsSquare
export IsComposite, IsSquareFree, IsPrimePower, IsPowerOfPrimes, IsPerfectPower
export PrimeList

zero(::Nemo.fmpz) = fmpz(0)

doc"""
Return true if n is prime false otherwise.
"""
function IsPrime(n) Nemo.isprime(fmpz(n)) end
#function IsPrime(n)::Bool Nemo.isprime(fmpz(n)) end

doc"""
Return factors of ``n``.
"""
Factors(n) = Nemo.factor(fmpz(n))

doc"""
Return the positive integers dividing ``n``.
"""
function Divisors(n)
    n == fmpz(0) && return fmpz[]
    IsPrime(n) && return [fmpz(1), n]
    d = [fmpz(1)]
    for (p, e) in Factors(n)
        c = [p^i for i in 0:e]
        d = d * c'
        d = reshape(d, length(d))
    end
    sort!(d)
end

doc"""
Return the prime numbers dividing ``n``.
"""
function PrimeDivisors(n)
    IsPrime(n) && return [fmpz(n)]
    f = Factors(n)
    sort!([p for (p, e) in f])
end

doc"""
Return the radical of n which is the product of the prime numbers dividing ``n``
(also called the squarefree kernel of ``n``).
"""
Radical(n) = ∏(PrimeDivisors(n))

doc"""
Return ``Ω(n)``,  the number of prime divisors of ``n`` counted with multiplicity.
Cf. A001222.
"""
function Ω(n)
    n == fmpz(0) && return 0
    IsPrime(n) && return fmpz(1)
    f = Factors(n)
    sum([e for (p, e) in f])
end

doc"""
Return the number of prime divisors of ``n`` counted with multiplicity.
"""
A001222(n) = Ω(n)

doc"""
Return the result of replacing each prime factor of n with 2.
"""
A061142(n) = 1 << Int(Ω(n))

doc"""
Return ``ω(n)``,  the number of distinct prime divisors of ``n``. Cf. A001221.
"""
ω(n) = fmpz(length(PrimeDivisors(n)))

doc"""
Return the number of distinct prime divisors of ``n``.
"""
A001221(n) = ω(n)

doc"""
Return the number of unitary divisors of ``n``, ``d`` such that ``d`` divides ``n``
and ``d ⊥ n/d``.
"""
A034444(n::Int) = 1 << Int(ω(n))

doc"""
Return ``τ(n)`` (or ``σ_0(n)``), the number of divisors of ``n``. Cf. A000005.
"""
τ(n) = Nemo.sigma(fmpz(n), 0)

doc"""
Return the number of divisors of ``n``.
"""
A000005(n) = τ(n)

doc"""
Return ``σ(n)`` (or ``σ_1(n)``), the sum of the divisors of ``n``. Cf. A000203.
"""
σ(n) = Nemo.sigma(fmpz(n), 1)

doc"""
Return the Euler totient ``ϕ(n)``, numbers which are ``≤ n`` and prime to ``n``.
"""
ϕ(n) = Nemo.eulerphi(fmpz(n))

doc"""
Return the number of integers ``≤ n`` and prime to ``n``.
"""
A000010(n) = ϕ(n)

doc"""
Return the value of the Möbius function ``μ(n)`` which is the sum of the
primitive n-th roots of unity.
"""
μ(n) = Nemo.moebiusmu(fmpz(n))

doc"""
Return the value of the Möbius function ``μ(n)`` which is the sum of the
primitive n-th roots of unity.
"""
A008683(n) = μ(n)

doc"""
Return the sum of the divisors of ``n``.
"""
A000203(n) = σ(n)

doc"""
Query if ``n`` is primeto to ``k``.
"""
IsPrimeTo(n, k) = Nemo.gcd(fmpz(n), fmpz(k)) == fmpz(1)

doc"""
Query if ``m`` is primeto to ``n``.

Knuth, Graham and Patashnik write in "Concrete Mathematics":
"Hear us, O mathematicians of the world! Let us not wait any longer!
We can make many formulas clearer by defining a new notation now!
Let us agree to write m ⊥ n, and to say "m is prime to n", if m and n are
relatively prime."
"""
⊥(m, n) = IsPrimeTo(m, n)

doc"""
Query if ``n`` is strong prime to ``k``.
"""
IsStrongPrimeTo(n, k) = IsPrimeTo(n, k) && k ∉ Divisors(n - 1)

doc"""
Query if ``n`` is strong prime to ``k``.
"""
⍊(n, k) = IsStrongPrimeTo(n, k)

doc"""
Return the number of integers ``≤ n`` which are strong primeto to ``n``.
"""
A181830(n) = n == 0 ? 0 : ϕ(n) - τ(n - 1)

doc"""
Is ``n`` a cyclic number? ``n`` such that there is just one group of order ``n``.
"""
IsCyclic(n) = n == 0 ? false : ⊥(n, ϕ(n))

doc"""
Return list of cyclic numbers of length len.
"""
L003277(len) = List(len, IsCyclic)

doc"""
Return the least absolute remainder mods uses the symmetric representation for
integers modulo m, i.e. remainders will be reduced to integers in the range
``[-``div``(|m| - 1, 2),``div``(|m|, 2)]``.
"""
function mods(b, a)
    b == 0 && return a
    h = a >> 1
    (q, r) = divrem(b, a)
    if h <  r  r -= a end
    if h < -r  r += a end
    r
end

doc"""
Is the integer ``n`` nonnegative?
"""
IsNonnegative(n) = n ≥ 0

doc"""
Is the integer ``n`` positive?
"""
IsPositive(n) = n > 0

doc"""
Is the integer ``n`` a square number?
"""
IsSquare(n) = Nemo.issquare(fmpz(n))

doc"""
Is the integer ``n`` a composite number?
"""
IsComposite(n) = 1 < Ω(n) && (n > 0)

doc"""
Is the integer ``n`` a squarefree number?
"""
IsSquareFree(n) = ω(n) == Ω(n) && (n > 0)

doc"""
Is the integer ``n`` a prime power?
"""
IsPrimePower(n) = ω(n) == 1

doc"""
Is the integer ``n`` a power of primes?
"""
IsPowerOfPrimes(n) = (n == 1) || (ω(n) == 1)

doc"""
Is the integer ``n`` a perfect powers?
"""
IsPerfectPower(n) = ω(n) == 1 && Ω(n) ≠ 1

doc"""
Return `true` if b is divisible by a, otherwise return `false`.
"""
Divides(a, b) = a ≠ 0 && rem(fmpz(b), fmpz(a)) == fmpz(0)

# Defined in the module Base.
doc"""
Is n divisble by 2?
"""
IsEven(n) = Base.iseven(n)
# n % 2 == 0 ? true : false

doc"""
Is n indivisble by 2?
"""
IsOdd(n) = Base.isodd(n)
# n % 2 != 0 ? true : false

doc"""
Return the primorial of n, the product of the primes ≤ n.
"""
A034386(n) = Nemo.primorial(n)

doc"""
Return a list of the first n primes.
"""
PrimeList(len::Int) = SeqArray(len, IsPrime)

doc"""
Return the product of first n primes.
"""
A002110(n) = ∏(PrimeList(n))

# In the module GaussFactorial are the definitions of
# HasPrimitiveRoot
# HasNoPrimitiveRoot

# In the module Abundant
# IsAbundant

# Further indicators, less suited for computations.

# isA008578(n::Int) = all(⊥(k, n)     for k in 1:n-1)
# isA002182(n::Int) = all(τ(k) < τ(n) for k in 1:n-1)
# isA002110(n::Int) = all(ω(k) < ω(n) for k in 1:n-1)
# isA131577(n::Int) = all(Ω(k) < Ω(n) for k in 1:n-1)

end # module

module NumberTheoryTest
using Base.Test, OEISUtils, SeqBase, SeqTests, Nemo, NumberTheory, Products

function test()

    # 0-based version of sequences
    Data = Dict{Int, Array{fmpz}}(
    034386 => [1, 1, 2, 6, 6, 30, 30, 210, 210, 210],
    061142 => [1, 1, 2, 2, 4, 2, 4, 2, 8, 4],
    002110 => [1, 2, 6, 30, 210, 2310, 30030, 510510, 9699690, 223092870],
    000005 => [0, 1, 2, 2, 3, 2, 4, 2, 4, 3],
    000010 => [0, 1, 1, 2, 2, 4, 2, 6, 4, 6],
    000203 => [0, 1, 3, 4, 7, 6, 12, 8, 15, 13],
    001222 => [0, 0, 1, 1, 2, 1, 2, 1, 3, 2],
    001221 => [0, 0, 1, 1, 1, 1, 2, 1, 1, 1],
    008683 => [0, 1, -1, -1, 0, -1, 1, -1, 0, 0],
    181830 => [0, 1, 0, 0, 0, 1, 0, 2, 2, 2],
    034444 => [1, 1, 2, 2, 2, 2, 4, 2, 2, 2]
    )

    @testset "NumTheory" begin
        @test τ(7560) == 64
        @test τ(46080) == 66
        @test τ(25920) == 70

        @test σ(7560) == 28800
        @test σ(46080) == 159666
        @test σ(25919) == 25920
        @test σ(25920) == 92202

        @test ω(7560) == 4
        @test ω(46080) == 3
        @test ω(25919) == 1
        @test ω(25920) == 3

        @test Ω(7560) == 8
        @test Ω(46080) == 13
        @test Ω(25919) == 1
        @test Ω(25920) == 11

        @test Radical(58564) == 22
        @test Radical(58565) == 58565

        FB(n::Int) = (r = 1; for k in 1:n ⊥(n, k) && (r = mods(r * k, n)) end; r)
        FA(n::Int) = mods(∏([j for j in 1:n if ⊥(j, n)]), n)

        for n in 1:20
            @test FA(n) == FB(n)
        end

        if oeis_isinstalled()
            A = [A034386, A061142, A002110, A000005, A000010, A000203,
                A001222, A001221, A008683, A181830, A034444]

            for a in A
                S = SeqArray([a(i) for i in 0:9])
                anum = SeqNum(a)
                data = SeqArray(Data[anum])
                # println(anum); println(S); println(data)
                @test all(S[0:9] .== data[0:9])
            end

            L = [L003277]
            SeqTest(L, 'L')
        end
    end

    composita = [false, false, false, false, true, false, true, false]
    @testset "Queries" begin
        for n in 0:7
            @test IsComposite(n) == composita[n + 1]
        end
    end
end

function demo()
    for n in 390:400
        println(n, " ---")
        println(Factors(n))
        println(Divisors(n))
        println(τ(n), ", ", σ(n))
        println(PrimeDivisors(n))
        println(Radical(n))
    end
end

doc"""
[Divisors(n) for n in 1:10000] :: 0.124251 seconds (480.92 k allocations: 26.297 MB, 14.44% gc time)
[Radical(n)  for n in 1:10000] :: 0.015362 seconds (62.13 k allocations: 6.409 MB)
Timing is unstable, some chaching somewhere?
[Divisors(n) for n in 1:10000] :: 0.281049 seconds (1.87 M allocations: 75.677 MB)
[Radical(n)  for n in 1:10000] :: 0.081750 seconds (250.84 k allocations: 15.650 MB)
"""
function perf()
    @time [Divisors(n) for n in 1:10000]
    @time [Radical(n)  for n in 1:10000]
end

function main()
    test()
    demo()
    perf()
end

main()

 end # module

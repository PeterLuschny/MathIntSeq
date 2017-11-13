# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Counts
using Nemo, SeqBase, NumberTheory

export L000961, L002808, L005117, L013928, L025528, L065515
export L065855, L069637, L246547, L246655, L000720
export A007917, A151800, A257993
export PreviousPrime, NextPrime, PrimePiList

doc"""
Return a list of composite numbers of length len.
(Numbers which have more than one prime divisor.)

julia> L002808(8)
"""
L002808(len) = List(len, IsComposite)

doc"""
Return a list of the number of composite numbers ``≤ n``.

julia> L065855(8)
"""
L065855(len) = CountList(len, IsComposite)

doc"""
Return a list of squarefree numbers of length len.
(Numbers which are not divisible by a square greater than 1.)

julia> L005117(8)
"""
L005117(len) = List(len, IsSquareFree)

doc"""
Return a list of the number of squarefree numbers ``< n``.

julia> L013928(8)
"""
L013928(len) = CountList(len, IsSquareFree)

doc"""
Return a list of powers of primes of length len.
(Numbers of the form ``p^k`` where ``p`` is a prime and ``k ≥ 0``.)

julia> L000961(8)
"""
L000961(len) = List(len, IsPowerOfPrimes)

doc"""
Return the number of powers of primes ``≤ n``. (Powers of primes are numbers
of the form ``p^k`` where ``p`` is a prime and ``k ≥ 0``.)

julia> L065515(8)
"""
L065515(len) = CountList(len, IsPowerOfPrimes)

doc"""
Return a list of prime powers of length len.
(Numbers of the form ``p^k`` where ``p`` is a prime and ``k ≥ 1``.)

julia> L246655(8)
"""
L246655(len) = List(len, IsPrimePower)

doc"""
Return a list of the number of prime powers ``≤ n`` with exponents ``k ≥ 1``.

julia> L025528(8)
"""
L025528(len) = CountList(len, IsPrimePower)

doc"""
Return a list of perfect powers of length len.
(Numbers of the form ``p^k`` where ``p`` is a prime and ``k ≥ 2``.

julia> L246547(8)
"""
L246547(len) = List(len, IsPerfectPower)

doc"""
Return a list of the number of prime powers ``≤ n`` with exponents ``k ≥ 2``.

julia> L069637(8)
"""
L069637(len) = CountList(len, IsPerfectPower)

# cf. also:
# A067535 Smallest squarefree number >= n.
# A070321 Largest squarefree number <= n.
# A025528 Number of prime powers <= n with exponents > 0.
# A000015 Smallest prime power >= n.
# A031218 Largest prime power <= n
# A167184 Smallest prime power >= n that is not prime.
# A081676 Largest perfect power <= n

doc"""
Return the largest prime in ``N`` (the semiring of natural numbers including zero)
less than n for ``n ≥ 0``. (The "prev_prime" function of Mathematica, Maple,
Magma and SageMath.)
"""
A007917(n::Int) = Previous(n, IsPrime)

doc"""
Return the largest prime in ``Z`` (the ring of all integers)
less than ``n`` for ``n ≥ 0``.
"""
PreviousPrime(n::Int) = n ∈ [0, 1, 2] ? -2 : Previous(n - 1, IsPrime)

doc"""
Return least prime ``> n``. (The "next_prime" function of Mathematica, Maple,
Magma and SageMath.)
"""
NextPrime(n::Int) = Next(n, IsPrime)

doc"""
Return least prime ``> n``. (The "next_prime" function of Mathematica, Maple,
Magma and SageMath.)
"""
A151800(n::Int) = Next(n, IsPrime)

doc"""
Return the list of number of primes ``≤ n`` for ``n ≥ 0``.

julia> PrimePiList(8)
[0, 0, 1, 2, 2, 3, 3, 4]
"""
PrimePiList(len::Int) = CountList(len, IsPrime)

doc"""
Return the list of number of primes ``≤ n`` for ``n ≥ 0``.

julia> L000720(8)
[0, 0, 1, 2, 2, 3, 3, 4]
"""
L000720(len::Int) = PrimePiList(len)

doc"""
Return the index of the least prime not dividing n.
"""
function A257993(n::Int)
    c, p = 1, 2
    while n % p == 0
        p = NextPrime(p)
        c += 1
    end
    c
end

end # module

module CountsTest
using Base.Test, SeqBase, SeqTests, Counts, NumberTheory, OEISUtils, Nemo

function test()

    indicators = [IsPositive, IsEven, IsSquare, IsPrime]
    indicatorNames = ["IsPositive", "IsEven", "IsSquare", "IsPrime"]

    len = 14
    @testset "Counts" begin

        @test Nth(96, IsPrime) == 503
        @test Nth(97, IsPrime) == 509
        @test Nth(98, IsPrime) == 521

        # In other words: 97 is the 25-th prime.
        @test Count(96, IsPrime) == 24
        @test Count(97, IsPrime) == 25
        @test Count(98, IsPrime) == 25

        @test Last(List(24, IsPrime)) == 89
        @test Last(List(25, IsPrime)) == 97

        for (i, isA) in enumerate(indicators)
            # This test shows that the logic behind 'Nth' and 'Count' is OK.
            for n in 1:len
                @test isA(n) == (Nth(Count(n, isA), isA) == n)
                @test     n  ==  Count(Nth(n, isA), isA)
            end
        end

        a = [A257993(n) for n in 1:10]
        b = [1, 2, 1, 2, 1, 3, 1, 2, 1, 2]
        @test all(a .== b)

        if oeis_isinstalled()
            L = [
            L000961, L002808, L005117, L013928, L246547, L246655
            # L025528, L065515, L065855, L069637, L000720
            # Since the OEIS offset is 1 and our offset is 0 there is an
            # additional 0 in our version of these sequences. Everything OK here
            # but as usual references to the OEIS are - as always - approximative.
            ]
            SeqTest(L, 'L')
        end
    end
end

function demo()

    indicators = [IsNonnegative, IsPositive, IsEven, IsComposite,
    IsSquare, IsSquareFree, IsPrimePower, IsPowerOfPrimes,
    IsPerfectPower, IsPrime
    ]
    indicatorNames = ["IsNonnegative", "IsPositive", "IsEven", "IsComposite",
    "IsSquare", "IsSquareFree", "IsPrimePower", "IsPowerOfPrimes",
    "IsPerfectPower", "IsPrime"
    ]

    len = 14
    for (i, isA) in enumerate(indicators)

        println("Predicate      ", indicatorNames[i])
        println("First          ", First(isA))
        println("---")

        println("Nth            ", [Nth(n, isA) for n in 0:len])
        println("List           ", List(len, isA))
        println("FindUpTo       ", FindUpTo(len, isA))
        println("IterateUpTo    ", [k for k in IterateUpTo(len, isA)])
        println("---")

        println("IndicatorsFind ", IndicatorsFind(len, isA))
        println("Indicators     ", [k for k in Indicators(len, isA)])
        println("---")

        println("IndexIn (list) ", [IndexIn(fmpz(n), List(len, isA)) for n in 0:len])

        println("Count   (list) ", [Count(n, isA) for n in 0:len - 1])
        println("CountList      ", CountList(len, isA))

        println("Previous(n)    ", [Previous(n, isA) for n in 0:len])
        println("Next(n)        ", [Next(n, isA)  for n in 0:len])
        println("Nth(Count(n))  ", [Nth(Count(n, isA), isA) for n in 0:len])
        println("Count(Nth(n))  ", [Count(Nth(n, isA), isA) for n in 0:len])

        println("---------------------------")
    end
end

doc"""
[A257993(n) for n in 1:10000] :: 0.001372 seconds (16.18 k allocations: 331.016 KB)
PrimePiList(10000) :: 0.004299 seconds (42.15 k allocations: 736.781 KiB)
"""
function perf()
    @time [A257993(n) for n in 1:10000]
    gc()
    @time PrimePiList(10000)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

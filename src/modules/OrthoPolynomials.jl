# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module OrthoPolynomials
using Nemo, SeqBase

export OrthoPoly, InvOrthoPoly
export T053121, T216916, T217537, T064189, T202327, T111062, T026300, T099174
export T066325, T049310, T137338, T104562, T037027, T049218, T159834, T137286
export T053120, T053117, T111593, T059419
export L217924, L005773, L108624, L005425, L000085, L001464, L003723, L006229

# Cf. http://oeis.org/wiki/User:Peter_Luschny/AignerTriangles

doc"""
By the theorem of Favard an orthogonal polynomial systems  ``p_{n}(x)`` is a
sequence of real polynomials with deg``(p_{n}(x)) = n`` for all ``n`` if and only if
``$ p_{n+1}(x) = (x - s_n)p_n(x) - t_n p_{n-1}(x) $``
with ``p_{0}(x)=1`` for some pair of seq's ``s_k`` and ``t_k``. Return the
coefficients of the polynomials as a triangular array with dim rows.
"""
function OrthoPoly(dim::Int, s::Function, t::Function)
    T = SeqTriangle(dim)
    T[0] = 1
    lo = hi = 0

    for row in 2:dim
        olo = lo
        lo = hi + 1
        hi += row
        j = lo - row

        for k in lo:hi - 1
            T[k] = (j < olo ? 0 : T[j]) + s(k - 1 - lo) * T[j + 1] + (j + 2 < lo ? t(k - lo) * T[j + 2] : 0)
            j += 1
        end
        T[hi] = 1
    end
    T
end

doc"""
Return the inverse of the coefficients of the orthogonal polynomials generated
by ``s`` and ``t`` as a triangular array with dim rows.
"""
function InvOrthoPoly(dim::Int, s::Function, t::Function)
    dim ≤ 0 && return fmpz[]
    T = zeros(ZZ, dim, dim)
    #T = SeqMatrix(dim)
    for n in 1:dim T[n, n] = 1 end
    for n in 1:dim - 1
        for k in 1:n + 1
            T[n + 1, k] = (k > 1 ? T[n,k - 1] : 0) - s(n - 1) * T[n,k] - (n > 1 ? t(n - 1) * T[n - 1,k] : 0)
        end
    end
    # TODO eliminate □toT.
    SeqArray(□toΔ(T))
end

doc"""
Return the Catalan triangle (with 0's) read by rows.
"""
T053121(dim::Int) = OrthoPoly(dim, n -> 0, n -> 1)

# doc"""
# binomial(n, floor(n/2)).
# """
# L001405(len::Int) = RowSums(T053121(len))

doc"""
Return the coefficients of some orthogonal polynomials related to set partitions
without singletons (A000296).
"""
T216916(dim::Int) = OrthoPoly(dim, n -> n + 2, n -> n + 2)

doc"""
Return the triangle ``T(n,k)`` of tangent numbers, coefficient of ``x^n/n!``
in the expansion of ``(\tan x)^k/k!``.
"""
T059419(dim::Int) = OrthoPoly(dim, n -> 0, n -> (n + 1) * (n + 2))

doc"""
Return the expansion of ``\exp(\tan(x))``.
"""
L006229(len::Int) = RowSums(T059419(len))

doc"""
Return the first len integers defined as ``a(n) = n! [x^n] \exp(2 \exp(x)-x-2)``.
"""
L217924(len::Int) = RowSums(T217537(len))

doc"""
Return the coefficients of some orthogonal polynomials related to indecomposable
set partitions without singletons (A098742).
"""
T217537(dim::Int) = OrthoPoly(dim, n -> n + 1, n -> n + 1)

doc"""
Return the Motzkin triangle read in row reverse order.
"""
T064189(dim::Int) = OrthoPoly(dim, n -> 1, n -> 1)

doc"""
Return the Motzkin triangle.
"""
T026300(dim::Int) = RowReverse(OrthoPoly(dim, n -> 1, n -> 1))

doc"""
Return the number of directed animals of size n as an array of length len.
"""
L005773(len::Int) = RowSums(T064189(len))

doc"""
Return the coefficients of ``x^n`` in the expansion of ``((-1-x+√(1+2x+5x^2))/2)^k``
as a triangle with dim rows.
"""
T202327(dim::Int) = OrthoPoly(dim, n -> -1, n -> -1)

doc"""
Return the sequence with generating function satisfying
``x = (A(x)+(A(x))^2)/(1-A(x)-(A(x))^2)``.
"""
L108624(len::Int) = RowSums(T202327(len))

doc"""
Return the triangle T(n, k) = ``\binom{n}{k}`` involutions``(n - k)``.
"""
T111062(dim::Int) = OrthoPoly(dim, n -> 1, n -> n + 1)

doc"""
Return the number of self-inverse partial permutations.
"""
L005425(len::Int) = RowSums(T111062(len))

doc"""
Return the coefficients of the modified Hermite polynomials.
"""
T099174(dim::Int) = OrthoPoly(dim, n -> 0, n -> n + 1)

# Also
# T099174(dim::Int) = InvOrthoPoly(dim, n -> 0, n -> -n)

doc"""
Return the number of involutions.
"""
L000085(len::Int) = RowSums(T099174(len))

doc"""
Return the coefficients of unitary Hermite polynomials He``_n(x)``.
"""
T066325(dim::Int) = InvOrthoPoly(dim, n -> 0, n -> n)

doc"""
Return the sequence defined by ``a(n) = n! [x^n] \exp(-x-(x^2)/2)``.
"""
L001464(len::Int) = RowSums(T066325(len), true)

doc"""
Return the triangle of tanh numbers.
"""
T111593(dim::Int) = OrthoPoly(dim, n -> 0, n -> -n * (n + 1))

doc"""
Return the sequence defined by ``A(n) = n! [x^n] \exp{\tan{x}}`` as an array
of length len.
"""
L003723(len::Int) = RowSums(T111593(len))

doc"""
Return the coefficients of Chebyshev's U``(n, x/2)`` polynomials.
"""
T049310(dim::Int) = InvOrthoPoly(dim, n -> 0, n -> 1)

doc"""
Return the coefficients of the Charlier polynomials with parameter ``a=1``.
"""
T137338(dim::Int) = InvOrthoPoly(dim, n -> n + 1, n -> n + 1)

doc"""
Return the inverse of the Motzkin triangle (A064189).
"""
T104562(dim::Int) = InvOrthoPoly(dim, n -> 1, n -> 1)

doc"""
Return the skew Fibonacci-Pascal triangle with dim rows.
"""
T037027(dim::Int) = InvOrthoPoly(dim, n -> -1, n -> -1)

doc"""
Return the arctangent numbers (expansion of arctan``(x)^n/n!``).
"""
T049218(dim::Int) = InvOrthoPoly(dim, n -> 0, n -> n * (n + 1))

doc"""
Return the coefficients of Hermite polynomials ``H(n, (x-1)/√(2))/(√(2))^n``.
"""
T159834(dim::Int) = InvOrthoPoly(dim, n -> 1, n -> n)

doc"""
Return the coefficients of a variant of the Hermite polynomials.
"""
T137286(dim::Int) = InvOrthoPoly(dim, n -> 0, n -> n + 1)

doc"""
Return the coefficients of the Chebyshev-T polynomials.
"""
function T053120(len)
    T = SeqTriangle(len)
    R, x = PolynomialRing(ZZ, "x")
    m = 0
    for n in 0:len - 1
        f = chebyshev_t(n, x)
        for k in 0:n
            T[m] = coeff(f, k)
            m += 1
        end
    end
    T
end

doc"""
Return the coefficients of the Chebyshev-U polynomials.
"""
function T053117(len)
    T = SeqTriangle(len)
    R, x = PolynomialRing(ZZ, "x")
    m = 0
    for n in 0:len - 1
        f = chebyshev_u(n, x)
        for k in 0:n
            T[m] = coeff(f, k)
            m += 1
        end
    end
    T
end

end # module

module OrthoPolynomialsTest
using Base.Test, SeqBase, SeqTests, Nemo, OEISUtils, OrthoPolynomials

function test()
    @testset "OrthoPoly" begin

        @test isa(OrthoPoly(10, n -> 1, n -> n + 1)[end], fmpz)
        @test isa(InvOrthoPoly(10, n -> 1, n -> n + 1)[end], fmpz)

        @test RowSums(T217537(8)) == L217924(8)

        if oeis_isinstalled()

            T = [T066325, T049310, T137338, T104562, T037027, T049218, T159834,
                 T137286, T053120, T053117, T053121, T216916, T217537, T064189,
                 T202327, T111062, T099174, T111593, T059419, T026300]
            SeqTest(T, 'T')

            L = [L217924, L005425, L000085, L001464, L003723, L108624]
            SeqTest(L, 'L')

            # R = [L006229, L005773]; SeqTest(R, 'L')
        end
    end
end

function demo()
    T = T111593(8)
    Show(T, ',')
    println(RowSums(T))

    T = T217537(8)
    Show(T, ',')
    println(RowSums(T))

    T = T053117(8)
    Show(T)
    println(RowSums(T))
end

doc"""
T111062(1000) :: 1.452374 seconds (2.00 M allocations: 34.306 MB, 7.60% gc time)
T066325(1000) :: 1.939496 seconds (5.33 M allocations: 92.740 MB, 35.44% gc time)
T053120(1000) :: 0.249195 seconds (626.33 k allocations: 18.565 MB)

"""
function perf()
    T111062(2)
    T066325(2)
    T053120(2)
    gc()
    @time T111062(1000)
    gc()
    @time T066325(1000)
    gc()
    @time T053120(1000)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module BellNumbers
using Memoize, Nemo, SeqBase, Products, NumberTheory

export BellTrans, BellTriangle, BellNumberList, BellNumber
export A000110, L000110
export T137452, T264428, T137513, T104556, T001497, T132062, T039683, T203412
export T004747, T051141, T265606, T119274, T000369, T051142

# Note that the triangles in the OEIS often are missing the first column 1,0,0,0,...
# (the triangles are (1,1)-based). On the other hand the Bell transforme is
# (0,0)-based and always has as first column 1,0,0,0,...

# Other Stirling number related sequences are implemented in the module StirlingLah.

doc"""
Return a list of the first m Bell numbers (a.k.a. exponential numbers).
"""
function BellNumberList(m::Int)
    m == 0 && return fmpz[]
    R = SeqArray(m)
    R[0] = 1; m == 1 && return R
    R[1] = 1; m == 2 && return R

    A = SeqArray(m)
    A[0] = fmpz(1)
    for n in 1:m - 2
        A[n] = A[0]
        for k in n:-1:1
            A[k - 1] += A[k]
        end
        R[n + 1] = A[0]
    end
    R
end

# much slower:
#function BellNumberList(m::Int)
#    m == 0 && return fmpz[]
#    SeqArray(m, n -> Nemo.bell(n))
#end

doc"""
Return the n-th Bell number. Bell numbers count the ways to partition a set of
``n`` labeled elements.

julia> BellNumber(10)
115975
"""
BellNumber(n::Int) = Nemo.bell(n)

doc"""
Return the n-th Bell number ``B_n``.

julia> A000110(11)
678570
"""
A000110(n::Int) = Nemo.bell(n)

doc"""
Return a list of Bell numbers of length len.

julia> L000110(10)
[1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147]
"""
L000110(len::Int) = BellNumberList(len)

doc"""
The Bell transform transforms an integer sequence into an integer triangle;
also known as incomplete Bell polynomials.

Let ``X`` be an integer sequence, then
``$B_{n,k}(X) = \sum_{m=1}^{n-k+1} \binom{n-1}{m-1} X[m] B_{n-m,k-1}(X)$`` where
``B_{0,0} = 1, B_{n,0} = 0`` for ``n≥1, B_{0,k} = 0`` for ``k≥1``.
"""
BellTrans(n::Int, k::Int, X) = Belltrans(n, k, X)

@memoize function Belltrans(n::Int, k::Int, X::Array)
    a = fmpz(1); s = fmpz(0)

    if (n == 0) && (k == 0) return a end
    if (n == 0) || (k == 0) return s end

    for m in 1:n-k+1
        s += a * Belltrans(n - m, k - 1, X) * X[m]
        a = div(a * (n - m), m)
    end
    s
end

@memoize function Belltrans(n::Int, k::Int, X::Function)
    a = fmpz(1); s = fmpz(0)

    if (n == 0) && (k == 0) return a end
    if (n == 0) || (k == 0) return s end

    for m in 1:n-k+1
        s += a * Belltrans(n - m, k - 1, X) * X(m - 1)
        a = div(a * (n - m), m)
    end
    s
end

doc"""
The Bell triangle gathers the results of the Bell transform applied to the
initial segments of the input sequence.

Famously the sequence (1,1,1,...) is mapped to the triangle of the Stirling
set numbers.

julia> Show(BellTriangle(5, k -> 1))
1
0 1
0 1 1
0 1 3 1
0 1 7 6 1
"""
function BellTriangle(n::Int, seq)
    M = SeqTriangle(n)
    i = 0
    for j in 0:n - 1, k in 0:j
        M[i] = BellTrans(j, k, seq)
        i += 1
    end
    M
end

doc"""
Return the coefficients of the first ``n`` Abel polynomials.

julia> Show(T137452(5))
1
0 1
0 -2 1
0 9 -6 1
0 -64 48 -12 1
"""
T137452(n::Int) = BellTriangle(n, k -> (-k - 1)^k)

doc"""
Return the Bell transform of the Bell numbers.

julia> Show(T264428(5))
1
0 1
0 1 1
0 2 3 1
0 5 11 6 1
"""
T264428(n::Int) = BellTriangle(n, BellNumber)

# Keep the comments below.
# Stirling number related sequences are now implemented in the module StirlingLah.
# doc"""
# Return the triangle of Stirling set numbers (a.k.a. Stirling numbers of 2nd kind).
# """
# StirlingSetTriangle(n::Int) = BellTriangle(n, k -> 1)
# T048993(n::Int) = StirlingSetTriangle(n)

doc"""
Return the triangle of the coefficients of the Mittag-Leffler polynomials.

julia> Show(T137513(5))
1
0 2
0 0 4
0 4 0 8
0 0 32 0 16
"""
T137513(n::Int) = BellTriangle(n, k -> IsOdd(k) ? 0 : 2fac(k))

doc"""
Return the matrix inverse of coefficients of Bessel polynomials; essentially the
same as coefficients of modified Hermite polynomials T096713.

julia> Show(T104556(5))
1
0 1
0 -1 1
0 0 -3 1
0 0 3 -6 1
"""
T104556(n::Int) = BellTriangle(n, k -> k < 2 ? (-1)^k : 0)

doc"""
Return a triangle of coefficients of Bessel polynomials (better use A132062).
"""
T001497(n::Int) = BellTriangle(n, MultiFactorial(2, 1))

# Stirling number related sequences are now implemented in the module StirlingLah.
# Keep the comments below.
# doc"""
# Return the triangle of unsigned Stirling numbers of the first kind (a.k.a.
# Stirling cycle numbers).
# """
# StirlingCycleTriangle(n::Int) = BellTriangle(n, MultiFactorial(1,1))
# T132393(n::Int) = StirlingCycleTriangle(n)

doc"""
Return the triangle of coefficients of Bessel polynomials, also the Sheffer
triangle ``(1, 1 - √(1 - 2x))`` (Cf. A001497).

julia> Show(T132062(5))
1
0 1
0 1 1
0 3 3 1
0 15 15 6 1
"""
T132062(n::Int) = BellTriangle(n, MultiFactorial(2, 1))

doc"""
Return the signed double Pochhammer triangle: expansion of ``x(x-2)(x-4)..(x-2n+2)``.

julia> Show(T039683(5))
1
0 1
0 2 1
0 8 6 1
0 48 44 12 1
"""
T039683(n::Int) = BellTriangle(n, MultiFactorial(2, 2))

doc"""
Return the Bell transform of the MultiFactorial numbers of type (3,1).

julia> Show(T203412(5))
1
0 1
0 1 1
0 4 3 1
0 28 19 6 1
"""
T203412(n::Int) = BellTriangle(n, MultiFactorial(3, 1))

doc"""
Return the Bell transform of the MultiFactorial numbers of type (3,2).

julia> Show(T004747(5))
1
0 1
0 2 1
0 10 6 1
0 80 52 12 1
"""
T004747(n::Int) = BellTriangle(n, MultiFactorial(3, 2))

doc"""
Return the triangle ``3^{n-m}S1(n, m)`` where S1 are the signed Stirling numbers
of first kind.

julia> Show(T051141(5))
1
0 1
0 3 1
0 18 9 1
0 162 99 18 1
"""
T051141(n::Int) = BellTriangle(n, MultiFactorial(3, 3))

doc"""
Return the Bell transform of the quartic factorial numbers.

julia> Show(T265606(5))
1
0 1
0 1 1
0 5 3 1
0 45 23 6 1
"""
T265606(n::Int) = BellTriangle(n, MultiFactorial(4, 1))

doc"""
Return the triangle of coefficients of numerators in Pade approximation to ``e^x``.

julia> Show(T119274(5))
1
0 1
0 2 1
0 12 6 1
0 120 60 12 1
"""
T119274(n::Int) = BellTriangle(n, MultiFactorial(4, 2))

doc"""
Return the Bell transform of the MultiFactorial numbers of type (4,3).

julia> Show(T000369(5))
1
0 1
0 3 1
0 21 9 1
0 231 111 18 1
"""
T000369(n::Int) = BellTriangle(n, MultiFactorial(4, 3))

doc"""
Return the Bell transform of the MultiFactorial numbers of type (4,4).

julia> Show(T051142(5))
1
0 1
0 4 1
0 32 12 1
0 384 176 24 1 
"""
T051142(n::Int) = BellTriangle(n, MultiFactorial(4, 4))

end # module

module BellNumbersTest
using Base.Test, SeqBase, SeqTests, BellNumbers, OEISUtils, NumberTheory

function test()
    @testset "Bell" begin
        seq = [1, 1, 2, 5, 15, 52]
        len = size(seq)[1]
        a = [0, 52, 205, 210, 85, 15, 1]
        b = [BellTrans(len, k, seq) for k in 0:len]
        @test all(a .== b)

        seq = [1, 1, 1, 1, 1, 1]
        len = size(seq)[1]
        a = SeqArray([1, 0, 1, 0, 1, 1, 0, 1, 3, 1, 0, 1, 7, 6, 1, 0, 1, 15, 25,
                      10, 1, 0, 1, 31, 90, 65, 15, 1])
        b = BellTriangle(len, seq)
        @test all(a[0:len - 1] .== b[0:len - 1])

        a = [1, 1, 2, 5, 15, 52, 203]
        b = [BellNumber(n) for n in 0:6]
        @test all(a .== b)

        if oeis_isinstalled()
            SeqTest([A000110], 'A')
            SeqTest([L000110], 'L')

            T = [T264428, T137452, T132062, T265606]
            SeqTest(T, 'T')

            # These triangles have an additional column 1,0,0,... on the left
            # or have different signs. But they are essentially the same.
            t = [T104556, T001497, T039683, T203412, T004747, T051141,
                 T119274, T000369, T051142, T137513]
            SeqTest(t, 'P')
        end
    end
end

function demo()
    seq = [1, 1, 2, 5, 15, 52]
    len = size(seq)[1]

    for n in 0:len
        println([BellTrans(n, k, seq) for k in 0:n])
    end

    M = BellTriangle(len, seq)
    Show(M)

    for n in 0:6
        println(n, " ", BellNumberList(n))
    end

    T = T137452(8)
    Show(T)
end

doc"""
BellNumberList(1000) :: 0.539088 seconds (855.00 k allocations: 13.062 MiB)
BellTriangle(100, [1 for _ in 1:100]) :: 0.576952 seconds (1.22 M allocations: 24.571 MB, 14.23% gc time)
"""
function perf()
    gc()
    BellNumberList(5)
    @time BellNumberList(1000)
    @time BellTriangle(100, [1 for _ in 1:100])
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

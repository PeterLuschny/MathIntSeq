# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Fibonacci
using Nemo

export JacobsthalNumber, JacobsthalNumberList
export FibonacciNumber, FibonacciNumberList
export FibonacciGeneralized, FibonacciGeneralizedList
export A015445, A083856, A000045, A001045, A006130, A006131, A015440, A015441
export A015442, A015443, L000045, L001045, L006130, L006131, L015440, L015441
export L015442, L015443, L015445

doc"""
Return the generalized Fibonacci numbers which are defined recursively
``$F(n,k) = F(n-1,k) + kF(n-2,k), F(0,k) = 0, F(1,k) = 1.$``
"""
function FibonacciGeneralized(n::Int, k::Int)
    F = BigInt[1 k; 1 0]
    Fn = F^n
    fmpz(Fn[2, 1])
end

doc"""
Return the Jacobsthal numbers ``J(n) = F(n,2)`` where ``F(n,k)`` are the
generalized Fibonacci numbers.
"""
JacobsthalNumber(n::Int) = FibonacciGeneralized(n, 2)

doc"""
Return the Fibonacci numbers.
"""
FibonacciNumber(n::Int) = FibonacciGeneralized(n, 1)

doc"""
Return the Fibonacci numbers.
"""
A000045(n::Int) = FibonacciGeneralized(n, 1)

doc"""
Return the Jacobsthal numbers.
"""
A001045(n::Int) = FibonacciGeneralized(n, 2)

doc"""
Return the number of length-n strings with letters ``{0,1,2,3}`` where no two
consecutive letters are nonzero. (Joerg Arndt)
"""
A006130(n::Int) = FibonacciGeneralized(n, 3)

doc"""
Return the number of length-n strings with letters ``{0,1,2,3,4}`` where no two
consecutive letters are nonzero. (Joerg Arndt)
"""
A006131(n::Int) = FibonacciGeneralized(n, 4)

doc"""
Return the number of length-n strings with letters ``{0,1,2,3,4,5}`` where no
two consecutive letters are nonzero. (Joerg Arndt)
"""
A015440(n::Int) = FibonacciGeneralized(n, 5)

doc"""
Return the number of length-n strings with letters ``{0,1,2,...,6}`` where no
two consecutive letters are nonzero. (Joerg Arndt)
"""
A015441(n::Int) = FibonacciGeneralized(n, 6)

doc"""
Return the number of length-n strings with letters ``{0,1,2,...,7}`` where no
two consecutive letters are nonzero. (Joerg Arndt)
"""
A015442(n::Int) = FibonacciGeneralized(n, 7)

doc"""
Return the number of length-n strings with letters ``{0,1,2,...,8}`` where no
two consecutive letters are nonzero. (Joerg Arndt)
"""
A015443(n::Int) = FibonacciGeneralized(n, 8)

doc"""
Return the number of length-n strings with letters ``{0,1,2,...,9}`` where no
two consecutive letters are nonzero. (Joerg Arndt)
"""
A015445(n::Int) = FibonacciGeneralized(n, 9)

doc"""
Return the square array of generalized Fibonacci numbers, read by antidiagonals.
"""
A083856(r::Int, n::Int) = FibonacciGeneralized(n, r)

doc"""
Return the list of the first ``n`` terms of the generalized Fibonacci numbers
of order r.
"""
function FibonacciGeneralizedList(len::Int, r::Int)
    len ≤ 0 && return fmpz[]
    len < 2 && return [ZZ(i) for i in 0:len - 1]
    R = Array{fmpz}(len)
    R[1] = ZZ(0); R[2] = ZZ(1)
    for i in 3:len
        R[i] = R[i - 1] + r * R[i - 2]
    end
    R
end

doc"""
Return the list of the first ``n`` terms of the Fibonacci numbers.
"""
FibonacciNumberList(len::Int) = FibonacciGeneralizedList(len, 1)

doc"""
Return the list of the first ``n`` terms of the generalized Fibonacci numbers
of order 2.
"""
JacobsthalNumberList(len::Int) = FibonacciGeneralizedList(len, 2)

doc"""
Return the list of the first ``n`` terms of the Fibonacci numbers.
"""
L000045(len::Int) = FibonacciGeneralizedList(len, 1)

doc"""
Return the list of the first ``n`` terms of the Jacobsthal numbers.
"""
L001045(len::Int) = FibonacciGeneralizedList(len, 2)

doc"""
Return the list of the first ``n`` terms of the generalized Fibonacci numbers
of order 3.
"""
L006130(len::Int) = FibonacciGeneralizedList(len, 3)

doc"""
Return the list of the first ``n`` terms of the generalized Fibonacci numbers
of order 4.
"""
L006131(len::Int) = FibonacciGeneralizedList(len, 4)

doc"""
Return the list of the first ``n`` terms of the generalized Fibonacci numbers
of order 5.
"""
L015440(len::Int) = FibonacciGeneralizedList(len, 5)

doc"""
Return the list of the first ``n`` terms of the generalized Fibonacci numbers
of order 6.
"""
L015441(len::Int) = FibonacciGeneralizedList(len, 6)

doc"""
Return the list of the first ``n`` terms of the generalized Fibonacci numbers
of order 7.
"""
L015442(len::Int) = FibonacciGeneralizedList(len, 7)

doc"""
Return the list of the first ``n`` terms of the generalized Fibonacci numbers
of order 8.
"""
L015443(len::Int) = FibonacciGeneralizedList(len, 8)

doc"""
Return the list of the first ``n`` terms of the generalized Fibonacci numbers
of order 9.
"""
L015445(len::Int) = FibonacciGeneralizedList(len, 9)

end # module

module FibonacciTest
using Base.Test, SeqTests, SeqBase, Nemo, Fibonacci, OEISUtils

# References to the OEIS A-numbers are only approximative. Some sequences
# in the OEIS are missing the initial value '0' and starting at offset '1'.

# 0, 1, 1,  1,  1,   1,   1,    1,    1,     1, ... [A057427]
# 0, 1, 1,  2,  3,   5,   8,   13,   21,    34, ... [A000045]
# 0, 1, 1,  3,  5,  11,  21,   43,   85,   171, ... [A001045]
# 0, 1, 1,  4,  7,  19,  40,   97,  217,   508, ... [A006130]
# 0, 1, 1,  5,  9,  29,  65,  181,  441,  1165, ... [A006131]
# 0, 1, 1,  6, 11,  41,  96,  301,  781,  2286, ... [A015440]
# 0, 1, 1,  7, 13,  55, 133,  463, 1261,  4039, ... [A015441]
# 0, 1, 1,  8, 15,  71, 176,  673, 1905,  6616, ... [A015442]
# 0, 1, 1,  9, 17,  89, 225,  937, 2737, 10233, ... [A015443]
# 0, 1, 1, 10, 19, 109, 280, 1261, 3781, 15130, ... [A015445]

function test()
    @testset "Fibonacci" begin

        @test isa(FibonacciGeneralized(3, 30), fmpz)
        @test isa(FibonacciGeneralizedList(2, 20)[end], fmpz)

        @test FibonacciGeneralized(0, 1) == fmpz(0)
        @test FibonacciGeneralized(30, 2) == 357913941
        @test FibonacciGeneralized(2000, 1) == fmpz("4224696333392304878706725602341482782579852840250681098010280137314308584370130707224123599639141511088446087538909603607640194711643596029271983312598737326253555802606991585915229492453904998722256795316982874482472992263901833716778060607011615497886719879858311468870876264597369086722884023654422295243347964480139515349562972087652656069529806499841977448720155612802665404554171717881930324025204312082516817125")
        @test FibonacciGeneralized(9999, 1) == FibonacciGeneralizedList(10000, 1)[10000]

        if oeis_isinstalled()

            A = [ A015445, A000045, A001045, A006130, A006131, A015440,
                A015441, A015442, A015443 ]

            L = [ L000045, L001045, L006130, L006131, L015440, L015441, L015442,
                L015443, L015445 ]

            for C in [A, L], seq in C
                O = oeis_local(SeqName(seq), 10)
                # Per definition a Fibonacci sequence starts 0, 1, ...
                # This corrects some OEIS sequences:
                O[0] ≠ 0 ? s = 1 : s = 0
                C == A ? S = SeqArray(10, seq) : S = SeqArray(seq(10))
                #println(O); println(S)
                @test all(O[0:8] .== S[s:8+s])
            end
        end
    end
end

function demo()
    for r in 0:9
        println([FibonacciGeneralized(n, r) for n in 0:9])
    end

    for len in 0:6
        println(FibonacciGeneralizedList(len, 1))
    end
end

doc"""
for n in 1:10000 FibonacciGeneralized(n,1) end :: 1.139748 seconds (7.13 M allocations: 313.951 MB, 28.01% gc time)
FibonacciGeneralizedList(10000, 1) :: 0.011026 seconds (20.00 k allocations: 390.672 KB)
"""
function perf()
    gc()
    @time (for n in 1:10000 FibonacciGeneralized(n, 1) end)
    @time FibonacciGeneralizedList(10000, 1)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

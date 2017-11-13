# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module JacobiTheta
using Nemo, SeqBase

export JacobiTheta3Powers, JacobiTheta4Powers
export L000122, L002448, L004018, L104794, L005875, L213384, L000118, L035016
export L008452, L096727, L000132, L000141, L008451, L000143, L000144, L008453
export L000145, L276285, L276286, L276287, L004402, L004406, L004407, L015128
export L004403, L001934, L004404, L004405, L004408, L004409, L004410, L004411
export L004412, L004413, L004414, L004420, L004421, L004415, L004416, L004417
export L004418, L004419, L004422, L004423, L004424, L004425

doc"""
The ``q``-expansion to length len of the Jacobi theta function raised
to the power ``r``, i.e. ``ϑ(q)^r`` where ``ϑ(q) = 1 + ∑_{k ≥ 1} q^{k^2}``.
Number of ways of writing ``n`` as a sum of ``r`` squares.
"""
function JacobiTheta3Powers(len::Int, r::Int)
    len ≤ 0 && return fmpz[]
    R, x = PolynomialRing(ZZ, "x")
    e = theta_qexp(r, len, x)
    SeqArray([fmpz(coeff(e, j)) for j in 0:len - 1])
end

doc"""
Return the ``q``-expansion to length ``len`` of the Jacobi theta function raised
to the power ``r``, i.e. ``ϑ(-q)^r`` where ``ϑ(q) = 1 + ∑_{k ≥ 1} q^{k^2}``.
"""
function JacobiTheta4Powers(len::Int, r::Int)
    len ≤ 0 && return fmpz[]
    R, x = PolynomialRing(ZZ, "x")
    e = theta_qexp(r, len, -x)
    SeqArray([fmpz(coeff(e, j)) for j in 0:len - 1])
end

doc"""
Return the number of ways of writing a nonnegative integer n as a square.
"""
L000122(len::Int) = JacobiTheta3Powers(len, 1)

doc"""
Return the expansion of Jacobi theta function ``ϑ(-q)``.
"""
L002448(len::Int) = JacobiTheta4Powers(len, 1)

doc"""
Return the number of ways of writing a nonnegative integer n as a sum of 2 squares.
"""
L004018(len::Int) = JacobiTheta3Powers(len, 2)

doc"""
Return the expansion of ``ϑ_4(q)^2`` in powers of q.
"""
L104794(len::Int) = JacobiTheta4Powers(len, 2)

doc"""
Return the number of ways of writing a nonnegative integer n as a sum of 3 squares.
"""
L005875(len::Int) = JacobiTheta3Powers(len, 3)

doc"""
Return the expansion of ``ϑ_4(q)^3`` in powers of q.
"""
L213384(len::Int) = JacobiTheta4Powers(len, 3)

doc"""
Number of ways of writing a nonnegative integer n as a sum
of 4 squares.
"""
L000118(len::Int) = JacobiTheta3Powers(len, 4)

doc"""
Return the expansion of ``ϑ_4(q)^4`` in powers of q.
"""
L096727(len::Int) = JacobiTheta4Powers(len, 4)

doc"""
Return the number of ways of writing a nonnegative integer n as a sum
of 5 squares.
"""
L000132(len::Int) = JacobiTheta3Powers(len, 5)

doc"""
Return the number of ways of writing a nonnegative integer n as a sum
of 6 squares.
"""
L000141(len::Int) = JacobiTheta3Powers(len, 6)

doc"""
Return the nnumber of ways of writing a nonnegative integer n as a sum
of 7 squares.
"""
L008451(len::Int) = JacobiTheta3Powers(len, 7)

doc"""
Return the number of ways of writing a nonnegative integer n as a sum
of 8 squares.
"""
L000143(len::Int) = JacobiTheta3Powers(len, 8)

doc"""
Return the expansion of ``ϑ_4(q)^8`` in powers of q.
"""
L035016(len::Int) = JacobiTheta4Powers(len, 8)

doc"""
Return the number of ways of writing a nonnegative integer n as a sum
of 9 squares.
"""
L008452(len::Int) = JacobiTheta3Powers(len, 9)

doc"""
Return the number of ways of writing a nonnegative integer n as a sum
of 10 squares.
"""
L000144(len::Int) = JacobiTheta3Powers(len, 10)

doc"""
Return the number of ways of writing a nonnegative integer n as a sum
of 11 squares.
"""
L008453(len::Int) = JacobiTheta3Powers(len, 11)

doc"""
Return the number of ways of writing a nonnegative integer n as a sum
of 12 squares.
"""
L000145(len::Int) = JacobiTheta3Powers(len, 12)

doc"""
Return the number of ways of writing a nonnegative integer n as a sum
of 13 squares.
"""
L276285(len::Int) = JacobiTheta3Powers(len, 13)

doc"""
Return the number of ways of writing a nonnegative integer n as a sum
of 14 squares.
"""
L276286(len::Int) = JacobiTheta3Powers(len, 14)

doc"""
Return the nnumber of ways of writing a nonnegative integer n as a sum
of 15 squares.
"""
L276287(len::Int) = JacobiTheta3Powers(len, 15)

doc"""
Return the expansion of ``1/ϑ_3(q)`` in powers of q.
"""
L004402(len::Int) = JacobiTheta3Powers(len, -1)

doc"""
Return the expansion of ``1/ϑ_4(q)`` in powers of q.
"""
L015128(len::Int) = JacobiTheta4Powers(len, -1)

doc"""
Return the expansion of ``1/ϑ_3(q)^2`` in powers of q.
"""
L004403(len::Int) = JacobiTheta3Powers(len, -2)

doc"""
Return the expansion of ``1/ϑ_4(q)^2`` in powers of q.
"""
L001934(len::Int) = JacobiTheta4Powers(len, -2)

doc"""
Return the expansion of ``1/ϑ_3(q)^3`` in powers of q.
"""
L004404(len::Int) = JacobiTheta3Powers(len, -3)

doc"""
Return the expansion of ``1/ϑ_3(q)^4`` in powers of q.
"""
L004405(len::Int) = JacobiTheta3Powers(len, -4)

doc"""
Return the expansion of ``1/ϑ_3(q)^5`` in powers of q.
"""
L004406(len::Int) = JacobiTheta3Powers(len, -5)

doc"""
Return the expansion of ``1/ϑ_3(q)^6`` in powers of q.
"""
L004407(len::Int) = JacobiTheta3Powers(len, -6)

doc"""
Return the expansion of ``1/ϑ_3(q)^7`` in powers of q.
"""
L004408(len::Int) = JacobiTheta3Powers(len, -7)

doc"""
Return the expansion of ``1/ϑ_3(q)^8`` in powers of q.
"""
L004409(len::Int) = JacobiTheta3Powers(len, -8)

doc"""
Return the expansion of ``1/ϑ_3(q)^9`` in powers of q.
"""
L004410(len::Int) = JacobiTheta3Powers(len, -9)

doc"""
Return the expansion of ``1/ϑ_3(q)^{10}`` in powers of q.
"""
L004411(len::Int) = JacobiTheta3Powers(len, -10)

doc"""
Return the expansion of ``1/ϑ_3(q)^{11}`` in powers of q.
"""
L004412(len::Int) = JacobiTheta3Powers(len, -11)

doc"""
Return the expansion of ``1/ϑ_3(q)^{12}`` in powers of q.
"""
L004413(len::Int) = JacobiTheta3Powers(len, -12)

doc"""
Return the expansion of ``1/ϑ_3(q)^{13}`` in powers of q.
"""
L004414(len::Int) = JacobiTheta3Powers(len, -13)

doc"""
Return the expansion of ``1/ϑ_3(q)^{14}`` in powers of q.
"""
L004415(len::Int) = JacobiTheta3Powers(len, -14)

doc"""
Return the expansion of ``1/ϑ_3(q)^{15}`` in powers of q.
"""
L004416(len::Int) = JacobiTheta3Powers(len, -15)

doc"""
Return the expansion of ``1/ϑ_3(q)^{16}`` in powers of q.
"""
L004417(len::Int) = JacobiTheta3Powers(len, -16)

doc"""
Return the expansion of ``1/ϑ_3(q)^{17}`` in powers of q.
"""
L004418(len::Int) = JacobiTheta3Powers(len, -17)

doc"""
Return the expansion of ``1/ϑ_3(q)^{18}`` in powers of q.
"""
L004419(len::Int) = JacobiTheta3Powers(len, -18)

doc"""
Return the expansion of ``1/ϑ_3(q)^{19}`` in powers of q.
"""
L004420(len::Int) = JacobiTheta3Powers(len, -19)

doc"""
Return the expansion of ``1/ϑ_3(q)^{20}`` in powers of q.
"""
L004421(len::Int) = JacobiTheta3Powers(len, -20)

doc"""
Return the expansion of ``1/ϑ_3(q)^{21}`` in powers of q.
"""
L004422(len::Int) = JacobiTheta3Powers(len, -21)

doc"""
Return the expansion of ``1/ϑ_3(q)^{22}`` in powers of q.
"""
L004423(len::Int) = JacobiTheta3Powers(len, -22)

doc"""
Return the expansion of ``1/ϑ_3(q)^{23}`` in powers of q.
"""
L004424(len::Int) = JacobiTheta3Powers(len, -23)

doc"""
Return the expansion of ``1/ϑ_3(q)^{24}`` in powers of q.
"""
L004425(len::Int) = JacobiTheta3Powers(len, -24)

end # module

module JacobiThetaTest
using Base.Test, SeqTests, SeqBase, Nemo, JacobiTheta, OEISUtils

function test()
    @testset "JacobiTheta" begin

        @test JacobiTheta3Powers(0, 1) == fmpz[]
        @test isa(JacobiTheta3Powers(30, 1)[end], fmpz)
        @test isa(JacobiTheta4Powers(30, 1)[end], fmpz)

        @test L035016(999 + 1)[end] == -16565884160
        @test L035016(1000 + 1)[end] == 18365675328

        if oeis_isinstalled()

            A = [L000122, L002448, L004018, L104794, L005875, L213384, L000118,
                L035016, L008452, L096727, L000132, L000141, L008451, L000143,
                L000144, L008453, L000145, L276285, L276286, L276287, L004402,
                L004406, L004407, L015128, L004403, L001934, L004404, L004405,
                L004408, L004409, L004410, L004411, L004412, L004413, L004414,
                L004420, L004421, L004415, L004416, L004417, L004418, L004419,
                L004422, L004423, L004424, L004425 ]

            SeqTest(A, 'L')
        end
    end
end

function demo()
    for n in 1:6 println(n, ": ", JacobiTheta3Powers(8, n)) end
    for n in 0:6 println(-n, ": ", JacobiTheta3Powers(8, -n)) end
    println(L000143(10))
    println(L035016(10))
end

doc"""
L000143(100000) :: 0.077559 seconds (199.50 k allocations: 3.807 MB)
L035016(100000) :: 0.078684 seconds (199.51 k allocations: 3.807 MB)
"""
function perf()
    @time L000143(100000)
    @time L035016(100000)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

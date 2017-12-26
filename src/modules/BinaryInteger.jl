# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module BinaryInteger

export A001855, A003314, A033156, A054248, A061168, A083652, A097383, A123753
export A295513, BinaryIntegerLength

doc"""
Return the length of the binary extension of an integer ``n``, which is defined
as ``0`` if ``n = 0`` and for ``n > 0`` as ``⌊ log2(n) ⌋ + 1``.
"""
BinaryIntegerLength(n) = n == 0 ? 0 : Int(floor(log2(n))) + 1

doc"""
Alias for the function BinaryIntegerLength.
"""
bil(n) = BinaryIntegerLength(n)

doc"""
Return ``n ``BinaryIntegerLength``(n) - 2^``BinaryIntegerLength``(n)``.
"""
A295513(n) = n*bil(n) - 2^bil(n)

doc"""
Maximal number of comparisons for sorting ``n`` elements by binary insertion.
"""
A001855(n) = A295513(n) + 1

doc"""
Return the sum of lengths of binary expansions of ``0`` through ``n``.
"""
A083652(n) = A295513(n+1) + 2

doc"""
Recurrence ``a(n) = a(n-1) + ⌊ a(n-1)/(n-1) ⌋ + 2`` for ``m ≥ 2`` and ``a(1) = 1``.
"""
A033156(n) = A295513(n) + 2n

doc"""
Binary entropy function: ``a(n) = n + min { a(k) + a(n-k) : 1 ≤ k ≤ n-1 }`` for
``n > 1,`` and ``a(1) = 0``.
"""
A003314(n) = A295513(n) + n

doc"""
Binary entropy: ``a(n) = n + min { a(k) + a(n-k) : 1 ≤ k ≤ n-1 }.``
"""
A054248(n) = A295513(n) + n + rem(n, 2)

doc"""
Minimum total number of comparisons to find each of the values ``1`` through ``n``
using a binary search with ``3``-way comparisons.
"""
A097383(n) = A295513(n+1) - div(n-1, 2)

doc"""
Partial sums of the sequence ``⌊ log_2(n) ⌋``.
"""
A061168(n) = A295513(n+1) - n + 1

doc"""
Partial sums of the sequence of length of the binary expansion of ``2n+1``.
"""
A123753(n) = A295513(n+1) + n + 2

end # module

module BinaryIntegerTest
using BinaryInteger, Base.Test, SeqTests, SeqBase, OEISUtils

function test()
    @testset "BinaryInteger" begin

        @test A295513(0) == -1
        @test A295513(1) == -1
        @test A295513(2) == 0
        @test A295513(3) == 2

        if oeis_isinstalled()
            # A295513
            #SeqTest([A001855, A003314, A033156, A054248, A061168, A083652,
            #         A097383, A123753], 'A')
        end
    end
end

function demo()
    println([A295513(n) for n in 0:12])
    println([A123753(n) for n in 0:12])
    println([A001855(n) for n in 1:12])
    println([A083652(n) for n in 1:12])
    println([A033156(n) for n in 1:12])
    println([A003314(n) for n in 1:12])
    println([A054248(n) for n in 1:12])
    println([A097383(n) for n in 1:12])
    println([A061168(n) for n in 1:12])
end

doc"""
"""
function perf()
    @time [A295513(k) for k in 0:100000]
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

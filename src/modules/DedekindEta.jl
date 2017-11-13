# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module DedekindEta
using Nemo, SeqBase

export DedekindEtaPowers, RamanujanTau, RamanujanTauList, PartitionNumberList
export L010815, L002107, L010816, L000727, L000728, L000729, L000730, L000731
export L010817, L010819, L000735, L010820, L010821, L010822, L000739, L010823
export L010824, L010825, L010826, L010827, L010828, L010829, L000594, L010830
export L010831, L010832, L010833, L010834, L010835, L010836, L010837, L010840
export L010838, L010839, L010841, L000041, L000712, L000716, L023003, L023004
export L023005, L023006, L023007, L023008, L023009, L023010, L005758, L023011
export L023012, L023013, L023014, L023015, L023016, L023017, L023018, L023019
export L023020, L023021, L006922, L082556, L082557, L082558, L082559

doc"""
Compute the ``q``-expansion to length len of the Dedekind ``η`` function (without
the leading factor ``q^{1/24}``) raised to the power ``r``, i.e.
```math
(q^{-1/24} η(q))^r = ∏_{k ≥ 1} (1 - q^k)^r.
```
In particular, ``r = -1`` returns the generating function of the Partition
function ``p(k)`` and ``r = 24`` gives the Ramanujan tau function ``τ(k)``.
"""
function DedekindEtaPowers(len::Int, r::Int)
    len ≤ 0 && return fmpz[]
    R, x = PolynomialRing(ZZ, "x")
    e = eta_qexp(r, len, x)
    SeqArray([coeff(e, j) for j in 0:len - 1])
end

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)``.
"""
L010815(len::Int) = DedekindEtaPowers(len, 1)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^2``.
"""
L002107(len::Int) = DedekindEtaPowers(len, 2)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^3``.
"""
L010816(len::Int) = DedekindEtaPowers(len, 3)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^4``.
"""
L000727(len::Int) = DedekindEtaPowers(len, 4)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^5``.
"""
L000728(len::Int) = DedekindEtaPowers(len, 5)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^6``.
"""
L000729(len::Int) = DedekindEtaPowers(len, 6)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^7``.
"""
L000730(len::Int) = DedekindEtaPowers(len, 7)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^8``.
"""
L000731(len::Int) = DedekindEtaPowers(len, 8)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^9``.
"""
L010817(len::Int) = DedekindEtaPowers(len, 9)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{10}``.
"""
L010818(len::Int) = DedekindEtaPowers(len, 10)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{11}``.
"""
L010819(len::Int) = DedekindEtaPowers(len, 11)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{12}``.
"""
L000735(len::Int) = DedekindEtaPowers(len, 12)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{13}``.
"""
L010820(len::Int) = DedekindEtaPowers(len, 13)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{14}``.
"""
L010821(len::Int) = DedekindEtaPowers(len, 14)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{15}``.
"""
L010822(len::Int) = DedekindEtaPowers(len, 15)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{16}``.
"""
L000739(len::Int) = DedekindEtaPowers(len, 16)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{17}``.
"""
L010823(len::Int) = DedekindEtaPowers(len, 17)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{18}``.
"""
L010824(len::Int) = DedekindEtaPowers(len, 18)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{19}``.
"""
L010825(len::Int) = DedekindEtaPowers(len, 19)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{20}``.
"""
L010826(len::Int) = DedekindEtaPowers(len, 20)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{21}``.
"""
L010827(len::Int) = DedekindEtaPowers(len, 21)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{22}``.
"""
L010828(len::Int) = DedekindEtaPowers(len, 22)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{23}``.
"""
L010829(len::Int) = DedekindEtaPowers(len, 23)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{24}``.
"""
L000594(len::Int) = DedekindEtaPowers(len, 24)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{25}``.
"""
L010830(len::Int) = DedekindEtaPowers(len, 25)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{26}``.
"""
L010831(len::Int) = DedekindEtaPowers(len, 26)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{27}``.
"""
L010832(len::Int) = DedekindEtaPowers(len, 27)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{28}``.
"""
L010833(len::Int) = DedekindEtaPowers(len, 28)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{29}``.
"""
L010834(len::Int) = DedekindEtaPowers(len, 29)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{30}``.
"""
L010835(len::Int) = DedekindEtaPowers(len, 30)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{31}``.
"""
L010836(len::Int) = DedekindEtaPowers(len, 31)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{32}``.
"""
L010837(len::Int) = DedekindEtaPowers(len, 32)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{40}``.
"""
L010840(len::Int) = DedekindEtaPowers(len, 40)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{44}``.
"""
L010838(len::Int) = DedekindEtaPowers(len, 44)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{48}``.
"""
L010839(len::Int) = DedekindEtaPowers(len, 48)

doc"""
Compute the expansion of ``∏_{m≥1} (1 - q^m)^{64}``.
"""
L010841(len::Int) = DedekindEtaPowers(len, 64)

doc"""
List of the first values of the Ramanujan tau function.
"""
RamanujanTauList(len::Int) = DedekindEtaPowers(len, 24)

doc"""
Return Ramanujan's tau(n).
"""
RamanujanTau(n::Int) = DedekindEtaPowers(n, 24)[end]

doc"""
Return the number of partitions of n.
"""
PartitionNumberList(len::Int) = DedekindEtaPowers(len, -1)

doc"""
Return the number of partitions of n.
"""
L000041(len::Int) = DedekindEtaPowers(len, -1)

doc"""
Return the number of partitions of n into parts of 2 kinds.
"""
L000712(len::Int) = DedekindEtaPowers(len, -2)

doc"""
Return the number of partitions of n into parts of 3 kinds.
"""
L000716(len::Int) = DedekindEtaPowers(len, -3)

doc"""
Return the number of partitions of n into parts of 4 kinds.
"""
L023003(len::Int) = DedekindEtaPowers(len, -4)

doc"""
Return the number of partitions of n into parts of 5 kinds.
"""
L023004(len::Int) = DedekindEtaPowers(len, -5)

doc"""
Return the number of partitions of n into parts of 6 kinds.
"""
L023005(len::Int) = DedekindEtaPowers(len, -6)

doc"""
Return the number of partitions of n into parts of 7 kinds.
"""
L023006(len::Int) = DedekindEtaPowers(len, -7)

doc"""
Return the number of partitions of n into parts of 8 kinds.
"""
L023007(len::Int) = DedekindEtaPowers(len, -8)

doc"""
Return the number of partitions of n into parts of 9 kinds.
"""
L023008(len::Int) = DedekindEtaPowers(len, -9)

doc"""
Return the number of partitions of n into parts of 10 kinds.
"""
L023009(len::Int) = DedekindEtaPowers(len, -10)

doc"""
Return the number of partitions of n into parts of 11 kinds.
"""
L023010(len::Int) = DedekindEtaPowers(len, -11)

doc"""
Return the number of partitions of n into parts of 12 kinds.
"""
L005758(len::Int) = DedekindEtaPowers(len, -12)

doc"""
Return the number of partitions of n into parts of 13 kinds.
"""
L023011(len::Int) = DedekindEtaPowers(len, -13)

doc"""
Return the number of partitions of n into parts of 14 kinds.
"""
L023012(len::Int) = DedekindEtaPowers(len, -14)

doc"""
Return the number of partitions of n into parts of 15 kinds.
"""
L023013(len::Int) = DedekindEtaPowers(len, -15)

doc"""
Return the number of partitions of n into parts of 16 kinds.
"""
L023014(len::Int) = DedekindEtaPowers(len, -16)

doc"""
Return the number of partitions of n into parts of 17 kinds.
"""
L023015(len::Int) = DedekindEtaPowers(len, -17)

doc"""
Return the number of partitions of n into parts of 18 kinds.
"""
L023016(len::Int) = DedekindEtaPowers(len, -18)

doc"""
Return the number of partitions of n into parts of 19 kinds.
"""
L023017(len::Int) = DedekindEtaPowers(len, -19)

doc"""
Return the number of partitions of n into parts of 20 kinds.
"""
L023018(len::Int) = DedekindEtaPowers(len, -20)

doc"""
Return the number of partitions of n into parts of 21 kinds.
"""
L023019(len::Int) = DedekindEtaPowers(len, -21)

doc"""
Return the number of partitions of n into parts of 22 kinds.
"""
L023020(len::Int) = DedekindEtaPowers(len, -22)

doc"""
Return the number of partitions of n into parts of 23 kinds.
"""
L023021(len::Int) = DedekindEtaPowers(len, -23)

doc"""
Return the number of partitions of n into parts of 24 kinds.
"""
L006922(len::Int) = DedekindEtaPowers(len, -24)

doc"""
Return the number of partitions of n into parts of 30 kinds.
"""
L082556(len::Int) = DedekindEtaPowers(len, -30)

doc"""
Return the number of partitions of n into parts of 32 kinds.
"""
L082557(len::Int) = DedekindEtaPowers(len, -32)

doc"""
Return the number of partitions of n into parts of 48 kinds.
"""
L082558(len::Int) = DedekindEtaPowers(len, -48)

doc"""
Return the number of partitions of n into parts of 64 kinds.
"""
L082559(len::Int) = DedekindEtaPowers(len, -64)

end # module

module DedekindEtaTest
using Base.Test, SeqBase, SeqTests, Nemo, DedekindEta, OEISUtils

function test()
    @testset "DedekindEta" begin

        @test DedekindEtaPowers(0, 1) == fmpz[]
        @test isa(DedekindEtaPowers(30, 1)[10], fmpz)
        @test isa(DedekindEtaPowers(30, -1)[10], fmpz)

        @test RamanujanTau(20) == -7109760
        P0 = SeqArray([1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101])
        P1 = PartitionNumberList(14)
        @test all(P0 .== P1)

        if oeis_isinstalled()

            L = [L010815, L002107, L010816, L000727, L000728, L000729, L000730,
                L000731, L010817, L010819, L000735, L010820, L010821, L010822,
                L000739, L010823, L010824, L010825, L010826, L010827, L010828,
                L010829, L000594, L010830, L010831, L010832, L010833, L010834,
                L010835, L010836, L010837, L010840, L010838, L010839, L010841,
                L000041, L000712, L000716, L023003, L023004, L023005, L023006,
                L023007, L023008, L023009, L023010, L005758, L023011, L023012,
                L023013, L023014, L023015, L023016, L023017, L023018, L023019,
                L023020, L023021, L006922, L082556, L082557, L082558, L082559]

            SeqTest(L, 'L')
        end
    end
end

function demo()
    for n in 0:6
        println(n, ": ", DedekindEtaPowers(8, n))
    end

    for n in 0:6
        println(-n, ": ", DedekindEtaPowers(8, -n))
    end
end

doc"""
PartitionNumberList(10000) :: 0.083250 seconds (19.50 k allocations: 383.156 KB)
L000731(10000) :: 0.003809 seconds (19.50 k allocations: 383.156 KB)
RamanujanTauList(10000) :: 0.015598 seconds (19.50 k allocations: 383.125 KB)
"""
function perf()
    @time PartitionNumberList(10000)
    @time L000731(10000)
    @time RamanujanTauList(10000)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

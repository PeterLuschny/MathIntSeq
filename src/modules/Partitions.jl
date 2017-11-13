# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Partitions

export Partition, A080577

doc"""
Generates the integer partitions of ``n`` in lexicographic order.
Ported from Wilf/Nijenhuis "Combinatorial Algorithms". Cf. A080577.
"""
function NEXPAR(N::Int)
    R = Dict{Int,Int}()
    M = Dict{Int,Int}()
    NLAST = 0
@label(L10)
    N == NLAST && @goto(L20)
    NLAST = N
@label(L30)
    S = N
    D = 0
@label(L50)
    D = D + 1
    R[D] = S
    M[D] = 1
@label(L40)
    MTC = M[D] ≠ N
    # comment out when benchmarking
    PRINTPAR(D, R, M)
    ! MTC && return
    @goto(L10)
@label(L20)
    ! MTC && @goto(L30)
    SUM = 1
    R[D] > 1 && @goto(L60)
    SUM = M[D] + 1
    D = D - 1
@label(L60)
    F = R[D] - 1
    M[D] == 1 && @goto(L70)
    M[D] = M[D] - 1
    D = D + 1
@label(L70)
    R[D] = F
    M[D] = 1 + div(SUM, F)
    S = SUM % F
    S == 0 && @goto(L40)
    @goto(L50)
end

doc"""
Prints the partitions given in the format used in function NEXPAR.
"""
function PRINTPAR(D::Int, R::Dict{Int,Int}, M::Dict{Int,Int})
    F = Int64[]
    for i in 1:D
        append!(F, [R[i] for j in 1:M[i]])
    end
    println(F)
end

doc"""
Generates the integer partitions of ``n`` in lexicographic order.
"""
Partition(n::Int) = NEXPAR(n)

doc"""
Generates the integer partitions of ``n`` in lexicographic order.
"""
A080577(n::Int) = NEXPAR(n)

end # module

module PartitionsTest
using Partitions # Combinatorics

function test()
end

function demo()
    for i in 1:6
        Partition(i)
        println()
    end

    #for i in 1:6
    #    for p in Combinatorics.partitions(i)
    #        println(p)
    #    end
    #end
end

doc"""
Partition(10)  ::  0.000022 seconds (8 allocations: 1.188 KB)
Partition(20)  ::  0.000645 seconds (8 allocations: 1.188 KB)
Partition(30)  ::  0.002075 seconds (8 allocations: 1.188 KB)
Partition(40)  ::  0.014640 seconds (8 allocations: 1.188 KB)
Partition(50)  ::  0.078197 seconds (8 allocations: 1.188 KB)
Partition(60)  ::  0.351380 seconds (8 allocations: 1.188 KB)
Partition(70)  ::  1.482960 seconds (14 allocations: 3.938 KB)
Partition(80)  ::  5.608707 seconds (14 allocations: 3.938 KB)
Partition(90)  :: 20.271781 seconds (14 allocations: 3.938 KB)
Partition(100) :: 67.797670 seconds (14 allocations: 3.938 KB)
"""
function perf()
    # -- first comment out PRINTPAR(D,R,M) in NEXPAR (line 25)

    #for i in 10:10:100
    #    print("i=", i, ": ")
    #    @time Partition(i)
    #end

    # -- For comparison (in particular note the allocations!)
    #for i in 10:10:100
    #    print(i," : ")
    #    @time (for p in Combinatorics.partitions(i) p end)
    #end
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

# 10 :   0.000010 seconds (84 allocations: 5.688 KB)
# 20 :   0.000126 seconds (1.25 k allocations: 97.469 KB)
# 30 :   0.001182 seconds (11.21 k allocations: 973.563 KB)
# 40 :   0.015650 seconds (74.68 k allocations: 6.934 MB, 39.64% gc time)
# 50 :   0.075141 seconds (408.45 k allocations: 40.978 MB, 18.19% gc time)
# 60 :   0.316078 seconds (1.93 M allocations: 207.615 MB, 15.50% gc time)
# 70 :   1.407447 seconds (8.18 M allocations: 933.866 MB, 13.81% gc time)
# 80 :   5.728053 seconds (31.59 M allocations: 3.728 GB, 12.84% gc time)
# 90 :  20.943281 seconds (113.27 M allocations: 14.075 GB, 11.88% gc time)
# 100 : 73.561937 seconds (381.14 M allocations: 49.697 GB, 11.42% gc time)

# https://www.math.upenn.edu/~wilf/website/CombinatorialAlgorithms.pdf

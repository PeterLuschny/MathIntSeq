# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module SubFactorial
using Nemo
export Derangement, Subfactorial, A000166, A000255

function Derangement(n)
    prec = 100
    while prec < 10000
        CC = AcbField(prec)
        c = exp(CC(-1)) * Nemo.gamma(CC(n+1), CC(-1))
        b, i = unique_integer(c)
        if b return i end
        prec *= 2
    end
    println("n = $n gives an InexactError!")
    # throw(InexactError()) error()
end

Subfactorial(n) = Derangement(n)
A000166(n) = Derangement(n)

# A000255(n) =  GAMMA(n+3,-1)*exp(-1)/(n+1)
# Do not use Subfactorial(n+2)/(n+1) !!
function A000255(n)
    CC = AcbField(500)
    c = exp(CC(-1)) * Nemo.gamma(CC(n+3), CC(-1)) / CC(n+1)
    b, i = unique_integer(c)
    b && return i
    println("n = $n gives an InexactError!")
    # throw(InexactError()) error()
end

# todo A105927, A105928

end # module

module SubFactorialTest
using Base.Test, SeqBase, SeqTests, Nemo, OEISUtils, SubFactorial

function test()
end

function demo()
    for n in 0:30 println(Subfactorial(n)) end
    println(Derangement(450))
    for n in 0:30 println(A000255(n)) end
end

doc"""
(for n in 0:450 Subfactorial(n) end) :: 0.323053 seconds (74.02 k allocations: 2.359 MiB)
"""
function perf()
    gc()
    @time (for n in 0:450 Subfactorial(n) end)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

#=
sf = [
1,
0,
1,
2,
9,
44,
265,
1854,
14833,
133496,
1334961,
14684570,
176214841,
2290792932,
32071101049,
481066515734,
7697064251745,
130850092279664,
2355301661033953,
44750731559645106,
895014631192902121,
18795307255050944540,
413496759611120779881,
9510425471055777937262,
228250211305338670494289,
5706255282633466762357224,
148362637348470135821287825,
4005791208408693667174771274,
112162153835443422680893595673,
3252702461227859257745914274516,
97581073836835777732377428235481 ]
=#

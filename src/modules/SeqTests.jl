# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module SeqTests

using Base.Test, OEISUtils, SeqBase
export SeqTest

const ShowTest = false

function SeqTest(seqarray, kind)
    if kind == 'A' return SeqATest(seqarray) end
    if kind == 'B' return SeqBTest(seqarray) end
    if kind == 'L' return SeqLTest(seqarray) end
    if kind == 'T' return SeqTTest(seqarray) end
    if kind == 'Q' return SeqQTest(seqarray) end
    if kind == 'P' return SeqPTest(seqarray) end
end

function SeqATest(seqarray)
    for seq in seqarray
        name = SeqName(seq)
        O = oeis_local(name, 10)
        S = SeqArray(10, seq)
        if ShowTest
            println("A --> ", name); println(O); println(S)
        end
        AssertSeqArray(S)
        @test all(S .== O)
    end
end

function SeqBTest(seqarray)
    for seq in seqarray
        name = SeqName(seq)
        # the parameter is not 'length' but 'search bound'.
        O = oeis_local(name, 12)
        S = seq(300)
        if ShowTest
            println("B --> ", name); println(O); println(S)
        end
        AssertSeqArray(S)
        all(S[0:11] .== O[0:11])
    end
end

function SeqLTest(seqarray)
    for seq in seqarray
        name = SeqName(seq)
        O = oeis_local(name, 12)
        S = seq(12)
        if ShowTest
            println("L --> ", name); println(O); println(S)
        end
        AssertSeqArray(S)
        @test all(S .== O)
    end
end

function xSeqLTest(seqarray)
    for seq in seqarray
        name = SeqName(seq)
        O = oeis_local(name, 6)
        S = seq(6)
        if ShowTest
            println("L --> ", name); println(O); println(S)
        end
        AssertSeqArray(S)
    end
end

function SeqTTest(seqarray)
    for seq in seqarray
        name = SeqName(seq)
        O = oeis_local(name, 21)
        S = seq(6)
        if ShowTest
            println("T --> ", name); println(O); println(S)
        end
        AssertSeqArray(S)
        @test all(S .== O)
    end
end

function SeqPTest(seqarray)
    for seq in seqarray
        name = SeqName(seq)
        O = oeis_local(name, 28)
        S = seq(7)
        if ShowTest
            println("P --> ", name); Show(O); Show(S)
        end
        AssertSeqArray(S)
    end
end

end # module

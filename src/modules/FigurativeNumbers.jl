module FigurativeNumbers

using SeqBase

export PolygonalNumber, PyramidalNumber
export A014107, A095794, A067998, A080956, A001477, A000217, A000290, A000326
export A000384, A000566, A000567, A001106, A001107
export A005564, A058373, A254749, A000292, A000330, A002411, A002412
export A002413, A002414, A007584, A007585

doc"""
Return the polygonal number with shape k.
"""
function PolygonalNumber(n, k)
    s = div(n^2 * (k - 2) - n * (k - 4), 2)
    k < 2 ? -s : s
end

doc"""
Return the pyramidal number with shape k.
"""
function PyramidalNumber(n, k)
    s = div(3 * n^2 + n^3 * (k - 2) - n * (k - 5), 6)
    k < 2 ? -s : s
end

doc"""
Return the polygonal numbers of shape -2.
"""
A014107(n) = PolygonalNumber(n, -2)
doc"""
Return the polygonal numbers of shape -1.
"""
A095794(n) = PolygonalNumber(n, -1)
doc"""
Return the polygonal numbers of shape 0.
"""
A067998(n) = PolygonalNumber(n, 0)
doc"""
Return the polygonal numbers of shape 1.
"""
A080956(n) = PolygonalNumber(n, 1)
doc"""
Return the polygonal numbers of shape 2 (these are the natural numbers).
"""
A001477(n) = PolygonalNumber(n, 2)
doc"""
Return the polygonal numbers of shape 3 (the triangular numbers).
"""
A000217(n) = PolygonalNumber(n, 3)
doc"""
Return the polygonal numbers of shape 4 (the squares).
"""
A000290(n) = PolygonalNumber(n, 4)
doc"""
Return the polygonal numbers of shape 5 (the pentagonal numbers).
"""
A000326(n) = PolygonalNumber(n, 5)
doc"""
Return the polygonal numbers of shape 6 (the hexagonal numbers).
"""
A000384(n) = PolygonalNumber(n, 6)
doc"""
Return the polygonal numbers of shape 7 (the heptagonal numbers).
"""
A000566(n) = PolygonalNumber(n, 7)
doc"""
Return the polygonal numbers of shape 8 (the octagonal numbers).
"""
A000567(n) = PolygonalNumber(n, 8)
doc"""
Return the polygonal numbers of shape 9 (the nonagonal numbers).
"""
A001106(n) = PolygonalNumber(n, 9)
doc"""
Return the polygonal numbers of shape 10 (decagonal numbers).
"""
A001107(n) = PolygonalNumber(n, 10)

doc"""
Return the pyramidal numbers of shape -1.
"""
A005564(n) = PyramidalNumber(n, -1)
doc"""
Return the pyramidal numbers of shape 0.
"""
A058373(n) = PyramidalNumber(n, 0)
doc"""
Return the pyramidal numbers of shape 1.
"""
A254749(n) = PyramidalNumber(n, 1)
#doc"""
#Return the pyramidal numbers of shape 2 (triangular numbers).
#"""
#A000217(n) = PyramidalNumber(n, 2)
doc"""
Return the pyramidal numbers of shape 3 (tetrahedral numbers).
"""
A000292(n) = PyramidalNumber(n, 3)
doc"""
Return the pyramidal numbers of shape 4 (square pyramidal numbers).
"""
A000330(n) = PyramidalNumber(n, 4)
doc"""
Return the pyramidal numbers of shape 5 (pentagonal pyramidal numbers).
"""
A002411(n) = PyramidalNumber(n, 5)
doc"""
Return the pyramidal numbers of shape 6 (hexagonal pyramidal numbers).
"""
A002412(n) = PyramidalNumber(n, 6)
doc"""
Return the pyramidal numbers of shape 7 (heptagonal pyramidal numbers).
"""
A002413(n) = PyramidalNumber(n, 7)
doc"""
Return the pyramidal numbers of shape 8 (octagonal pyramidal numbers).
"""
A002414(n) = PyramidalNumber(n, 8)
doc"""
Return the pyramidal numbers of shape 9 (enneagonal pyramidal numbers).
"""
A007584(n) = PyramidalNumber(n, 9)
doc"""
Return the pyramidal numbers of shape 10 (decagonal pyramidal numbers).
"""
A007585(n) = PyramidalNumber(n, 10)

end # module

module FigurativeNumbersTest

using FigurativeNumbers, Base.Test, SeqBase, SeqTests, OEISUtils

function test()

    @testset "Figurative" begin
        if oeis_isinstalled()

            A = [A014107, A067998, A001477, A000217, A000290, A000326,
                A000384, A000566, A000567, A001106, A001107, A000292,
                A000330, A002411, A002412, A002413, A007584, A007585]
                # A080956, A095794, A005564, A058373, A254749, A002414
            SeqTest(A, 'A')
        end
    end
end

function demo()
    for k in -2:10
        A = [PolygonalNumber(n, k) for n in 0:9]
        println(k, " ", A)
    end

    for k in -2:10
        A = [PyramidalNumber(n, k) for n in 0:9]
        println(k, " ", A)
    end
end

doc"""
"""
function perf()
end

function main()
    test()
    demo()
    perf()
end

main()

#-2 [0, -1,  2, 9, 20, 35,  54,  77, 104, 135]
#-1 [0, -1,  1, 6, 14, 25,  39,  56,  76,  99]
# 0 [0, -1,  0, 3,  8, 15,  24,  35,  48,  63]
# 1 [0, -1, -1, 0,  2,  5,   9,  14,  20,  27]
# --------------------------------------------
# 2 [0, 1,  2,  3,  4,  5,   6,   7,   8,   9]
# 3 [0, 1,  3,  6, 10, 15,  21,  28,  36,  45]
# 4 [0, 1,  4,  9, 16, 25,  36,  49,  64,  81]
# 5 [0, 1,  5, 12, 22, 35,  51,  70,  92, 117]
# 6 [0, 1,  6, 15, 28, 45,  66,  91, 120, 153]
# 7 [0, 1,  7, 18, 34, 55,  81, 112, 148, 189]
# 8 [0, 1,  8, 21, 40, 65,  96, 133, 176, 225]
# 9 [0, 1,  9, 24, 46, 75, 111, 154, 204, 261]
#10 [0, 1, 10, 27, 52, 85, 126, 175, 232, 297]

# ==============================================
# -2 [0,-1,  1, 10, 30,  65, 119, 196, 300,  435]
# -1 [0,-1,  0,  6, 20,  45,  84, 140, 216,  315] A005564
#  0 [0,-1, -1,  2, 10,  25,  49,  84, 132,  195] A058373
#  1 [0,-1, -2, -2,  0,   5,  14,  28,  48,   75] A254749
# -----------------------------------------------
#  2 [0, 1,  3,  6, 10,  15,  21,  28,  36,   45] A000217
#  3 [0, 1,  4, 10, 20,  35,  56,  84, 120,  165] A000292
#  4 [0, 1,  5, 14, 30,  55,  91, 140, 204,  285] A000330
#  5 [0, 1,  6, 18, 40,  75, 126, 196, 288,  405] A002411
#  6 [0, 1,  7, 22, 50,  95, 161, 252, 372,  525] A002412
#  7 [0, 1,  8, 26, 60, 115, 196, 308, 456,  645] A002413
#  8 [0, 1,  9, 30, 70, 135, 231, 364, 540,  765] A002414
#  9 [0, 1, 10, 34, 80, 155, 266, 420, 624,  885] A007584
# 10 [0, 1, 11, 38, 90, 175, 301, 476, 708, 1005] A007585

end # module

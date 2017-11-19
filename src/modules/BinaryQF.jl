# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module BinaryQF
using Nemo, NumberTheory, SeqBase

export B002476, B008784, B031363, B034017, B035251, B038872, B038873
export B042965, B057126, B057127, B068228, B084916, B089270, B141158, B242660
export B243655, B244779, B244780, B244819, B243168, B244291, B007522

# export binaryQF
# export B001481, B002144, B002313, B002479, B003136
# export B007645, B028951, B033200, B033203, B033207
# export B035247, B045373, B047486, B056874, B244713

# Ported from [BinaryQuadraticForms](http://oeis.org/wiki/User:Peter_Luschny/BinaryQuadraticForms)
# which is also the documentation of this function.
# See also https://oeis.org/wiki/Binary_Quadratic_Forms_and_OEIS
#
# Note that we list the *positive* numbers of a form  whereas in the OEIS often
# the '0' is included. The reason for our choice is explained in Andrew Sutherland's
# "Introduction to Arithmetic Geometry", MIT Open Course Ware, Fall 2013.

doc"""
Return integers that are represented by the binary quadratic form
``a x^2 + b xy + c y^2`` over Z. Parameter 'subset' is in {"positive", "primitively",
"prime"}. Use it only as a 'private' function.
"""
function binaryQF(a::Int, b::Int, c::Int, bound = 100::Int, subset = "positive", verbose = false)

    α, β, γ = a, b, c

    discriminant() = fmpz(β^2 - 4 * α * γ)
    isreduced() = (-α < β ≤ α < γ) || (fmpz(0) ≤ β ≤ α == γ)

    function roots(a::Int, b::Int, c::Int, n::Int, y::Int)
#################################
    throw(ErrorException("not yet implemented"))
################################
#    x = var('x')
#    eq = a * x * x + b * x * y + c * y * y
#    (eq - n).roots(multiplicities = false, ring = ZZ)
end

    function sqr_disc(M, primitively = false)

    d = discriminant()
    d == 0 && throw(ValueError("discriminant must not be zero"))

    a, b, c = α, β, γ
    (a == 0 && c == 0) && return [b * n for n in 1:div(M, abs(b))]

    D = isqrt(d)
    # a must be ≠ 0
    if a == 0
        a, c = c, 0
    end
    k = 2 * D; m = 4 * a * D
    u = b + D; v = b - D
    S = fmpz[]

    # Solvability in Z.
    for n in 1:M
        h = fmpz(4 * a * n)  # a <> 0 and n <> 0
        for t in Divisors(h) # returns fmpz
            g = fmpz(div(h, t))
            # if divides(k, g - t) && divides(m, g * u - t * v)
            if rem(g - t, k) == 0 && rem(g * u - t * v, m) == 0

                if primitively
                    y = div(g - t, k)
                    R = roots(a, b, c, n, y)
                    if IsPrimeTo(R[1], y)
                        push!(S, n)
                        break
                    end
                else
                    push!(S, n)
                    break
                end
            end
        end
    end
    sort([s for s in Set(S)])
end

    function imag_prime(M)
#################################
    throw(ErrorException("not yet implemented"))
################################
#    solve = pari('qfbsolve')
#    Q = pari('Qfb')(α, β, γ)
#    p = 1
#    R = []
#    while true
#        p = next_prime(p)
#        p > M && break
#        solve(Q, p) && push!(R,p)
#    end
#    R
end

    function imag_primitively(M)

    a, b, c = α, β, γ
    d = c - div(b * b, 4 * a)
    A = []

    for y in 0:isqrt(div(M, d))
        r = y * b / (2 * a)
        s = sqrt((M - d * y * y) / a)
        for x in Int(round(ceil(-s - r))):Int(round(floor(s - r)))
            IsPrimeTo(x, y) && push!(A, a * x^2 + b * x * y + c * y^2)
        end
    end
    sort([s for s in Set(A)])
end

    function imag_positive(M)
#################################
    throw(ErrorException("not yet implemented"))
################################
#   L = [2*ZZ(α), ZZ(β), ZZ(β), 2*ZZ(γ)]
#   G = Matrix(ZZ, 2, 2, L)
#   A = pari('qfrep')(G, M, 1)
#   [k+1 for k in 0:M-1 if A[k] > 0]
end

    function primitive_reps(a, h, b, M, S)

    if a ≤ M
        push!(S, a)
        if b ≤ M
            push!(S, b)
            if a ≤ (M - b) && h ≤ (M - a - b)
                a ≤ (M - a - h) && primitive_reps(a, h + 2 * a, a + b + h, M, S)
                b ≤ (M - b - h) && primitive_reps(a + b + h, h + 2 * b, b, M, S)
            end
        end
    end
end

    function positive_primitives(bound, primitively)

    a, b, c = α, β, γ
    S = fmpz[]

    while true
        new_val = a + b + c
        if new_val > 0
            primitive_reps(a, b + 2 * a, new_val, bound, S)
            b += 2 * c
            a = new_val
        else
            if new_val < 0
                b += 2 * a
                c = new_val
            end
        end
        if a == α && b == β && c == γ break end
    end

    if ! primitively
        X = fmpz[]
        for p in S
            q = t = 1
            while q ≤ bound
                push!(X, q)
                q = t * t * p
                t += 1
            end
        end
        S = X
    end

    sort([s for s in Set(S)])
end

    function reduce_real()

    d = discriminant()
    IsSquare(d) && throw(ValueError("form must not have square discriminant"))
    droot = isqrt(d)
    a, b, c = α, β, γ

    while a ≤ 0 || c ≥ 0 || b ≤ abs(a + c)

        cAbs = c
        if cAbs < 0  cAbs *= -1 end

        # cAbs = 0 will not happen for a non square form
        delta = div(b + droot, 2 * cAbs)
        if c < 0  delta *= -1 end
        aa = c
        bb = 2 * c * delta - b
        cc = c * delta * delta - b * delta + a
        a, b, c = aa, bb, cc
    end

    return [a, b, c]
end

    function reduce_imag()

    a, b, c = α, β, γ
    if a < 0
        a, b, c = -a, -b, -c
    end
    d = discriminant()

    while true
        A = (a == c && b < 0) || (c < a)
        B = (-a == b && a < c) || (a < abs(b))

        !(A || B) && break

        if A
            a, b, c = c, -b, a
        end

        if B
            b -= 2 * a * div(b, 2 * a)
            if abs(b) > a  b -= 2 * a end
            c = div(b * b - d, 4 * a)
        end
    end

    return [a, b, c]
end

# Return the unique reduced form equivalent to BQF(a,b,c)
    function reduced_form()

    isreduced() && return [α, β, γ]

    if discriminant() ≥ 0
        return reduce_real()
    else
        return reduce_imag()
    end
end

   ### --- eval part of function starts here --- ###

    prime = false || subset == "prime"
    primitively = false || subset == "primitively"

    d = discriminant()
    d == 0 && throw(ValueError("discriminant must not be 0"))

    a, b, c = α, β, γ
    if verbose
        println("Original form [", a, ", ", b, ", ", c, "] with discriminant ", d)
    end

    if IsSquare(d)

        verbose && println("Square discriminant!")
        if prime primitively = false end # for efficiency
        pp = sqr_disc(bound, primitively)
        if prime pp = [m for m in pp if IsPrime(m)] end

    else

        α, β, γ = reduced_form()
        verbose && println("Reduced form  [", α, ", ", β, ", ", γ, "]")

        if d < 0
            if prime
                pp = imag_prime(bound)
            else
                if primitively
                    pp = imag_primitively(bound)
                else
                    pp = imag_positive(bound)
                end
            end
        ### real case, indefinite form ###
        else # d > 0 and not square
            if prime
                primitively = true
            end # for efficiency
            pp = positive_primitives(bound, primitively)
            if prime pp = [m for m in pp if IsPrime(m)] end
        end
    end

    if verbose
        msg0 = prime ? " primes " : " positive integers "
        msg1 = primitively ? "primitively " : ""
        msg2 = "represented up to "
        println("There are ", length(pp), msg0, msg1, msg2, bound)
    end
    SeqArray(pp)
end # binary_QF

doc"""
Return positive numbers of the form ``n = x^2-3y^2`` of discriminant 12.
"""
B084916(bound::Int) = binaryQF(1, 0, -3, bound)

doc"""
Return positive numbers that are primitively represented by the indefinite quadratic
form ``x^2 - 3y^2`` of discriminant 12.
"""
B243655(bound::Int) = binaryQF(1, 0, -3, bound, "primitively")

doc"""
Return primes congruent to 1 (mod 12).
"""
B068228(bound::Int) = binaryQF(1, 0, -3, bound, "prime")

doc"""
Return positive numbers of the form ``x^2+xy-2y^2``.
"""
B242660(bound::Int) = binaryQF(1, 1, -2, bound)

doc"""
Return positive numbers primitively represented by the binary quadratic form (1,1,-2).
"""
B244713(bound::Int) = binaryQF(1, 1, -2, bound, "primitively")

doc"""
Return primes of the form ``6m + 1``.
"""
B002476(bound::Int) = binaryQF(1, 1, -2, bound, "prime")

doc"""
Return positive numbers of the form ``x^2 - 2y^2`` with integers ``x, y``
(discriminant is 8).
"""
B035251(bound::Int) = binaryQF(1, 0, -2, bound)

doc"""
Return mumbers n such that 2 is a square mod n.
"""
B057126(bound::Int) = binaryQF(1, 0, -2, bound, "primitively")

doc"""
Return primes p such that 2 is a square mod p; or, primes congruent to
``{1, 2, 7}`` mod ``8``.
"""
B038873(bound::Int) = binaryQF(1, 0, -2, bound, "prime")

#doc"""
#Return primes congruent to ±1 mod 8 apart from the first term.
#"""
#B001132(bound::Int) = binaryQF(1, 0, -2, bound, "prime")

doc"""
Return positive numbers of the form ``x^2+xy-y^2``; or, of the form ``5x^2-y^2``.
"""
B031363(bound::Int) = binaryQF(1, 1, -1, bound) # "positive"

doc"""
Return positive numbers represented by the integer binary quadratic form ``x^2+xy-y^2``
with ``x`` and ``y`` relatively prime.
"""
B089270(bound::Int) = binaryQF(1, 1, -1, bound, "primitively")

doc"""
Return primes represented by the integer binary quadratic form ``x^2+xy-y^2``.
"""
B141158(bound::Int) = binaryQF(1, 1, -1, bound, "prime")

doc"""
Return primes congruent to ``{0, 1, 4}`` mod 5. (Cf. also A141158)
"""
B038872(bound::Int) = binaryQF(1, 1, -1, bound, "prime")

doc"""
Return positive integers not congruent to 2 mod 4; regular numbers modulo 4.
"""
B042965(bound::Int) = binaryQF(1, 0, -1, bound) # "positive"

doc"""
Return numbers that are congruent to ``{0, 1, 3, 5, 7}`` mod 8.
Positive integers represented by the binary quadratic form ``x^2-y^2``
with ``x`` and ``y`` relatively prime.
"""
B047486(bound::Int) = binaryQF(1, 0, -1, bound, "primitively")

#doc"""
#Return odd primes; primes represented by the binary quadratic form ``x^2-y^2``.
#"""
#B065091(bound::Int) = binaryQF(1, 0, -1, bound, "prime")

doc"""
Return positive integers of the form ``x^2 + xy + y^2`` (Loeschian numbers).
"""
B003136(bound::Int) = binaryQF(1, 1, 1, bound) # "positive"

doc"""
Return positive integers that are primitively represented by ``x^2 + xy + y^2``.
"""
B034017(bound::Int) = binaryQF(1, 1, 1, bound, "primitively")

doc"""
Return primes that are represented by ``x^2 + xy + y^2`` (generalized cuban primes).
"""
B007645(bound::Int) = binaryQF(1, 1, 1, bound, "prime")

doc"""
Return positive integers that are the sum of 2 squares.
"""
B001481(bound::Int) = binaryQF(1, 0, 1, bound) # "positive"

doc"""
Return numbers ``n`` that are primitively represented by ``x^2 + y^2``.
Also numbers n such that ``√(-1)`` mod ``n`` exists.
"""
B008784(bound::Int) = binaryQF(1, 0, 1, bound, "primitively")

doc"""
Return primes of form ``x^2 + y^2``; or primes congruent to 1 or 2 modulo 4.
"""
B002313(bound::Int) = binaryQF(1, 0, 1, bound, "prime")

doc"""
Return Pythagorean primes: primes of form ``4n + 1``.
"""
B002144(bound::Int) = binaryQF(1, 0, 1, bound, "prime")

doc"""
Return positive integers of the form ``x^2+xy+2y^2`` with ``x`` and ``y`` integers.
See also A035248.
"""
B028951(bound::Int) = binaryQF(1, 1, 2, bound) # "positive"

doc"""
Return positive numbers primitively represented by the binary quadratic form (1, 1, 2).
"""
B244779(bound::Int) = binaryQF(1, 1, 2, bound, "primitively")

doc"""
Return primes represented by the binary quadratic form (1, 1, 2). Primes congruent
to ``{0, 1, 2, 4}`` mod 7.
"""
B045373(bound::Int) = binaryQF(1, 1, 2, bound, "prime")

doc"""
Return primes of form ``x^2+7*y^2``.
"""
B033207(bound::Int) = binaryQF(1, 1, 2, bound, "prime")

doc"""
Return integers of form ``x^2 + 2y^2``.
"""
B002479(bound::Int) = binaryQF(1, 0, 2, bound)

doc"""
Return positive integers primitively represented by ``x^2 + 2y^2``.
"""
B057127(bound::Int) = binaryQF(1, 0, 2, bound, "primitively")

doc"""
Return primes of form ``x^2+2*y^2``. Primes congruent to ``{1, 2, 3}`` mod 8.
"""
B033203(bound::Int) = binaryQF(1, 0, 2, bound, "prime")

doc"""
Return odd primes of form ``x^2+2*y^2``. Primes congruent to ``{1, 3}`` mod 8.
"""
B033200(bound::Int) = binaryQF(1, 0, 2, bound, "prime")

doc"""
Return positive numbers represented by the binary quadratic form (1,1,3).
(See also A028954.)
"""
B035247(bound::Int) = binaryQF(1, 1, 3, bound) # "positive"

doc"""
Return positive numbers primitively represented by the binary quadratic form (1,1,3).
"""
B244780(bound::Int) = binaryQF(1, 1, 3, bound, "primitively")

doc"""
Return primes of form ``x^2+xy+3y^2``, discriminant -11.
"""
B056874(bound::Int) = binaryQF(1, 1, 3, bound, "prime")

doc"""
Return positive numbers primitively represented by the binary quadratic form (1,0,3).
"""
B244819(bound::Int) = binaryQF(1, 0, 3, bound, "primitively")

doc"""
Return positive integers of the form ``x^2+6xy-3y^2``.
"""
B243168(bound::Int) = binaryQF(1, 6, -3, bound)

doc"""
Return positive numbers primitively represented by the binary quadratic form (1,6,-3).
"""
B244291(bound::Int) = binaryQF(1, 6, -3, bound, "primitively")

doc"""
Return primes of the form ``8n+7``, that is, primes congruent to -1 mod 8.
"""
B007522(bound::Int) = binaryQF(-1, 4, 4, bound, "prime")

end # module

module BinaryQFTest
using Base.Test, SeqTests, SeqBase, Nemo, BinaryQF, OEISUtils

function test()

    # Note that we list the *positive* numbers of a form  whereas in the OEIS often
    # the '0' is included.

    Data = Dict{Int, Array{fmpz}}(
    002476 => [7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103],
    008784 => [1, 2, 5, 10, 13, 17, 25, 26, 29, 34, 37, 41],
    031363 => [1, 4, 5, 9, 11, 16, 19, 20, 25, 29, 31, 36],
    034017 => [1, 3, 7, 13, 19, 21, 31, 37, 39, 43, 49, 57],
    035251 => [1, 2, 4, 7, 8, 9, 14, 16, 17, 18, 23, 25],
    038872 => [5, 11, 19, 29, 31, 41, 59, 61, 71, 79, 89, 101],
    038873 => [2, 7, 17, 23, 31, 41, 47, 71, 73, 79, 89, 97],
    042965 => [1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16],
    057126 => [1, 2, 7, 14, 17, 23, 31, 34, 41, 46, 47, 49],
    057127 => [1, 2, 3, 6, 9, 11, 17, 18, 19, 22, 27, 33],
    068228 => [13, 37, 61, 73, 97, 109, 157, 181, 193, 229, 241, 277],
    084916 => [1, 4, 6, 9, 13, 16, 22, 24, 25, 33, 36, 37],
    089270 => [1, 5, 11, 19, 29, 31, 41, 55, 59, 61, 71, 79],
    141158 => [5, 11, 19, 29, 31, 41, 59, 61, 71, 79, 89, 101],
    242660 => [1, 4, 7, 9, 10, 13, 16, 18, 19, 22, 25, 27],
    243655 => [1, 6, 13, 22, 33, 37, 46, 61, 69, 73, 78, 94],
    244779 => [1, 2, 4, 7, 8, 11, 14, 16, 22, 23, 28, 29],
    244780 => [1, 3, 5, 9, 11, 15, 23, 25, 27, 31, 33, 37],
    244819 => [1, 3, 4, 7, 12, 13, 19, 21, 28, 31, 37, 39],
    243168 => [1, 4, 9, 13, 16, 24, 25, 33, 36, 37, 49, 52],
    244291 => [1, 4, 13, 24, 33, 37, 52, 61, 69, 73, 88, 97],
    007522 => [7, 23, 31, 47, 71, 79, 103, 127, 151, 167, 191, 199]
    )

    B = [B002476, B008784, B031363, B034017, B035251, B038872,
        B038873, B042965, B057126, B057127, B068228, B084916, B089270,
        B141158, B242660, B243655, B244779, B244780, B244819, B243168,
        B244291, B007522]

    @testset "BinaryQF" begin

        for seq in B
            S = seq(300)
            AssertSeqArray(S)
            anum = SeqNum(seq)
            data = SeqArray(Data[anum])
            # println(anum); println(S); println(data)
            @test all(S[0:11] .== data[0:11])
        end
    end
end

function demo()
    # binaryQF(1, 1, 1, bound)
    # binaryQF(1, 1, 1, bound, "primitively")
    # binaryQF(1, 1, 1, bound, "prime")
    println(B034017(50))
end

doc"""
B034017(10000) :: 0.006860 seconds (62.91 k allocations: 1.374 MB)
"""
function perf()
    gc()
    @time B034017(10000)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

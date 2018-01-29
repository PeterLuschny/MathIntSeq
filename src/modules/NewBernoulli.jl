module NewBern

using DataStructures, Nemo
import Base.show

export Sequence, show, numerate, vectorize

struct Sequence
    fun::Function
    count::Int
    A::Array{BigInt}
    function Sequence(fun, count)
        A = fill(BigInt(0), count)
        new(fun, count, A)
    end
end

NewBern.start(::Sequence) = 0
NewBern.eltype(::Sequence) = Nemo.fmpq
NewBern.length(S::Sequence) = S.count

function NewBern.next(S::Sequence, n)
    (S.fun(S.A, n), n + 1)
end

function NewBern.done(S::Sequence, n)
    n < S.count && return false
    for i in 1:S.count S.A[i] = BigInt(0) end
    return true
end

function show(S::Sequence)
    for s in S print(s, ", ") end
    println("...")
end

function vectorize(S::Sequence)
    [s for s in S]
end

function numerate(S::Sequence, offset=0)
    n = offset
    for s in S
        println(n, " ↦ ", s)
        n += 1
    end
end

function BernoulliNumber(A::Array{BigInt}, n::Int)
    if n < 2 A[2] = 1; return fmpq(1, n + 1) end
    n % 2 == 1 && return fmpq(0, 1)

    h = 1 + div(n, 2)
    for k in h:-1:2 A[k] += A[k+1] end
    for k in 2:1:h A[k] += A[k-1] end
    println(A[h])

    denom = BigInt(1) << (n+1) - BigInt(2)
    B = fmpq(A[h], denom)
    return n % 4 == 0 ? -B : B
end

bnum = Sequence(BernoulliNumber, 20)
show(bnum)
numerate(bnum)
BNUM = vectorize(bnum)
println(BNUM)

gc()
#@time Sequence(BernoulliNumber, 1000)
#@time vectorize(Sequence(BernoulliNumber, 1000))

# 0.000086 seconds (8 allocations: 8.195 KiB)
# 0.258567 seconds (767.70 k allocations: 141.253 MiB, 21.44% gc time)

function GenocchiNumber(A::Array{BigInt}, n::Int)
    if n < 2 A[2] = 1; return fmpz(1) end
    n % 2 == 1 && return fmpz(0)

    h = 1 + div(n, 2)
    for k in h:-1:2 A[k] += A[k+1] end
    for k in 2:1:h A[k] += A[k-1] end
    println(A)
    return A[h]  # A[2] Genocchi Medians A005439
end

gnum = Sequence(GenocchiNumber, 20)
show(gnum)
numerate(gnum)
GNUM = vectorize(gnum)
println("that")
println(GNUM)

function HammersleyPolynomial(A::Array{BigInt}, n::Int)
    if n < 2 A[1] = 1; return fmpz(1) end
    #n % 2 == 1 && return fmpz(0)

    h = 1 + div(n, 2)
    for k in h:-1:2 A[k] += A[k+1] end
    for k in 2:1:h A[k] += A[k-1] end
    return A[h]
end

hnum = Sequence(HammersleyPolynomial, 20)
show(hnum)
numerate(hnum)
HNUM = vectorize(hnum)
println(HNUM)

# 1, 1, 2, 7, 41, 376, 5033, 92821, 2257166, 69981919, 2694447797, 126128146156,
# 7054258103921, 464584757637001, 35586641825705882,
function A006846list(len::Int)
    R = Array{BigInt}(len)
    A = fill(BigInt(0), len+1); A[1] = 1

    for n in 1:len
        for k in n:-1:2 A[k] += A[k+1] end
        for k in 2:1:n A[k] += A[k-1] end
        R[n] = A[n]
        #R[n] = A[2]
        println(A)
    end
    return R
end

println("this")
println(A006846list(12))

doc"""
Return a list of the first len Bernoulli numbers ``B_n``. Cf. A027641/A027642.

julia> BernoulliList(10)
[1, 1//2, 1//6, 0, -1//30, 0, 1//42, 0, -1//30, 0]
"""
function BernoulliList(len::Int)
    if len ≤ 0 return fmpq[] end
    R = Array{fmpq}(len)
    R[1] = fmpq(1, 1); len == 1 && return R
    R[2] = fmpq(1, 2); len == 2 && return R

    A = Dict{Int,fmpz}(0 => 1, -2 => 0, -1 => 1, 1 => 0)
    a = fmpz(12); b = fmpz(240)
    k = e = 1

    for i in 2:len - 1
        Am = 0; A[k + e] = 0; e = -e
        for j in 0:i
            Am += A[k]; A[k] = Am; k += e
        end
        if e > 0
            R[i + 1] = fmpq(0, 1)
        else
            d = i >> 1
            R[i + 1] = d % 2 == 0 ? fmpq(-A[-d], a) : fmpq(A[-d], a)
            a, b = b, b << 4 + b << 2 - a << 6
        end
    end
    R
end

gc()
#BernoulliList(10)
#@time BernoulliList(1000)

end # module

module SeqBase

using Nemo, OffsetArrays

export SeqRange, SeqSize, SeqArray, SeqShow, SeqPrint, SeqNum, SeqName
export List, FillArray, IsSeqArray, AssertSeqArray
export TriangularNumber, IsTriangular, AssertTriangular, SeqTriangle
export Show, ShowAsMatrix, Row, RowSums, RowReverse
export ι, ιι, Δto□, □toΔ, Start, Last, HilbertHotel
export FindUpTo, FindInInterval, Iterator, Indicators
export Nth, IterateUpTo, Count, IndexIn, CountList, Enumerator, Accumulate
export IndicatorsFind, First, Previous, Next

# Base *everything* on Iterators
# Two flavors: 'List of lenght n' and 'List up to'.
# Do we need 0/1 idicators or are boolean indicators 'isA' enough?
# https://docs.julialang.org/en/latest/manual/interfaces.html

# keep private
ZeroRange(n) = 0:n - 1

doc"""
Return the size of a SeqArray.
"""
SeqSize(A) = Base.length(linearindices(A))

doc"""
Return the range of a SeqArray.
"""
SeqRange(A) = ZeroRange(SeqSize(A))

doc"""
Return a SeqArray of size n. (Preset with zeros.)
"""
SeqArray(n::Int) = n ≤ 0 ? fmpz[] : OffsetArray(fill(fmpz(0), n), ZeroRange(n))

doc"""
Return a SeqArray of size n. (Preset with zeros.)
"""
function SeqArray(len::Int, offset::Int)
    offset == 0 && return SeqArray(len)
    (offset < 0 || offset >= len) && error("Argument error offset: $(offset).")
    A = OffsetArray(Array{fmpz}(len), ZeroRange(len))
    for i in offset:len - 1 A[i] = 0 end
    A
end

doc"""
Return a SeqArray of size n preset with m.
"""
FillArray(m, n::Int) = OffsetArray(fill(fmpz(m), n), ZeroRange(n))

# private
const SeqType = OffsetArrays.OffsetArray{Nemo.fmpz,1,Array{Nemo.fmpz,1}}

doc"""
Return true if A is a SeqArray, otherwise false.
"""
function IsSeqArray(A)
    # ::Bool
    typeof(A) == SeqType
end

doc"""
Throws an ArgumentError if A is not a SeqArray.
"""
function AssertSeqArray(A)
    if !IsSeqArray(A)
        throw(ArgumentError(string(A) * " is not a SeqArray."))
    end
end

doc"""
Convert a tuple or a 1-based array A into a SeqArray.
"""
function SeqArray(A)
    IsSeqArray(A) && return A
    n = length(A)
    R = OffsetArray(Array{fmpz}(n), ZeroRange(n))
    for k in SeqRange(R) R[k] = fmpz(A[k + 1]) end
    R
end

# doc"""
# Given a boolean predicate 'isA' the function returns integers ``n`` which are
# isA for ``0 ≤ n ≤`` bound.
#
# julia> FindUpTo(7, IsPrime)
# [2, 3, 5, 7]
# """
# function FindUpTo(bound, isA::Function)
#    bound < 0 && return fmpz[]
#    SeqArray(filter(isA, 0:bound))
# end

doc"""
Return a iterator listing the values satisfying the predicate isA for arguments
in ``0 ≤ n ≤ bound .``
"""
function IterateUpTo(bound, isA)
    (i for i in 0:bound if isA(i))
end

doc"""
Return a SeqArray listing the values satisfying the predicate isA for arguments
``0 ≤ x ≤ `` bound.

julia> FindUpTo(7, IsPrime)
[2, 3, 5, 7]
"""
function FindUpTo(bound, isA)
    bound < 0 && return fmpz[]
    SeqArray(filter(isA, 0:bound))
end

#function FindUpTo(bound, isA::Function)
#    s(n) = n-1
#    A = find(!iszero, Indicators(bound+1, isA))
#    SeqArray(s.(A))
#end

doc"""
Given a boolean predicate 'isA' the function returns integers ``n`` which are
isA for ``a < n ≤ b``. This supports convenient partitioning of intervals.

julia> FindInInterval(7, 13, IsPrime)
[11, 13]
FindInInterval(13, 23, IsPrime)
[17, 19, 23]
"""
function FindInInterval(a, b, isA::Function)
    b < 0 && return fmpz[]
    SeqArray(filter(isA, a + 1:b))
end

doc"""
Return a list of length len of integers ``≥ 0`` which are isA.

julia> List(7, IsPrime)
[2, 3, 5, 7, 11, 13, 17]
"""
function List(len, isA::Function)
    len ≤ 0 && return fmpz[]
    j, c = Int(0), Int(0)
    A = OffsetArray(Array{fmpz}(len), ZeroRange(len))
    while c < len
        if isA(j)
            A[c] = fmpz(j)
            c += 1
        end
        j += 1
    end
    A
end

# 'hasparameters' and 'returnsboolean' from Chris Rackauckas.
hasparameters(f) = 0 < maximum((length(m.sig.parameters) for m in methods(f))...)
returnsboolean(f) = first((methods(f)...)).specializations.func.rettype == Bool

doc"""
Return a SeqArray of size n filled by the function f. If f is a boolean indicator
then List(n, f) is returned otherswise A[i] = f(i).
"""
function SeqArray(n::Int, f)
    n ≤ 0 && return fmpz[]

    if returnsboolean(f)
        return List(n, f)
    end

    if hasparameters(f)
        A = OffsetArray(Array{fmpz}(n), ZeroRange(n))
        for i in SeqRange(A) A[i] = f(i) end
        return A
    end

    error("Function not supported by this constructor!")
end

doc"""
Return a SeqArray of size n filled by the next n values taken from the channel.
"""
function SeqArray(n::Int, c::Channel)
    n ≤ 0 && return fmpz[]
    A = OffsetArray(Array{fmpz}(n), ZeroRange(n))
    for i in SeqRange(A) A[i] = take!(c) end
    A
end

doc"""
Return the n-th triangular number.
"""
TriangularNumber(n) = div(n * (n + 1), 2)

doc"""
Is n a triangular number?
"""
IsTriangular(n) = n == TriangularNumber(isqrt(2n))

doc"""
Return the sqrt of 2n or throw an ArgumentError if n is not a triangular number.
"""
function AssertTriangular(n)
    dim = isqrt(2n)
    n ≠ TriangularNumber(dim) && throw(ArgumentError("This is not a triangular array!"))
    dim
end

# A 'triangle' in MathIntSeq is a SeqArray of length (n+1)(n+2)/2.
# Thus we use linear indexing when generating a triangle.
# To display a triangle the linear index is transformed to a two dimensional
# index which is (0,0)-based and displayed row-wise as a lower triangular array.
# In particular note that SeqTriangle(size) has rows and cols from 0 to size-1.
# T(0,0)                          row 0
# T(1,0) T(1,1)	                  row 1
# T(2,0) T(2,1)	T(2,2)            row 2
# T(3,0) T(3,1)	T(3,2) T(3,3)     row 3
# col 0  col 1  col 2  col 3

doc"""
Return a trianguler array with n rows which is (0,0)-based.
"""
SeqTriangle(n::Int) = SeqArray(TriangularNumber(n))

doc"""
Convert a tuple or a 1-based array A into a SeqTriangle provided its length
is triangular.
"""
function SeqTriangle(A)
    if IsSeqArray(A)
        n = SeqSize(A)
        AssertTriangular(n)
        return A
    end

    n = length(A)
    AssertTriangular(n)
    R = OffsetArray(Array{fmpz}(n), ZeroRange(n))
    for k in SeqRange(R) R[k] = fmpz(A[k + 1]) end
    R
end

doc"""
Return a SeqTriangle with r rows generated by a function f(n::Int, k::Int).
"""
function SeqTriangle(r::Int, f::Function)
    SeqTriangle([f(k, j) for k in 0:r - 1 for j in 0:k])
end
    #r == 0 && return fmpz[]
    #len = TriangularNumber(r)
    #m = 0
    #R = OffsetArray(Array{fmpz}(len), ZeroRange(len))
    #for n in 0:r-1
    #    for k in 0:n R[m] = fmpz(f(n, k)); m += 1 end
    #end
    #R

doc"""
Return a (1,1)-based quadratic matrix of dimension dim preset with 0.
"""
SeqMatrix(dim) = dim ≤ 0 ? (return fmpz[]) : fill(fmpz(0), dim, dim)

doc"""
Print the SeqArray in the format 'n ↦ A[n]' for n in the range first:last.
"""
SeqShow(A, first::Int, last::Int) = SeqShow(A, first:last)

doc"""
Print the SeqArray in the format 'n ↦ A[n]'.
"""
SeqShow(A) = SeqShow(A, SeqRange(A))
function SeqShow(A, R)
    for i in R
        if isassigned(A, i)
            println(i, " ↦ ", A[i])
        else
            println(i, " ↦ ", "undef")
        end
    end
end

function print_without_type(io, v::AbstractVector)
    print(io, "[")
    for (i, el) in enumerate(v)
        i > 1 && print(io, ", ")
        print(io, el)
    end
    println(io, "]")
end

doc"""
Print the array without typeinfo.
"""
Println(v::AbstractVector) = print_without_type(IOContext(STDOUT; :compact => true), v)

doc"""
Print the SeqArray with or without typeinfo.
"""
function SeqPrint(v::AbstractVector, typeinfo = false)
    typeinfo ? println(v) : Println(v)
end

doc"""
Return row n of a lower triangular matrix (0 ≤ n).
"""
function Row(T, n::Int)
    AssertSeqArray(T)
    s = SeqSize(T)
    (s == 0 || n < 0) && return
    AssertTriangular(s)

    t = TriangularNumber(n + 1)
    s < t && error("This row is not in the matrix.")
    SeqArray([T[k] for k in t - n - 1:t - 1])
end

doc"""
Display a lower triangular matrix.

julia> Show(T225478(6))

    1
    3      4
   21     40     16
  231    524    336    64
 3465   8784   7136  2304   256
65835 180756 170720 72320 14080 1024
"""
function Show(T, separator = ' ')
    AssertSeqArray(T)
    n = SeqSize(T)
    n == 0 && return
    AssertTriangular(n)

    i = k = 0
    while k < n - 1
        for j = 0:i
            print(T[k], separator)
            k += 1
        end
        i += 1
        println()
    end
end

doc"""
Display the row n of a lower triangular matrix (0 ≤ n).
"""
function Show(T, n::Int, separator = ' ')
    R = Row(T, n)
    for r in R print(r, separator) end
    println()
end

doc"""
Display element k of row n of a lower triangular matrix (0 ≤ k ≤ n).
"""
function Show(T, n::Int, k::Int)
    n < 0 && return
    AssertSeqArray(T)
    s = SeqSize(T)
    t = TriangularNumber(n)
    index = k + t
    (k > n || s < index) && error("This element is not in the matrix.")
    println(T[index])
end

doc"""
Reverse the rows of a triangular array in place.
"""
function RowReverse(T)
    AssertSeqArray(T)
    n = SeqSize(T)
    n == 0 && return T
    dim = AssertTriangular(n)

    lo = hi = 0
    step = 1
    while hi ≤ n
        l = lo; h = hi
        while l < h
            T[l], T[h] = T[h], T[l]
            l += 1; h -= 1
        end
        step += 1
        lo, hi = hi + 1, hi + step
    end
    T
end

doc"""
Convert a lower triangular array to a square matrix.
"""
function Δto□(T)
    AssertSeqArray(T)
    n = SeqSize(T)
    n == 0 && return T
    dim = AssertTriangular(n)

    A = fill(fmpz(0), dim, dim)
    k = 0
    for r in 1:dim
        for j = 1:r
            A[r, j] = T[k]
            k += 1
        end
    end
    A
end

doc"""
Convert a square matrix to a triangular array.
"""
function □toΔ(M)
    n = length(M)
    n == 0 && return M
    dim = AssertSquare(n)
    len = TriangularNumber(dim)

    T = SeqArray(len)
    k = 0
    for r in 1:dim
        for j = 1:r
            T[k] = fmpz(M[r, j])
            k += 1
        end
    end
    T
end

doc"""
Print the triangle T as a matrix.
"""
ShowAsMatrix(T) = println(Δto□(T))

doc"""
Return the row sums of a triangle, if "alternate=true" the alternating row sums.
"""
function RowSums(T, alternate = false)
    AssertSeqArray(T)
    n = SeqSize(T)
    n == 0 && return T
    dim = AssertTriangular(n)
    S = SeqArray(dim)
    lo = hi = step = 0

    while hi < n
        if alternate
            s = sum((-1)^k * T[lo + k] for k in 0:hi - lo)
        else
            s = sum(T[k] for k in lo:hi)
        end
        S[step] = s
        step += 1
        hi += 1
        lo, hi = hi, hi + step
    end
    S
end

doc"""
Return the name of a OEIS sequence given a similar named function as a string.
"""
function SeqName(fun)
    aname = string(fun)
    for X in ['L', 'T', 'G', 'B', 'Q']
        aname = replace(aname, X, 'A')
    end

    if !ismatch(r"^A[0-9]{6}$", aname)
        fullname = split(aname, ".")
        aname = String(fullname[2])
        if !ismatch(r"^A[0-9]{6}$", aname)
            warn("Not a valid A-name!")
            return
        end
    end
    aname
end

doc"""
Return the A-number of a OEIS sequence given a similar named function as an integer.
"""
function SeqNum(seq)
    aname = SeqName(seq)
    parse(Int, aname[2:end])
end

doc"""
Iverson brackets.
"""
ι(b) = b ? 1 : 0

doc"""
Inverse Iverson brackets.
"""
ιι(n) = n == 0 ? true : false
 # The Unix way: A successful command returns a 0, while an unsuccessful one
 # returns a non-zero value that usually can be interpreted as an error code.

doc"""
Returns an integer which is the highest index in `b` for the value `a`.
Whenever `a` is not a member of `b` it returns -1.

julia> L = List(10, IsPrime)
IndexIn(13, L)
5
"""
function IndexIn(a, b::AbstractArray)
    bdict = Dict(zip(b, 0:length(linearindices(b))))
    get(bdict, fmpz(a), -1)
end

doc"""
Return the first element of the SeqArray A if A is not empty, 0 otherwise.
"""
Start(A) = A == [] ? "undef" : A[0]

doc"""
Return the element at the end of the list A if A is not empty, 0 otherwise.
"""
Last(A) = A == [] ? "undef" : A[SeqSize(A) - 1]  # A[end]

doc"""
Trick described by David Hilbert in a 1924 lecture "Über das Unendliche".
"""
HilbertHotel(guest, hotel) = prepend!(hotel, guest)

doc"""
Return the sqrt of ``n`` or throw an ArgumentError if ``n`` is not a square.
"""
function AssertSquare(n)
    dim = isqrt(n)
    n ≠ dim * dim && throw(ArgumentError("This is not a square!"))
    dim
end

doc"""
Enumerate the SeqArray A with linear index.
"""
function Enumerator(A)
    enumerate(IndexLinear(), A)
end

doc"""
Return a iterator listing the values of f for x in ``0≤x≤n``.
"""
function Iterator(n, f)
    (f(i) for i in 0:n)
end

doc"""
Return a iterator of length n reading from a channel.
"""
function Iterator(n, c::Channel)
    (take!(c) for i in 0:n - 1)
end

doc"""
Return a iterator of length n which has value 1 if isA(i) is true and otherwise 0.
"""
function Indicators(n, isA)
    (ι(isA(i)) for i in 0:n - 1)
end

doc"""
Return a list of length len which gives the numbers of integers ≤ n which
are isA.

julia> CountList(8, IsPrime)
[0, 0, 1, 2, 2, 3, 3, 4]
"""
CountList(len, isA) = Accumulate(Indicators(len, isA))

# Consider two sequences A and invA. We say invA is the left inverse of A iff
# invA(A(n)) = n and A(n) is the least number m such that A(invA(m)) = A(n).
# A(invA(n)) = n if and only if isA(n) = true where 'isA' is the indicator
# function of the sequence A.
# We introduce 'Count' as the left inverse of 'Nth'. For example
# if Nth(n) = NthPrime(n) then Count(n) = PrimePi(n) (A000720).

doc"""
Return the numbers of integers ≤ n which are isA.

julia> Count(8, IsPrime)
4
"""
function Count(n, isA)
    # +1 because of offset 0
    count(!iszero, Indicators(n + 1, isA))
end

doc"""
Return the Nth integer which is isA. (For N ≤ 0 return 0.)

julia> Nth(7, IsPrime)
17
"""
function Nth(N, isA::Function)
    N ≤ 0 && return 0
    n, c = Int(0), Int(0)
    while c < N
        i = isA(n)
        i && (c += 1)
        n += 1
    end
    n - 1
end

doc"""
Return the cumulative sum of an SeqArray.
"""
function Accumulate(A)
    R = SeqArray(SeqSize(A))
    i, acu = 0, 0
    for a in A
        acu += a
        R[i] = acu
        i += 1
    end
    R
end

doc"""
Return the smallest list of indicators of isA with sum(A) = count.

julia> IndicatorsFind(7, IsPrime)
[0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1]
"""
function IndicatorsFind(count, isA::Function)
    count ≤ 0 && return []
    n, c = Int(0), Int(0)
    A = Int[]
    while c < count
        i = isA(n)
        i && (c += 1)
        push!(A, ι(i))
        n += 1
    end
    A
end

doc"""
Return the first integer ``n ≥ 0`` such that isA(n) = true.

julia> First(IsPrime)
2
"""
function First(isA::Function)
    n = 0
    while !isA(n)
        n += 1
    end
    n
end

First(A::Array{Int}) = A == [] ? nothing : first(A)

doc"""
Return largest ``0 < k < n`` such that isA(k) = true or nothing if no such
``k`` exists.

julia> Previous(7, IsPrime)
5
"""
function Previous(n, isA::Function)
    n == nothing && return First(isA)
    while true
        n -= 1
        isA(n) && break
        # This definition avoids throwing an exception. Alternatively:
        # n < 0 && throw(ArgumentError("Not defined for $n."))
        n < 0 && return First(isA)
    end
    n
end

doc"""
Return least ``k > n ≥ 0`` such that isA(k) = true. NOTE: It is assumed that
such a ``k`` exists! (If not, the function will run forever.)

julia> Next(7, IsPrime)
11
"""
function Next(n, isA::Function)
    ((n ≤ 0) || (n == nothing)) && return First(isA)
    while true
        n += 1
        isA(n) && break
    end
    n
end

end # module

# ------------------------------------------------------

module SeqBaseTest
using Nemo, Base.Test, SeqBase, OffsetArrays

function test()
    BinTri(n) = SeqTriangle([binomial(k, j) for k in 0:n - 1 for j in 0:k])

    X(n, k) = ((k > n || k < 0) && return 0;
        (n == 0 && k == 0) && return 1;
        4 * X(n - 1, k - 1) + (4 * n - 1) * X(n - 1, k); )

    T(n) = SeqTriangle([X(j, k) for j in 0:n for k in 0:j])

    @testset "SeqBase" begin
        t = T(6)
        s = RowReverse(RowReverse(t))
        @test all(t .== s)

        a = SeqArray([3465, 8784, 7136, 2304, 256])
        b = Row(t, 4)
        @test all(a .== b)

        # RowSums
        t = BinTri(8)
        a = SeqArray([1, 2, 4, 8, 16, 32, 64, 128])
        b = RowSums(t)
        @test all(a .== b)

        a = SeqArray([1, 0, 0, 0, 0, 0, 0, 0])
        b = RowSums(t, true)  # alternating sum
        @test all(a .== b)
    end
end

function demo()
    IsPrime(n) = Nemo.isprime(fmpz(n))

    R = Array{fmpz}(7)
    println(IsSeqArray(R))

    R = Array{Int}(7)
    println(IsSeqArray(R))

    seq = SeqArray(7)
    println(IsSeqArray(seq))
    SeqPrint(seq)
    SeqShow(seq)

    seq = SeqArray([99, 2, 3, 4])
    SeqPrint(seq)
    SeqShow(seq)

    seq = SeqArray(4, x -> x^2)
    SeqPrint(seq)
    SeqShow(seq)

    seq = SeqArray(4, n -> IsPrime(n))
    SeqPrint(seq)
    SeqShow(seq)

    G() = Channel(csize = 2) do c
        for n in 0:1 put!(c, n) end
        n = 2
        while true
        IsPrime(n) && put!(c, fmpz(n))
        n += 1
    end
    end

    seq = SeqArray(4, G())
    SeqPrint(seq)
    SeqShow(seq)

    X(n, k) = ((k > n || k < 0) && return 0;
        (n == 0 && k == 0) && return 1;
        4 * X(n - 1, k - 1) + (4 * n - 1) * X(n - 1, k); )

    dim = 6
    T = SeqTriangle(dim, X)

    # The four ways to display an SeqTriangle.
    println(T) # as an array
    println("--")
    Show(T)   # as a triangle
    println("--")
    for n in 0:dim - 1 Show(T, n) end  # by rows
    println("--")
    for n in 0:dim - 1, k in 0:n Show(T, n, k) end  # by elements
    println("--")
    ShowAsMatrix(T)
    println("--")

    M = [1 0 0 0 0 0; 3 4 0 0 0 0; 21 40 16 0 0 0;
         231 524 336 64 0 0; 3465 8784 7136 2304 256 0;
         65835 180756 170720 72320 14080 1024]
    T = □toΔ(M)
    Show(T)
    println("--")
    for (a, b) in Enumerator(T)
        println(a, " ", b)
    end
    println("--")
    offset = 3
    len = 7
    SB = SeqArray(len, offset)
    println(SB)
    println("--")
    SeqShow(SB)
    println("--")
end

function perf()
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

* TODO

function roots(a::Int, b::Int, c::Int, n::Int, y::Int)
    throw(ErrorException("not yet implemented"))

function imag_prime(M)
    throw(ErrorException("not yet implemented"))

function imag_positive(M)
    throw(ErrorException("not yet implemented"))

Precision handling for Hyper1F1?

Add the inverse Bell transform to the Bell module
https://oeis.org/wiki/User:Peter_Luschny/BellTransform#.E2.99.A6.C2.A0The_inverse_Bell_transform

* Improve Workarounds

Look at the regression after Divisors(n::Int) -> Divisors(n::fmpz)

function FibonacciGeneralized(n::Int, k::Int)
    F = BigInt[1 k; 1 0]  # does not work with fmpz!
    Fn = F^n
    fmpz(Fn[2, 1])
end

Make it triangular from the start to eliminate □toT.

function InvOrthoPoly(s::Function, t::Function, dim::Int)
    dim ≤ 0 && return fmpz[]
    T = zeros(ZZ, dim, dim)
    for n in 1:dim T[n,n] = 1 end
    for n in 1:dim-1
        for k in 1:n+1
            T[n+1,k] = (k>1?T[n,k-1]:0)-s(n-1)×T[n,k]-(n>1?t(n-1)×T[n-1,k]:0)
        end
    end
□toΔ(T) end

* Better interface with Nemo

    - denominator(::fmpq)
    - numerator(::fmpq)
    - Divisors(::fmpz)
    - PrimeDivisors(::fmpz)

and

    E = Array{fmpz}()  # works
    E = fmpz[]         # works
    E = ZZ[]           # fails

    E = fmpz(1)       # works
    E = ZZ(1)         # works
    E = ones(fmpz,1)  # fails
    E = ones(ZZ, 1)   # fails

    E = fmpz(0)       # works
    E = ZZ(0)         # works
    E = zeros(fmpz,1) # fails
    E = zeros(ZZ, 1)  # works !

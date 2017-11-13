# ↓
# ↑ # uparrow
# ⤓ # DownArrowBar
# ⤒ # UpArrowBar

n = 8
k = 3

krono(n, k) = ∏(n - k + 1, n)
⊕(n, k) = krono(n, k)
c = 8 ⊕ 3
println(c)

#FallingFactorial(n::Int, k::Int) = ∏(n - k + 1, n)
↓(n, k) = FallingFactorial(n, k)
c = 8 ↓ 3
println(c)

#RisingFactorial(n::Int, k::Int) = ∏(n, n + k - 1)
↑(n, k) = RisingFactorial(n, k)
c = 8 ↑ 3
println(c)

using Printf
using Statistics

a11 = [n * m for n in 1:5, m in 1:5]
println(a11)
a12 = rand(0 :9, 5, 5)
for n = 1:5, m= 1:5
    @printf("%d | ",a12[n, m])
end
println()

d1 = Dict("mango"=>"tartar")
for kv in d1 println(kv) end
println(in("mango"=>"tartar"))

a = 10
A = 20
@printf("A : %s\na: %s\n", A, a)

ar = [1 2 3; 4 5 6; 7 8 9]
ar = ar[1:2,2:3]
println(ar)
using Printf
using Statistics
using LinearAlgebra

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
@printf("A : %s\na[ 0.02 ]: %s\n", A, a)

ar = [1 2 3; 4 5 6; 7 8 9]
ar = ar[1:2,2:3]
println(ar)

# Getting a diagonal matrix with a given size and given diagonal value:
println(Diagonal([0.02 for cat in 1:5]))
println("printing range")
println([1,2,3,4,5][1:2])
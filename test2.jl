
include("./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/BC2_eth_pred.jl")
include("./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/BC2.jl")
include("./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/ME.jl")
include("./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/time.jl")
include("./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/timeE.jl")
# Testing how to make deriavtives
# using Calculus 
using Symbolics, Printf

@variables x y

# j = differentiate("x^2+4y", :x)
# display(j)

f = x^2+4*y

# substitute(f, Dict(x=>2, y=>1))
f_num = eval(build_function(f,x,y))
some = f_num(2,1)

@printf("The result is %d\n", some)
Dx = Differential(x)
g = Dx(f)
# display(g)

myMatrix = [x^2 y; x^3 y^2]

# display(myMatrix)
# display(typeof(myMatrix))
# display(myMatrix[:])
# display(typeof(myMatrix[:]))

j = Symbolics.jacobian(myMatrix[:], [x,y])
# display(j)
@variables P[1:5, 1:5] # Creates a 5x5 symbolic matrix
# display(P)

@variables x, y ,z
x_num = 2
y_num = 3
f = x + y + z
f = substitute(f, Dict(x=>x_num, y=>y_num))
f_num = eval(build_function(f,z))
# display(f)

arr1 = [1 2 3]
mat1 = [2 3 4; 5 6 7]
matConc = [mat1[:];arr1[:]]
@printf("O vetor é:\n")
display(arr1)
@printf("\n")
display(arr1[:])
@printf("\n")
@printf("A matriz final é:\n")
display(mat1)
@printf("\n")
@printf("A matriz final é:\n")
display(mat1[:])
@printf("\n")
@variables P[1:5, 1:5]
display(P)
@printf("\n")
@printf("SomeMat is:\n")
someMat = [1,2,3]
display(length(someMat))
#@printf("Now for the matrices...\n")
#@variables matA[1:2, 1:2]
#@variables matB = [1:2, 1:2]
#matC = matA * matB
#@printf("matC is:\n")
#display(matC[1,1])
#@printf("\n")

somearr1 = [1 2 3 4 5;1 2 3 4 5;1 2 3 4 5;1 2 3 4 5;1 2 3 4 5]
println("The size of somearr1 is:")
println(size(somearr1))
println("The size of somearr1[:] is:")
println(size(somearr1[:]))

@variables somearr2[1:5,1:5]
somearr2 = 5*somearr2
# somearr = [1 2 3 4 5;1 2 3 4 5;1 2 3 4 5;1 2 3 4 5;1 2 3 4 5]
println("The size of somearr2 is:")
println(size(somearr2))
println("The size of somearr2[:] is:")
println(size(somearr2[:]))
println("somearr2 is:")
println(somearr2)
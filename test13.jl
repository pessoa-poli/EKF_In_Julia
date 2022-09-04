# Finding how to efficiently zip variables from a function and and elements from a vetor
# This is needed to build the Dict used in the substitute function from the symbolics library

using Symbolics

@variables x y z w
mat =[x + y + z;
      2x]
    
function mySolver(du,u,p,t)    
    du = [eval(build_function(cat,Symbolics.get_variables(cat)))()]
end

fun_callable = eval(build_function(f,array_of_symbolic_function_variables))
res = fun_callable(values_to_replace_symbols)
attributes_dictionary = Dict(x=>1, y=>1, z=>1,w=>1)
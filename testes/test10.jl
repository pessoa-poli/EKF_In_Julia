# Finding how to efficiently zip variables from a function and and elements from a vetor
# This is needed to build the Dict used in the substitute function from the symbolics library

using Symbolics

@variables x y z w
f = x + y + z
array_of_symbolic_function_variables = Symbolics.get_variables(f)
values_to_replace_symbols = [1 2 3]
array_of_2_tuples = zip(array_of_symbolic_function_variables,values_to_replace_symbols)
println(array_of_2_tuples)

fun_callable = eval(build_function(f,array_of_symbolic_function_variables))
res = fun_callable(values_to_replace_symbols)
attributes_dictionary = Dict(x=>1, y=>1, z=>1,w=>1)
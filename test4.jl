using Symbolics
using Latexify

@variables a b c
f =  a + b*c
display(latexify(f))

f1 = substitute(f, Dict([a=>1.2, b=>2.2]))
final = substitute(f1, Dict([c=>3]))
println(final)
println("Typeof final is:")
println(typeof(final.val))
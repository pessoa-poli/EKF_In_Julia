using symbolics
arr = [1 2 3 ;4 5 6]
lines, columns = size(arr)
println(lines,'-', columns)
println(arr[1:columns])

@variables i j
arr[i,j] = symbols("x$i$j")
display(arr)
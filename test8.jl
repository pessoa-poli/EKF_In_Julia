arr1 = [1 2 3]
arr2 = [cat*2 for cat=arr1]
println(arr2)
arr3 = [1:5,1:5]
arrRows=size(arr3)
arr4=size([1:5 1:5])[1]
println("Size $arr4")
# println(arrRows)
# println([cat for cat=[1:size([1:5 1:5])[1]]])
using OrbitAl.permutation

a = Perm([3, 8, 7, 2, 1, 4, 6, 9, 5])
a1 = inv(a)
b = Perm([3, 2, 4, 8, 7, 9, 5, 6, 1])
c = Perm([4, 6, 8, 3, 2, 5, 1, 7, 9])

@assert (a * b) * c == a * (b * c)
@assert (collect(1:9)^a)^b == collect(1:9)^(a*b)


using OrbitAl.syt

list = [1, 3, 6, 10]

@assert newtonDif(newtonSum(list)) == list
@assert newtonSum(newtonDif(list)) == list

@assert newtonDifR(newtonSumR(list)) == list
@assert newtonSumR(newtonDifR(list)) == list

set, composition = [4, 5, 7], [1,1,1,3,2,1]
@assert compositionSubset(9, set) == composition
@assert subsetComposition(9, composition) == set

@assert length(partitions(6)) == 11

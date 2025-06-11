using Test
using OrbitAl

@testset "Permutation Basics" begin
    p = Perm([2, 3, 1])
    q = p^2
    r = p * q

    @test degree(p) == 3
    @test domain(p) == 1:3
    @test isidentity(p^3)
    @test inv(p) * p == one(p)
    @test r == p * q
    @test p / q == p * inv(q)
    @test sign(one(p)) == 1
end

@testset "Cycle Operations" begin
    p = Perm([2, 3, 1, 5, 4])
    @test shape(p) == [3, 2]
    @test order(p) == 6
    @test sign(p) == -1
end

@testset "Point and Vector Action" begin
    p = Perm([3, 1, 2])
    v = [10, 20, 30]
    @test 1^p == 3
    @test [1,2,3]^p == [3,1,2]
    @test permuted(v, p) == v[inv(p).list]
end

@testset "Random and Identity Checks" begin
    p = rand(Perm, 10)
    @test length(p.list) == 10
    @test isidentity(one(p))
end

@testset "Transposition and Moved Point" begin
    p = transposition(5, 2, 4)
    @test p.list[2] == 4
    @test p.list[4] == 2
    @test last_moved(p) == 4
end

using OrbitAl.bfsdfs

##  a tree
nodes = [Node(i, []) for i in 1:7]
for (i,k) in pairs([3,4,4,5,6,6,6])
  i == k || push!(nodes[k].next, nodes[i])
end
root = nodes[6]

##  a visitor
pr(x) = print(x.id, ", ")

DFS(root, pr)
println()
BFS(root, pr)
println()
tree_print(root)

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

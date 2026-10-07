using Test
using OrbitAl

using OrbitAl.permutation
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
    p = Perm(5, [[2, 4]])
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


using OrbitAl.orbits

@testset "Orbit Basics" begin
    aaa = transpositions(3)
    o = orbit(aaa, 1, ^)
    @test length(o) == 3
    o = orbit(aaa, aaa[1], *)
    @test length(o) == 6
end

using OrbitAl.presentations
using OrbitAl.enumerator

@testset "Coset Enumeration" begin
    for (name, index, order) in [(:G4, 8, 24), (:G12, 24, 48), (:G333, 9, 54)]
        G = getfield(presentations, name)
        T = coset_table(G, G.sbgp)
        @test T.active == index == length(active_cosets(T))
        @test length(perms(T)) == length(G.gens)
        @test sizeOfGroup(PermGp(coset_table(G, Vector{Int}[]))) == order
    end
end

using OrbitAl.sparsevec

@testset "Sparse Vectors" begin
    e(i) = unitVec(Rational{Int}, i)
    v = SparseVec(3 => 2//1, 1 => 1//1, 3 => 1//1)
    @test v.poss == [1, 3] && v.vals == [1, 3]
    @test v == e(1) + 3 * e(3)
    @test SparseVec(Rational{Int}[1, 0, 3]) == v == SparseVec(Rational{Int}[1 0 3 0])
    @test SparseVec(Rational{Int}[0 0 0]) == zero(v)
    @test v[3] == 3 && v[2] == 0 && v[7] == 0
    @test length(v) == 2
    @test eltype(v) == eltype(SparseVec{Rational{Int}}) == Rational{Int}
    @test iszero(v - v) && v - v == zero(v)
    @test -v + v == zero(v)
    @test 0 * v == zero(v)
    @test v / 3 == SparseVec(1 => 1//3, 3 => 1//1)
    @test e(1) + e(2) == e(2) + e(1)
    @test length(Set([e(1) + e(2), e(2) + e(1)])) == 1
    @test sprint(show, v) == "e1 + (3)e3"
    @test sprint(show, zero(v)) == "0"
end

using OrbitAl.unionfind

@testset "Union-Find" begin
    n = 10
    pairs = [(1, 6), (4, 9), (6, 2), (8, 3), (9, 10), (2, 7), (3, 5)]
    forest = collect(1:n)
    @test [unite!(forest, i, j) for (i, j) in pairs] == [6, 9, 2, 8, 10, 7, 5]
    @test unite!(forest, 7, 1) == 0
    classes = [[i for i in 1:n if find(forest, i) == r] for r in 1:n if forest[r] == r]
    @test classes == [[1, 2, 6, 7], [3, 5, 8], [4, 9, 10]]

    # linear: x_4 = x_1 + x_3, x_3 = 2 x_2, x_4 + x_5 = x_1
    e(i) = unitVec(Rational{Int}, i)
    parent = Dict{Int, SparseVec{Rational{Int}}}()
    @test unite!(parent, e(4), e(1) + e(3)) == 4
    @test unite!(parent, e(3), 2 * e(2)) == 3
    @test unite!(parent, e(4) + e(5), e(1)) == 5
    @test parent[5] == -2 * e(2)
    @test find(parent, e(4)) == e(1) + 2 * e(2)
    @test unite!(parent, e(5), -2 * e(2)) == 0

    # Union-Find is linear Union-Find for the relations e_i = e_j
    parent = Dict{Int, SparseVec{Rational{Int}}}()
    for (i, j) in pairs
        unite!(parent, e(i), e(j))
    end
    @test all(find(parent, e(i)) == e(find(forest, i)) for i in 1:n)

    # every relation holds in canonical form, and canonical forms are reduced
    m = 30
    for _ in 1:20
        rels = [(SparseVec(rand(1:m) => rand([-2, -1, 1, 2])//1, rand(1:m) => 1//1),
                 SparseVec(rand(1:m) => rand(1:3)//1)) for _ in 1:20]
        parent = Dict{Int, SparseVec{Rational{Int}}}()
        for (u, w) in rels
            unite!(parent, u, w)
        end
        @test all(find(parent, u) == find(parent, w) for (u, w) in rels)
        @test all(!haskey(parent, j) for i in 1:m for j in find(parent, e(i)).poss)
    end
end

using OrbitAl.linear

@testset "Spinning" begin
    # companion matrix of x^5 + x^4 + x^3 + x^2 + x + 1, acting on row vectors
    A = Rational{Int}[0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; -1 -1 -1 -1 -1]
    list = spinning([A], Rational{Int}[0 0 1 0 0], onRight)
    @test vcat(list...) == Rational{Int}[0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; -1 -1 -1 -1 -1; 1 0 0 0 0]
    @test length(spinning([A], Rational{Int}[1 1 0 1 1], onRight)) == 2
    @test length(spinning([A], Rational{Int}[1 0 1 0 1], onRight)) == 1

    # companion matrix of x^20 - 1: e_1 + e_2 spans a subspace of dimension 19
    n = 20
    C = zeros(Rational{Int}, n, n)
    for i in 1:n-1
        C[i, i+1] = 1
    end
    C[n, 1] = 1
    x = zeros(Rational{Int}, 1, n); x[1] = x[2] = 1
    @test length(spinning([C], x, onRight)) == 19

    # spinning with images: the same basis, and the images express the action
    acts(res, aaa) = all(onRight(res.list[j], a) ==
                         sum(c * res.list[l] for (l, c) in zip(res.images[k][j].poss, res.images[k][j].vals))
                         for (k, a) in enumerate(aaa) for j in eachindex(res.list))
    for v in (Rational{Int}[0 0 1 0 0], Rational{Int}[1 1 0 1 1], Rational{Int}[1 0 1 0 1])
        res = spinning_with_images([A], v, onRight)
        @test res.list == spinning([A], v, onRight)
        @test acts(res, [A])
    end
    res = spinning_with_images([A], Rational{Int}[1 1 0 1 1], onRight)
    @test res.images[1] == [SparseVec(2 => 1//1), SparseVec(1 => -1//1, 2 => -1//1)]
    for _ in 1:10
        B = [Rational{Int}.(rand(-1:1, 8, 8)) for _ in 1:2]
        x = Rational{Int}.(rand(-1:1, 1, 8))
        res = spinning_with_images(B, x, onRight)
        @test res.list == spinning(B, x, onRight) && acts(res, B)
    end
end

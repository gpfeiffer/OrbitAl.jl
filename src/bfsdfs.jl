#############################################################################
##
#A  bfsdfs.jl                                                         OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Simple BFS and DFS on trees of nodes that know their children.
##
module bfsdfs

export gcd, gcdex
export Node, BFS, DFS, tree_print

##  Euclid's algorithm as a one-liner
gcd(a, b) = b == 0 ? a : gcd(b, a % b)

##  the extended Euclidean algorithm: (x, y) with gcd(a, b) = a x + b y.
##  recursive, like gcd: if d = b x + (a mod b) y, then d = a y + b (x - (a ÷ b) y)
function gcdex(a, b)
  b == 0 && return (1, 0)
  q, r = divrem(a, b)
  x, y = gcdex(b, r)
  return (y, x - q * y)
end

##  a tree type
struct Node
  id
  next::Array{Node}
end

##  BFS
function BFS(x, visit)
  Q = [x]
  for y in Q
    visit(y)
    append!(Q, y.next)
  end
end

##  DFS
function DFS(x, visit)
  visit(x)
  for z in x.next
    DFS(z, visit)
  end
end

##  print as a tree
function tree_print(x, indent = "", first = true)
    first || print("\n", indent)
    print("-", x.id);
    first = true
    for c in x.next
        tree_print(c, indent * "  ", first)
        first = false
    end
end

end # module

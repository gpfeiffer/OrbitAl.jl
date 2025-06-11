#############################################################################
##
#A  bfsdfs.jl                                                         OrbitAl
#B    by Götz Pfeiffer <goetz.pfeiffer@universityofgalway.ie>
##
#C  Simple BFS and DFS on trees of nodes that know their children.
##
module bfsdfs

export gcd
export Node, BFS, DFS, tree_print

##  Euclid's algorithm in as a one-liner
gcd(a, b) = b == 0 ? a : gcd(b, a % b)

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

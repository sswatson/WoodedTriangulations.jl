module SchnyderWoods

import Graphics2D,
       Stats,
       Base.show

export PlanarMap,
       SchnyderWood,
       schnyder,
       USW,
       draw
       
type PlanarMap
    n::Int64
    nbs::Array{Array{Int64,1},1}
end

type SchnyderWood
    M::PlanarMap
    bluetree::Dict{Int64,Int64}
    redtree::Dict{Int64,Int64}
    greentree::Dict{Int64,Int64}
end

function show(io::IO,S::SchnyderWood)
    print(io,"SchnyderWood(")
    print(io,string(S.M.n-3))
    print(io,")")
end

PlanarMap(n) = PlanarMap(n,[Int64[] for _=1:n])
import Base.length
length(M::PlanarMap) = M.n

lbinom(n,k) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
p0(a,c,n::BigInt) = mod(n+c-a,2) ≠ 0 ? 0.0 : exp(lbinom(n,(n+c-a)/2) + n * log(1/2))
p0float(a,c,n) = exp(lbinom(n,(n+c-a)/2) + n * log(1/2))
p0precise(a,c,n) = Float64(binomial(big(n),big(div(n+c-a,2))) * big(1/2)^n)
function p0(a,c,n::Int64) 
    if mod(n+c-a,2) ≠ 0 || abs(a-c) > n
        return 0.0
    elseif n < 20_000
        return p0float(a,c,n)
    else
        # Use binomial approxmation
        k = div(n+c-a,2)
        return sqrt(2/(pi*n))*exp(-2*(k-n/2)^2/n)
    end
end
p(a,c,n) = p0(a,c,n) - p0(a,-c,n)
q(a,b,n) = p(b,3,n)*p(a,1,n) - p(b,1,n)*p(a,3,n)

function normalize(v)
    return v / sum(v) 
end

function transweights(a,b,n)
    return normalize([q(i,j,n) for i=a-1:2:a+1, j=b-1:2:b+1])
end

function takestep(a,b,n)
    return [(i,j) for i=a-1:2:a+1, j=b-1:2:b+1][:][Stats.sample(Stats.WeightVec(transweights(a,b,n)[:]))] 
end

function contour_to_tree(v::Array{Int64,1})
    n = div(length(v)-1,2)
    M = PlanarMap(n+3)
    stack = Int64[n+1]
    heads = Int64[]
    ones = 0
    blue_edges = Dict{Int64,Int64}()
    for s in diff(v)
        if s == 1
            push!(stack,ones+1)
            ones += 1
        else
            e = pop!(stack)
            push!(heads,e)
            blue_edges[e] = stack[end]
            push!(M.nbs[e],stack[end])
            push!(M.nbs[stack[end]],e)
        end    
    end
    return M, blue_edges, heads
end

function pair_dyck_paths(n)

    X = Int64[1]
    Y = Int64[3]

    for k=1:n
        if mod(k,100_000) == 0 print(k, " ") end
        (a,b) = takestep(X[end],Y[end],(n-k))
        push!(X,a)
        push!(Y,b)
    end

    X = X-1 
    Y = Y-3 
    
    return X,Y
end

function cw(M::PlanarMap,u::Int64,v::Int64)
    nbs = M.nbs[u]
    nbs[(findfirst(nbs,v) % length(nbs))+1]
end

function ccw(M::PlanarMap,u::Int64,v::Int64)
    nbs = M.nbs[u]
    nbs[mod((findfirst(nbs,v)-2),length(nbs))+1]
end

function face(M::PlanarMap,u::Int64,v::Int64)
    f = Int64[u,v]
    while length(f) == 2 || (f[end-1],f[end]) != (u,v)
        push!(f,cw(M,f[end],f[end-1]))
    end
    return f[1:end-1]
end

function edges(M::PlanarMap)
    return vcat([Tuple{Int64,Int64}[(u,v) for v in l] for (u,l) in enumerate(M.nbs)]...)
end

function faces(M::PlanarMap)
    all_edges = edges(M)
    found_edges = zeros(Bool,length(all_edges))
    all_faces = Array{Int64,1}[]
    while sum(found_edges) < length(found_edges)
        f = face(M,all_edges[findfirst(~found_edges)]...)
        for e in zip(f[1:end-1],f[2:end])
            found_edges[findfirst(all_edges,e)] = true
        end
        push!(all_faces,f)
    end
    return all_faces
end

function add_red_tree(M::PlanarMap,
                      blue_edges::Dict{Int64,Int64},
                      Y::Array{Int64,1},
                      heads::Array{Int64,1})
    
    n = div(length(Y)-1,2)
    tails = Int64[]
    ones = 0
    red_edges = Dict{Int64,Int64}()
    
    for s in diff(Y)
        if s == 1
            ones += 1
        end
        if s == -1
            push!(tails, ones ≠ n ? ones+1 : n+2)
        end
    end
    
    hc = deepcopy(heads)

    for t in tails
        headmatch = t == n+2 ? hc[end] : hc[findfirst(hc .≥ t) - 1]
        hc = hc[ hc .≠ headmatch ]
        red_edges[headmatch] = t
        blue_parent = blue_edges[headmatch]
        insert!(M.nbs[headmatch],findfirst(M.nbs[headmatch],blue_parent),t)
        push!(M.nbs[t],headmatch)
    end
    return M, red_edges
end

function find_split(f::Array{Int64,1},blue_edges::Dict{Int64,Int64},red_edges::Dict{Int64,Int64})
    # finds which vertex in each face has two incident faces outgoing red and blue
    for (i,v) in enumerate(f[1:end-1])
        prv = f[i == 1 ? length(f)-1 : i - 1]
        nxt = f[i == length(f)-1 ? 1 : i + 1]
        if v in keys(red_edges) && red_edges[v] == prv && blue_edges[v] == nxt
            return v, prv, nxt
        end
    end
    error("Split not found for face "*string(f))
    return 0,0,0
end

function add_green_tree(M::PlanarMap,blue_edges::Dict{Int64,Int64},red_edges::Dict{Int64,Int64})
    n = maximum(keys(blue_edges))
    green_edges = Dict{Int64,Int64}()
    found = false
    for f in faces(M)            
        if length(f) > 4
            if ~found && (n+1,1) in zip(f[1:end-1],f[2:end]) # detects outer face
                leftedges = f[1:findfirst(f,n+2)-1]
                for (i,w) in enumerate(leftedges)
                    green_edges[w] = n+3
                    insert!(M.nbs[w],findfirst(M.nbs[w],f[i+1]),n+3)
                    unshift!(M.nbs[n+3],w)
                end
                found = true
            else
                v, prv, nxt = find_split(f,blue_edges,red_edges)
                for (i,w) in enumerate(f[1:end-1])
                    if w == v || w == prv || w == nxt
                        continue
                    else
                        green_edges[w] = v
                        insert!(M.nbs[w],findfirst(M.nbs[w],f[i+1]),v)
                        insert!(M.nbs[v],findfirst(M.nbs[v],prv)+1,w)
                    end
                end
            end
        end
    end
    return M, green_edges
end

function USW(n::Integer)
    X,Y = pair_dyck_paths(2n)
    M, blue_edges, heads = contour_to_tree(X)
    M, red_edges = add_red_tree(M,blue_edges,Y,heads)
    M, green_edges = add_green_tree(M,blue_edges,red_edges)
    return SchnyderWood(M,blue_edges,red_edges,green_edges)
end

function descendants(d::Dict{Int64,Int64})
    n = maximum(keys(d))
    desc = zeros(Int64,n)
    for j = 1:n
        k = j
        while true
            k = d[k]
            if k < n+1
                desc[k] += 1
            else
                break
            end
        end
    end
    return desc
end

function flowline(v::Int64,d::Dict{Int64,Int64})
    n = maximum(keys(d))
    branch = Int64[v]
    while branch[end] < n+1
        push!(branch,d[branch[end]])
    end
    return branch
end

function schnyder_one(v::Int64,
                      blue_tree::Dict{Int64,Int64},
                      red_tree::Dict{Int64,Int64},
                      green_tree::Dict{Int64,Int64},
                      Db::Array{Int64,1},
                      Dr::Array{Int64,1},
                      Dg::Array{Int64,1})
    boundary_verts = 0
    interior_verts = 0
    for w in flowline(v,red_tree)[1:end-1]
        interior_verts += Dg[w]
        boundary_verts += 1
    end
    for w in flowline(v,blue_tree)[2:end-1]
        interior_verts += Dg[w] 
        boundary_verts += 1
    end
    return 2*interior_verts + boundary_verts
end

function schnyder_pair(v::Int64,
                       blue_tree::Dict{Int64,Int64},
                       red_tree::Dict{Int64,Int64},
                       green_tree::Dict{Int64,Int64},
                       Db::Array{Int64,1},
                       Dr::Array{Int64,1},
                       Dg::Array{Int64,1})
    return schnyder_one(v,blue_tree,red_tree,green_tree,Db,Dr,Dg), 
           schnyder_one(v,red_tree,green_tree,blue_tree,Dr,Dg,Db) 
end

function schnyder(S::SchnyderWood)
    Db, Dr, Dg = map(descendants,(S.bluetree,S.redtree,S.greentree))
    return [schnyder_pair(v,S.bluetree,S.redtree,S.greentree,Db,Dr,Dg) for v=1:length(S.M)-3]
end
    
function draw(S::SchnyderWood;rot=0.0,linesize=0.125,pointsize=0.001)
    f(z) = cis(0.0)*(z[1] + 0.5*z[2] + im * sqrt(3)/2 * z[2])
    n = length(S.M) - 3
    coords = schnyder(S)
    push!(coords,(0,2n+1))
    push!(coords,(0,0))
    push!(coords,(2n+1,0))
    grlist = Graphics2D.GraphicElement[]
    for (tree,color) in zip((S.bluetree,S.redtree,S.greentree),(Graphics2D.blue,Graphics2D.red,Graphics2D.green))
        for k=1:n
            push!(grlist,Graphics2D.Line([f(coords[k]), f(coords[tree[k]])];color=color,linesize=linesize)) 
        end 
    end
    append!(grlist,Graphics2D.GraphicElement[Graphics2D.Point(f(z);pointsize=pointsize) for z in coords])
    return grlist
end

end # module

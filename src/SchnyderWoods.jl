module SchnyderWoods

import Graphics2D,
       Stats,
       Base.show

export PlanarMap,
       faces,
       face,
       cw,
       ccw,
       edges,
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
    elseif n < 15_000
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
    bluetree = Dict{Int64,Int64}()
    for s in diff(v)
        if s == 1
            push!(stack,ones+1)
            ones += 1
        else
            e = pop!(stack)
            push!(heads,e)
            bluetree[e] = stack[end]
            push!(M.nbs[e],stack[end])
            push!(M.nbs[stack[end]],e)
        end    
    end
    return M, bluetree, heads
end

function pair_dyck_paths(n::Integer;verbose=false)

    X = Int64[1]
    Y = Int64[3]

    power_of_ten = 10^floor(Integer,log10(n/2))

    for k=1:n
        if verbose && k % power_of_ten == 0
	    println(k,"/",n)
	end
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

function edges(S::SchnyderWood)
    all_edges = Tuple{Int64,Int64,String}[]
    for (d,c) in zip((S.bluetree,S.redtree,S.greentree),("blue","red","green"))
        append!(all_edges,sort([[(k,d[k],c) for k in keys(d)];
                                [(d[k],k,c) for k in keys(d)]]))
    end				
    return all_edges
end

function faces_old(M::PlanarMap)
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

function faces(M::PlanarMap)
    rem_edges = Set(edges(M))
    all_faces = Array{Int64,1}[]
    while length(rem_edges) > 0
        f = SchnyderWoods.face(M,pop!(rem_edges)...)
        for e in zip(f[1:end-1],f[2:end])
	    delete!(rem_edges,e)
        end
        push!(all_faces,f)
    end
    return all_faces
end

function faces(S::SchnyderWood)
    n = length(S.M)-3
    colors = Dict([((a,b),c) for (a,b,c) in edges(S)]);
    all_edges = edges(S.M)
    found_edges = zeros(Bool,length(all_edges))
    all_faces = Array{Int64,1}[]
    (u,v) = (n+1,1)
    while sum(found_edges) < length(found_edges)
        if colors[(u,v)] == "blue"
            (u,v) = v,cw(S.M,v,u)
        else
            (u,v) = u,cw(S.M,u,v)
        end
        if found_edges[findfirst(all_edges,(u,v))]
            continue
        end
        f = face(S.M,u,v)
        for e in zip(f[1:end-1],f[2:end])
            found_edges[findfirst(all_edges,e)] = true
        end
        push!(all_faces,f)
    end
    return all_faces
end

function add_red_tree(M::PlanarMap,
                      bluetree::Dict{Int64,Int64},
                      Y::Array{Int64,1},
                      heads::Array{Int64,1})
    
    n = div(length(Y)-1,2)
    tails = Int64[]
    ones = 0
    redtree = Dict{Int64,Int64}()
    
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
        redtree[headmatch] = t
        blue_parent = bluetree[headmatch]
        insert!(M.nbs[headmatch],findfirst(M.nbs[headmatch],blue_parent),t)
        push!(M.nbs[t],headmatch)
    end
    return M, redtree
end

function find_split(f::Array{Int64,1},bluetree::Dict{Int64,Int64},redtree::Dict{Int64,Int64})
    # finds which vertex in each face has two incident faces outgoing red and blue
    for (i,v) in enumerate(f[1:end-1])
        prv = f[i == 1 ? length(f)-1 : i - 1]
        nxt = f[i == length(f)-1 ? 1 : i + 1]
        if v in keys(redtree) && redtree[v] == prv && bluetree[v] == nxt
            return v, prv, nxt
        end
    end
    error("Split not found for face "*string(f))
    return 0,0,0
end

function rotate(a::Array{Int64,1},n::Integer)
    return vcat(a[n:end],a[1:n-1])
end

function findpair(a::Array{Int64,1},m::Integer,n::Integer)
    for k=2:length(a)
        if a[k-1] == m && a[k] == n
            return k
        end
    end   
    return 0 
end


function add_green_tree(M::PlanarMap,bluetree::Dict{Int64,Int64},redtree::Dict{Int64,Int64})
    n = maximum(keys(bluetree))
    greentree = Dict{Int64,Int64}()
    found = false
    for f in faces(M)            
        if length(f) > 4
            if ~found && (n+1,1) in zip(f[1:end-1],f[2:end]) # detects outer face
	        f = rotate(f[1:end-1],findpair(f,n+1,1))
		push!(f,f[1])
                leftedges = f[1:findfirst(f,n+2)-1]
                for (i,w) in enumerate(leftedges)
                    greentree[w] = n+3
                    insert!(M.nbs[w],findfirst(M.nbs[w],f[i+1]),n+3)
                    unshift!(M.nbs[n+3],w)
                end
                found = true
            else
                v, prv, nxt = find_split(f,bluetree,redtree)
		idx = findfirst(M.nbs[v],prv)+1
		fc = rotate(f[1:end-1],findfirst(f,v))[3:end]
                for (i,w) in enumerate(fc[1:end-1])
                    greentree[w] = v
                    insert!(M.nbs[w],findfirst(M.nbs[w],fc[i+1]),v)
                    insert!(M.nbs[v],idx,w)
                end
            end
        end
    end
    return M, greentree
end

function USW(n::Integer)
    X,Y = pair_dyck_paths(2n)
    M, bluetree, heads = contour_to_tree(X)
    M, redtree = add_red_tree(M,bluetree,Y,heads)
    M, greentree = add_green_tree(M,bluetree,redtree)
    return SchnyderWood(M,bluetree,redtree,greentree)
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
                      bluetree::Dict{Int64,Int64},
                      redtree::Dict{Int64,Int64},
                      greentree::Dict{Int64,Int64},
                      Db::Array{Int64,1},
                      Dr::Array{Int64,1},
                      Dg::Array{Int64,1})
    boundary_verts = 0
    interior_verts = 0
    for w in flowline(v,redtree)[1:end-1]
        interior_verts += Dg[w]
        boundary_verts += 1
    end
    for w in flowline(v,bluetree)[2:end-1]
        interior_verts += Dg[w] 
        boundary_verts += 1
    end
    return 2*interior_verts + boundary_verts
end

function schnyder_pair(v::Int64,
                       bluetree::Dict{Int64,Int64},
                       redtree::Dict{Int64,Int64},
                       greentree::Dict{Int64,Int64},
                       Db::Array{Int64,1},
                       Dr::Array{Int64,1},
                       Dg::Array{Int64,1})
    return schnyder_one(v,bluetree,redtree,greentree,Db,Dr,Dg), 
           schnyder_one(v,redtree,greentree,bluetree,Dr,Dg,Db) 
end

function schnyder(S::SchnyderWood)
    n = length(S.M) - 3
    Db, Dr, Dg = map(descendants,(S.bluetree,S.redtree,S.greentree))
    return [[schnyder_pair(v,S.bluetree,S.redtree,S.greentree,Db,Dr,Dg) for v=1:length(S.M)-3];
    	      [(0,2n+1),(0,0),(2n+1,0)]]
end
    
function draw(S::SchnyderWood;
	      rot=0.0,
	      linewidth=0.125,
	      pointsize=0.001,
	      pointcolor="black",
	      includefaces=false,
	      includelabels=false,
	      textsize=1.0)
    ϕ(z) = cis(rot)*(z[1] + 0.5*z[2] + im * sqrt(3)/2 * z[2])
    n = length(S.M) - 3
    coords = schnyder(S)
    colors = Dict([((a,b),c) for (a,b,c) in edges(S)])
    grlist = Graphics2D.GraphicElement[]
    for (tree,color) in zip((S.bluetree,S.redtree,S.greentree),("blue","red","green"))
        for k=1:n
            push!(grlist,Graphics2D.Line([ϕ(coords[k]), ϕ(coords[tree[k]])];color=color,linewidth=linewidth)) 
        end 
    end
    if includelabels
        append!(grlist,Graphics2D.GraphicElement[Graphics2D.Point(ϕ(z);pointsize=1.1*pointsize) for z in coords])
    end
    append!(grlist,Graphics2D.GraphicElement[Graphics2D.Point(ϕ(z);pointsize=pointsize,color=pointcolor) for z in coords])
    if includefaces
        append!(grlist,Graphics2D.GraphicElement[
	Graphics2D.Line(Complex128[ϕ(coords[k]) for k in fc];
                        fill=true,
			linewidth=0.001,	
			fillcolor=facecolor(sort(String[colors[p] 
		   			          for p in zip(fc[1:end-1],fc[2:end])])))
                                                  for fc in sort(faces(S.M),by=length)[1:end-1]])
    end
    if includelabels
        append!(grlist,Graphics2D.GraphicElement[
		 Graphics2D.GraphicText([real(ϕ(z)),imag(ϕ(z))],string(i);textsize=textsize) for (i,z) in enumerate(coords)])
    end
    return grlist
end

function facecolor(s::Array{String,1})
    return sum([1/3 * c for c in s])
end

end # module

__precompile__(true)

module SchnyderWoods

import Graphics2D,
       Stats,
       Base.show

export Cycle,
       PlanarMap,
       OrientedPlanarMap,
       DirectedEdge,
       ColoredEdge,    
       faces,
       face,
       cw,
       ccw,
       edges,
       SchnyderWood,
       schnyder,
       USW,
       draw

#--- TYPES -----------------------------
immutable ColoredEdge
    tail::Int64
    head::Int64
    color::String
end

immutable DirectedEdge
    tail::Int64
    head::Int64
end

immutable Cycle
    elements::Array{Int64,1}
    next::Dict{Int64,Int64} 
    prev::Dict{Int64,Int64}
    function Cycle(a::Array{Int64,1})
        if length(a) == 0
            return new(a,Dict{Int64,Int64}(),Dict{Int64,Int64}())
        end
        return new(a,
            Dict(zip(a,[a[2:end];a[1]])),
            Dict(zip([a[2:end];a[1]],a)))
    end
end

immutable PlanarMap
    n::Int64
    nbs::Array{Cycle,1}
end

immutable OrientedPlanarMap
    M::PlanarMap
    O::Array{Set{Int64},1}
end

immutable SchnyderWood
    M::PlanarMap
    bluetree::Dict{Int64,Int64}
    redtree::Dict{Int64,Int64}
    greentree::Dict{Int64,Int64}
end

function PlanarMap(n::Integer,a::Array{Array{Int64,1},1})
    return PlanarMap(n,map(Cycle,a))
end
#---------------------------------------


#--- 3-ORIENTATION <-> SCHNYDER_WOOD translation --------
DirectedEdge(e::ColoredEdge) = DirectedEdge(e.tail,e.head)

function OrientedPlanarMap(S::SchnyderWood)
    M = S.M
    O = [Set{Int64}() for i=1:M.n]
    for e in colorededges(S)
        push!(O[e.tail],e.head)
    end
    return OrientedPlanarMap(M,O)
end

function SchnyderWood(O::OrientedPlanarMap)
    n = length(O) - 3
    bluetree = Dict{Int64,Int64}()
    redtree = Dict{Int64,Int64}()
    greentree = Dict{Int64,Int64}()
    undiscovered = deepcopy(O.O[1:end-3])
    foundedges = Dict{Tuple{Int64,Int64},Int64}()
    remaining = 3n
    while true
        u = findfirst(x->length(x)>0,undiscovered)
        if u == 0
            break
        end
        v = collect(undiscovered[u])[1]
        p = middlepath(O,u,v,foundedges)
        if p[end] == n + 1 || (p[end] ≤ n && foundedges[(p[end-1],p[end])] == 1)
            for k=1:length(p)-(p[end] == n+1 ? 1 : 2)
                delete!(undiscovered[p[k]],p[k+1])
                remaining -= 1
                bluetree[p[k]] = p[k+1]
                foundedges[(p[k],p[k+1])] = 1
            end
        elseif p[end] == n + 2 || (p[end] ≤ n && foundedges[(p[end-1],p[end])] == 2)
            for k=1:length(p)-(p[end] == n+2 ? 1 : 2)
                delete!(undiscovered[p[k]],p[k+1])
                remaining -= 1
                redtree[p[k]] = p[k+1]
                foundedges[(p[k],p[k+1])] = 2
            end
        else
            for k=1:length(p)-(p[end] == n+3 ? 1 : 2)
                delete!(undiscovered[p[k]],p[k+1])
                remaining -= 1
                greentree[p[k]] = p[k+1]
                foundedges[(p[k],p[k+1])] = 3
            end
        end
    end
    return SchnyderWood(O.M,bluetree,redtree,greentree)
end

function middlepath(O::OrientedPlanarMap,
                    u::Int64,
                    v::Int64,
                    foundedges::Dict{Tuple{Int64,Int64},Int64}=Dict()
                    )
    path = [u,v]
    while path[end] ≤ length(O) - 3 && ~((path[end-1],path[end]) in keys(foundedges))
        C = Cycle(filter(x->x in union(O.O[path[end]],Set(path[end-1])),
                         O.M.nbs[path[end]].elements
                  )
            )
        push!(path,shift(C,path[end-1],2))
    end
    return path
end
#---------------------------------------


#--- CYCLE FUNCTIONS -------------------
import Base.next
function next(C::Cycle,u::Int64)
    return C.next[u]
end
function prev(C::Cycle,u::Int64)
    return C.prev[u]
end

import Base.==
function ==(C::Cycle,D::Cycle)
    a = C.elements
    b = D.elements
    k = findfirst(b,a[1])
    if k == 0
        return false
    else
        return a == [b[k:end]; b[1:k-1]]
    end
end

function pairs(C::Cycle)
    return zip(C.elements,[C.elements[2:end]; C.elements[1:1]])
end

function shift(C::Cycle,u::Integer,k::Integer)
    if k == 0
        return u
    end
    v = u
    if k > 0
        for i=1:k
            v = next(C,v)
        end
    else
        for i=1:-k
            v = prev(C,v)
        end
    end
    return v
end
#-----------------------------------------

#--- DISPLAY FUNCTIONS -------------------
function show(io::IO,S::SchnyderWood)
    print(io,"SchnyderWood(")
    print(io,string(S.M.n-3))
    print(io,")")
end

function show(io::IO,C::Cycle)
    print(io,"Cycle((")
    print(io,string(C.elements))
    print(io,")")
end

function show(io::IO,M::PlanarMap)
    print(io,"PlanarMap(")
    print(io,string(M.n))
    print(io,")")
end

function show(io::IO,O::OrientedPlanarMap)
    print(io,"OrientedPlanarMap(")
    print(io,string(O.M.n))
    print(io,")")
end
#----------------------------------------

###---- FUNCTIONS FOR GENERATING CONDITIONED RANDOM WALK ----------
import Base.length
length(M::PlanarMap) = M.n
length(O::OrientedPlanarMap) = length(O.M)
length(C::Cycle) = length(C.elements)

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
    nbs = [Int64[] for _=1:n+3]
    stack = Int64[n+1]
    heads = Int64[]
    num_ones = 0
    bluetree = Dict{Int64,Int64}()
    for s in diff(v)
        if s == 1
            push!(stack,num_ones+1)
            num_ones += 1
        else
            e = pop!(stack)
            push!(heads,e)
            bluetree[e] = stack[end]
            push!(nbs[e],stack[end])
            push!(nbs[stack[end]],e)
        end    
    end
    return PlanarMap(n+3,nbs), bluetree, heads
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
#----------------------------------------------------

#--- FUNCTIONS FOR GENERATING WOODS -----------------
function cw(M::PlanarMap,u::Int64,v::Int64)
    return next(M.nbs[u],v)
end
function ccw(M::PlanarMap,u::Int64,v::Int64)
    return prev(M.nbs[u],v) 
end

function face(M::PlanarMap,u::Int64,v::Int64)
    f = Int64[u,v]
    while length(f) == 2 || (f[end-1],f[end]) != (u,v)
        push!(f,cw(M,f[end],f[end-1]))
    end
    return f[1:end-1]
end

function edges(M::PlanarMap)
    return vcat([[(u,v) for v in C.elements] for (u,C) in enumerate(M.nbs)]...)
end

function edges(S::SchnyderWood)
    all_edges = Tuple{Int64,Int64,String}[]
    for (d,c) in zip((S.bluetree,S.redtree,S.greentree),("blue","red","green"))
        append!(all_edges,sort([[(k,d[k],c) for k in keys(d)];
                                [(d[k],k,c) for k in keys(d)]]))
    end				
    return all_edges
end

function colorededges(S::SchnyderWood)
    all_edges = ColoredEdge[]
    for (d,c) in zip((S.bluetree,S.redtree,S.greentree),("blue","red","green"))
        append!(all_edges,[ColoredEdge(k,d[k],c) for k in keys(d)]);
    end				
    return all_edges
end

function faces(M::PlanarMap)
    rem_edges = Set(edges(M))
    all_faces = Array{Int64,1}[]
    while length(rem_edges) > 0
        f = face(M,pop!(rem_edges)...)
        for e in zip(f[1:end-1],f[2:end])
	    delete!(rem_edges,e)
        end
        push!(all_faces,f)
    end
    return [Cycle(f[1:end-1]) for f in all_faces]
end

faces(O::OrientedPlanarMap) = faces(O.M)

ccw(O::OrientedPlanarMap,C::Cycle) = ccw(O,[C.elements;C.elements[1:1]])
cw(O::OrientedPlanarMap,C::Cycle) = ccw(O,[C.elements;C.elements[1:1]])

function ccw(O::OrientedPlanarMap,f::Array{Int64,1})
    if maximum(f) > O.M.n - 3
        return false
    end
    for i=1:length(f)-1
        if ~(f[i+1] in O.O[f[i]])
            return false
        end
    end
    return true
end

function cw(O::OrientedPlanarMap,f::Array{Int64,1})
    if maximum(f) > O.M.n - 3
        return false
    end
    for i=1:length(f)-1
        if ~(f[i] in O.O[f[i+1]])
            return false
        end
    end
    return true
end

flipcw(O::OrientedPlanarMap,C::Cycle) = flipcw(O,C.elements)
flipccw(O::OrientedPlanarMap,C::Cycle) = flipccw(O,C.elements)

function flipcw(O::OrientedPlanarMap,f::Array{Int64,1})
    a,b,c = f[1:3]
    delete!(O.O[a],b)
    push!(O.O[a],c)
    delete!(O.O[b],c)
    push!(O.O[b],a)
    delete!(O.O[c],a)
    push!(O.O[c],b)
end

function flipccw(O::OrientedPlanarMap,f::Array{Int64,1})
    a,b,c = f[1:3]
    delete!(O.O[a],c)
    push!(O.O[a],b)
    delete!(O.O[b],a)
    push!(O.O[b],c)
    delete!(O.O[c],b)
    push!(O.O[c],a)
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
    return [Cycle(f[1:end-1]) for f in all_faces]
end

function add_red_tree(M::PlanarMap,
                      bluetree::Dict{Int64,Int64},
                      Y::Array{Int64,1},
                      heads::Array{Int64,1})

    nbs = [a.elements for a in M.nbs]
    
    n = div(length(Y)-1,2)
    tails = Int64[]
    num_ones = 0
    redtree = Dict{Int64,Int64}()
    
    for s in diff(Y)
        if s == 1
            num_ones += 1
        end
        if s == -1
            push!(tails, num_ones ≠ n ? num_ones+1 : n+2)
        end
    end
    
    hc = deepcopy(heads)

    for t in tails
        headmatch = t == n+2 ? hc[end] : hc[findfirst(hc .≥ t) - 1]
        hc = hc[ hc .≠ headmatch ]
        redtree[headmatch] = t
        blue_parent = bluetree[headmatch]
        insert!(nbs[headmatch],findfirst(nbs[headmatch],blue_parent),t)
        push!(nbs[t],headmatch)
    end
    return PlanarMap(length(M),nbs), redtree
end

function find_split(C::Cycle, 
                    bluetree::Dict{Int64,Int64},
                    redtree::Dict{Int64,Int64})
    # finds which vertex in each face has two incident faces outgoing red and blue
    f = [C.elements; C.elements[1:1]]
    for (i,v) in enumerate(f)
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

import Base.getindex, Base.findfirst
getindex(Z::Base.Zip2,k::Integer) = (Z.a[k],Z.b[k])
findfirst(C::Cycle,t::Tuple) = findfirst(pairs(C),t)

function add_green_tree(M::PlanarMap,
                        bluetree::Dict{Int64,Int64},
                        redtree::Dict{Int64,Int64})
    n = maximum(keys(bluetree))
    nbs = [a.elements for a in M.nbs]
    greentree = Dict{Int64,Int64}()
    found = false
    for f in faces(M)            
        if length(f) > 3
            if ~found && (n+1,1) in pairs(f) # detects outer face
	        fa = rotate(f.elements,findfirst(f,(n+1,1)))
		push!(fa,fa[1])
                leftedges = fa[1:findfirst(fa,n+2)-1]
                for (i,w) in enumerate(leftedges)
                    greentree[w] = n+3
                    insert!(nbs[w],findfirst(nbs[w],fa[i+1]),n+3)
                    unshift!(nbs[n+3],w)
                end
                found = true
            else
                fa = f.elements
                v, prv, nxt = find_split(f,bluetree,redtree)
		idx = findfirst(nbs[v],prv)+1
		fc = rotate(fa,findfirst(fa,v))[3:end]
                for (i,w) in enumerate(fc[1:end-1])
                    greentree[w] = v
                    insert!(nbs[w],findfirst(nbs[w],fc[i+1]),v)
                    insert!(nbs[v],idx,w)
                end
            end
        end
    end
    return PlanarMap(length(M),nbs), greentree
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

function schnyder_coordinate(v::Int64,
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

function schnyder_coordinates_2_trees(v::Int64,
                                      bluetree::Dict{Int64,Int64},
                                      redtree::Dict{Int64,Int64},
                                      greentree::Dict{Int64,Int64},
                                      Db::Array{Int64,1},
                                      Dr::Array{Int64,1},
                                      Dg::Array{Int64,1})
    return schnyder_coordinate(v,bluetree,redtree,greentree,Db,Dr,Dg), 
           schnyder_coordinate(v,redtree,greentree,bluetree,Dr,Dg,Db) 
end

function schnyder(S::SchnyderWood)
    n = length(S.M) - 3
    Db, Dr, Dg = map(descendants,(S.bluetree,S.redtree,S.greentree))
    return [vcat((collect(schnyder_coordinates_2_trees(v,S.bluetree,S.redtree,
                                          S.greentree,Db,Dr,Dg))'
             for v=1:length(S.M)-3)...);
    	    [[0 2n+1];[0 0];[2n+1 0]]]
end

function average_schnyder(S::SchnyderWood,M::Int64,N::Int64)
    O = OrientedPlanarMap(S)
    F = filter(x->length(x) == 3, faces(O.M))
    coords = zeros(length(O),2)
    for ctr=1:M
        for j=1:N
            k = rand(1:length(F))
            if ccw(O,F[k])
                flipcw(O,F[k])
            end
            if cw(O,F[k])
                flipccw(O,F[k])
            end
        end
        coords .+= schnyder(SchnyderWood(O))
    end
    return coords/M
end

function facecolor(s::Array{String,1})
    return sum([1/3 * c for c in s])
end
    
function draw(S::SchnyderWood;
	      rot=π,
	      linewidth=1.0,
	      pointsize=0.001,
	      pointcolor="black",
	      includefaces=false,
	      includelabels=false,
              coords=schnyder(S),
	      textsize=1.0)
    ϕ(z) = cis(rot)*(z[1] + 0.5*z[2] + im * sqrt(3)/2 * z[2])
    n = length(S.M) - 3
    colors = Dict([((a,b),c) for (a,b,c) in edges(S)])
    grlist = Graphics2D.GraphicElement[]
    for (tree,color) in zip((S.bluetree,S.redtree,S.greentree),
                            ("blue","red",0.5*"green"))
        for k=1:n
            push!(grlist,Graphics2D.Line([ϕ(coords[k,:]), ϕ(coords[tree[k],:])];
                                         color=color,linewidth=linewidth))
        end 
    end
    if includelabels
        append!(grlist,Graphics2D.GraphicElement[
                       Graphics2D.Point(ϕ(coords[k,:]);
                                        pointsize=1.1*pointsize)
            for k=1:size(coords,1)])
    end
    append!(grlist,Graphics2D.GraphicElement[
        Graphics2D.Point(ϕ(coords[k,:]);
                         pointsize=pointsize,color=pointcolor) for k=1:size(coords,1)])
    if includefaces
        append!(grlist,
                Graphics2D.GraphicElement[
	            Graphics2D.Line(Complex128[ϕ(coords[k,:]) for k in
                                               [fc.elements;fc.elements[1:1]]];
                                    fill=true,
			            linewidth=0.001,	
			            fillcolor=facecolor(sort(String[colors[p]
                                                               for p in pairs(fc)
                                                             ]
                                                        )
                                               )
                     )
                    for fc in sort(faces(S.M),by=length)[1:end-1]
                ]
         )
    end
    if includelabels
        append!(grlist,Graphics2D.GraphicElement[
	    Graphics2D.GraphicText([real(ϕ(coords[i,:])),imag(ϕ(coords[i,:]))],
                                   string(i);textsize=textsize)
            for i=1:size(coords,1)])
    end
    return grlist
end

#------------------------------------------

end # module

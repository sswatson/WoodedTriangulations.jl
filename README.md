# SchnyderWoods

`SchnyderWoods` is a package for generating and display uniform 
graph realizers. 

To install, run 

`Pkg.clone("git://github.com/sswatson/Graphics2D.jl.git")`
`Pkg.clone("git://github.com/sswatson/SchnyderWoods.jl.git")`

An example session: 

```julia
using Graphics2D
using SchnyderWoods

S = USW(1000) # generates a uniform Schnyder wood
showgraphics(draw(S;linewidth=0.125,pointsize=0.001,includefaces=true))
```

![A Schnyder wood](https://github.com/sswatson/SchnyderWoods.jl/blob/master/images/sw.png)

[![Build Status](https://travis-ci.org/sswatson/SchnyderWoods.jl.svg?branch=master)](https://travis-ci.org/sswatson/SchnyderWoods.jl)

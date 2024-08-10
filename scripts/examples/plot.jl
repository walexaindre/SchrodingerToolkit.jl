using GLMakie
using Makie

f(x,y,z) = exp(-(x-2)^2-(y-2)^2-(z-2)^2)

res  = LinRange(-8,8,100)
resv = collect(res)
x = [f(i,j,k) for i in res, j in res, k in res]

colors = to_colormap(:balance)
n = length(colors)
alpha = [ones(n÷3);zeros(n-2*(n÷3));ones(n÷3)]

rs = map(x->(x.r,x.g,x.b),colors)
cmap_alpha = map(x->RGBAf(x[2][1],x[2][2],x[2][3],alpha[x[1]]),enumerate(rs))
volume(resv,resv,resv,x,algorithm=:mip,colormap=cmap_alpha,colorrange=(-1.99,1))
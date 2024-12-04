fig = Figure(size = (800, 600), fontsize = 25,dpi=300)

Ax2D = Axis(fig[1, 1],xlabel = "Points per dimension", ylabel = "Speedup")
x  = [0,30,40,50,80,130]

# Step time 
y1 = [0, 0.15, 2.90, 10.83] #LU CPU
y2 = [0, 0.55, 0.99, 1.58, 4.91, 23.80] #GMRES CPU
y3 = [0, 0.23, 0.56, 1.18] #LU GPU
y4 = [0, 0.22, 0.23, 0.23, 0.42, 1.44] #GMRES GPU

speedupgmres = y2./y4
println("GMRES: ",speedup)

speeduplu = y1./y3
println("LU: ",speedup)

lines!(Ax2D,x[1:length(speeduplu)],speeduplu,label="LU", linestyle = :solid, linewidth = 3.5)
lines!(Ax2D,x[1:length(speedupgmres)],speedupgmres,label="GMRES", linestyle = :dash, linewidth = 3.5)
Ax2D.yticks = WilkinsonTicks(6,k_min=5)
Ax2D.xticks = WilkinsonTicks(6,k_min=5)
axislegend(Ax2D, position = :lt)

save("./speedup3D.png", fig)
fig
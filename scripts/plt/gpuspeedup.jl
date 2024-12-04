fig = Figure(size = (800, 600), fontsize = 25,dpi=300)

Ax2D = Axis(fig[1, 1],xlabel = "Points per dimension", ylabel = "Speedup")
x  = [0,200,400,600,800,1000,1200]

# Step time 
y1 = [0, 0.15, 0.67, 1.62, 2.95, 4.78, 7.00] #LU CPU
y2 = [0, 0.21, 1.65, 6.37, 18.11, 46.54, 102.96] #GMRES CPU
y3 = [0, 0.07, 0.20, 0.40, 0.69, 0.94, 1.41] #LU GPU
y4 = [0, 0.11, 0.37, 1.03, 2.49, 5.14, 9.51] #GMRES GPU


speedupgmres = y2./y4
println("GMRES: ",speedup)

speeduplu = y1./y3
println("LU: ",speedup)

lines!(Ax2D,x,speeduplu,label="LU", linestyle = :dash, linewidth = 3.5)
lines!(Ax2D,x,speedupgmres,label="GMRES", linestyle = :solid, linewidth = 3.5)
Ax2D.yticks = WilkinsonTicks(6,k_min=5)
Ax2D.xticks = WilkinsonTicks(6,k_min=5)
axislegend(Ax2D, position = (0.02,0.99))

save("./speedup.png", fig)
fig
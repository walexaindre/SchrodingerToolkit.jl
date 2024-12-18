using GLMakie
raw""" CPU
        LU & 10.49 & 30 & 0.15 \\
        LU & 29.28 & 40 & 2.90 \\
        LU & 209.71 & 50 & 10.83 \\
        GMRES & 0.35 & 30 & 0.55 \\
        GMRES & 0.35 & 40 & 0.99 \\
        GMRES & 0.53 & 50 & 1.58 \\
        GMRES & 2.34 & 80 & 4.91 \\
        GMRES & 8.83 & 130 & 23.80 \\
"""

raw"""GPU
        LU & 1.88 & 30 & 0.23 \\
        LU & 7.16 & 40 & 0.56 \\
        LU & 24.71 & 50 & 1.18 \\
        GMRES & 0.32 & 30 & 0.22 \\
        GMRES & 0.35 & 40 &  0.23\\
        GMRES & 0.60 & 50 & 0.23 \\
        GMRES & 2.41 & 80 & 0.42 \\
        GMRES & 9.75 & 130 & 1.44 \\
"""

fig = Figure(size = (800, 600), fontsize = 25,dpi=300)

Ax2D = Axis(fig[1, 1],yscale=log2,xlabel = "Points per dimension", ylabel = "Time (s)")
x  = [0,30,40,50,80,130]

# Step time 
y1 = [0, 10.49, 29.28, 209.71] #LU CPU
y2 = [0, 0.35, 0.35, 0.53, 2.34, 8.83] #GMRES CPU
y3 = [0, 1.88, 7.16, 24.71] #LU GPU
y4 = [0, 0.32, 0.35, 0.60, 2.41, 9.75] #GMRES GPU


lines!(Ax2D,x[1:length(y1)],y1,label="LU CPU", linestyle = :dash, linewidth = 3.5)
lines!(Ax2D,x[1:length(y2)],y2,label="GMRES CPU", linestyle = :dot, linewidth = 3.5)
lines!(Ax2D,x[1:length(y3)],y3,label="LU GPU", linestyle = :solid, linewidth = 3.5)
lines!(Ax2D,x[1:length(y4)],y4,label="GMRES GPU", linestyle = :solid, linewidth = 3.5)
ylims!(Ax2D,0.1,200)
xlims!(Ax2D,20,150)
Ax2D.yticks = WilkinsonTicks(5,k_min=4)
Ax2D.xticks = WilkinsonTicks(5,k_min=4)
axislegend(Ax2D, position = (0.99,0.99))

save("./warmup_time3D.png", fig)
fig
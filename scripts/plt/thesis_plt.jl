using GLMakie
raw""" CPU
    \begin{tabular}{lccc}
        \hline
        Solver & Warm up time (s) & Points per dimension & Step time (s) \\\hline\hline
        LU & 1.30 & 200 & 0.15 \\
        LU & 4.98 & 400 & 0.67 \\
        LU & 14.31 & 600 & 1.62 \\
        LU & 29.79 & 800 & 2.95 \\
        LU & 43.09 & 1000 & 4.78 \\
        LU & 69.02 & 1200 & 7.00 \\
        GMRES & 0.10 & 200 & 0.21 \\
        GMRES & 0.27 & 400 & 1.65 \\
        GMRES & 0.41 & 600 & 6.37 \\
        GMRES & 0.84 & 800 & 18.11 \\
        GMRES & 1.67 & 1000 & 46.54 \\
        GMRES & 3.84 & 1200 & 102.96 \\
        \hline
    \end{tabular}
"""

raw"""GPU
    Solver & Warm up time (s) & Points per dimension & Step time (s) \\\hline\hline
    LU & 0.44 & 200 & 0.07 \\
    LU & 1.86 & 400 & 0.20 \\
    LU & 4.33 & 600 & 0.40 \\
    LU & 7.69 & 800 & 0.69 \\
    LU & 12.33 & 1000 & 0.94 \\
    LU & 18.72 & 1200 & 1.41 \\
    GMRES & 0.12 & 200 & 0.11 \\
    GMRES & 0.29 & 400 & 0.37 \\
    GMRES & 0.59 & 600 & 1.03 \\
    GMRES & 0.97 & 800 & 2.49 \\
    GMRES & 1.52 & 1000 & 5.14 \\
    GMRES & 1.87 & 1200 & 9.51 \\
"""


fig = Figure(size = (800, 600), fontsize = 25,dpi=300)

Ax2D = Axis(fig[1, 1],yscale=sqrt,xlabel = "Points per dimension", ylabel = "Time (s)")
x  = [0,200,400,600,800,1000,1200]

# Warm up time
y1 = [0.0, 1.30, 4.98, 14.31, 25.79, 43.09, 69.02 ] #LU CPU
y2 = [0,0.1,0.27,0.41,0.84,1.67,3.84] #GMRES CPU
y3 = [0,0.44,1.86,4.33,7.69,12.33,18.72] #LU GPU
y4 = [0,0.12,0.29,0.59,0.97,1.52,1.87] #GMRES GPU


lines!(Ax2D,x,y1,label="LU CPU", linestyle = :dashdot, linewidth = 3.5)
lines!(Ax2D,x,y2,label="GMRES CPU", linestyle = :dash, linewidth = 3.5)
lines!(Ax2D,x,y3,label="LU GPU", linestyle = :solid, linewidth = 3.5)
lines!(Ax2D,x,y4,label="GMRES GPU", linestyle = :solid, linewidth = 3.5)

Ax2D.yticks = WilkinsonTicks(6,k_min=5)
Ax2D.xticks = WilkinsonTicks(6,k_min=5)
axislegend(Ax2D, position = :lt)

save("./warmup_time.png", fig)

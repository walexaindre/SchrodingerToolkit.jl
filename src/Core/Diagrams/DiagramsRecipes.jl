#Unstable... Makie.jl isn't at release 1.0 yet.

@recipe(SystemND,Memory,Grid,idx) do scene
    Attributes(
        colorscale=identity,
        colormap=theme(scene,:colormap),
        inspectable = theme(scene, :inspectable),
        visible = theme(scene, :visible)
    )
end

function Makie.plot!(Sys::SystemND)
    Grid = Sys[2][]
    Memory = Sys[1][]
    idx = Sys[3][]

    current_state = Memory.current_state

    points = abs2.(current_state[:,idx]) 
    
    points = points |> Array

    z = reshape(points,size(Grid))
    x = get_range(Grid,1)
    y = get_range(Grid,2)

    surface!(Sys,x,y, z;colormap=Sys.colormap,inspectable=Sys.inspectable,visible=Sys.visible)
end


export systemnd
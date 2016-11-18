# using ThreeJS
# import ThreeJS
# using Compat

function ConvertToTuple(p,t)
    points = Array(Tuple{Float64,Float64,Float64},size(p,2))
    for i in 1:size(p,2)
        points[i] = tuple(p[:,i]...)
    end
    
    faces = Array(Tuple{Int64,Int64,Int64},size(t,2))
    for i in 1:size(t,2)
        faces[i] = tuple(t[:,i]...)
    end

    return points,faces
end

function elipspar(p)
    elpar = (0,0,0)
    try
        elpar = ellipsoid_parameters(p)
    catch
        elpar = (0.,0.,0.)
    end

    return [round(i,2) for i in elpar]
end


function geom(frame)
    if frame<=length(memory)       
        geometry(ConvertToTuple(memory[frame][2],memory[frame][3])...)
    else
        info("Please reload browser")
        geometry(ConvertToTuple(memory[end][2],memory[end][3])...)
    end
end


function getabc(points)

    ar = 0
    al = 0
    br = 0
    bl = 0
    cr = 0
    cl = 0
    
    for xkey in 1:size(points,2)
        x = points[:,xkey]
        x[1]<al && (al=x[1])
        x[1]>ar && (ar=x[1])
        x[2]<bl && (bl=x[2])
        x[2]>br && (br=x[2])
        x[3]<cl && (cl=x[3])
        x[3]>cr && (cr=x[3])
    end

    a = (ar - al)/2
    b = (br - bl)/2
    c = (cr - cl)/2

    return a,b,c
end


function main(window)
    
    push!(window.assets,("ThreeJS","threejs"))
    push!(window.assets,"widgets")

    convert(Int,round(3.6))
        
    ff = Signal(length(memory))
    #ThreeJS.title(2,"Viewer for time dependant calculation"),


    vbox(         
                  #"Viewer for time dependant calculation $fname",
                  # hbox(
                  #     ),
                  #vbox(
                  slider(1:length(memory),value=length(memory)) >>> ff,
                  #hbox("frame",),

                  #"data",
                  map(ff) do ff
                  vbox(
                       # Elipsoid parameters $(try [round(i,2) for i in FitEllipsoid(memory[ff][2])] catch NaN end)
                       outerdiv() <<
                       (
                        initscene() <<
                        [
                         ### Pozikridis veids ka attlot meshu
                         #ThreeJS.mesh(0.0, 0.0, 0.0) << [ geom(ff), material(Dict(:color=>"red",:kind=>"phong")) ], 
                         ThreeJS.mesh(0.0, 0.0, 0.0) << [ geom(ff), material(Dict(:color=>"red",:kind=>"lambert")) ], 
                         #ThreeJS.mesh(0.0, 0.0, 0.0) << [ geom(ff), material(Dict(:kind=>"basic",:color=>"black",:wireframe=>true,:wireframeLinewidth=>2)) ], # :wireframe=>true ],
                      pointlight(0.,10.,10.),
                      #ambientlight(0x404040),
                      #            pointlight(0.,-10.,0.),
                      #           pointlight(0.,0.,10.),
                      #             pointlight(0.,0.,-10.),
                      camera(0., 0., 35.,far=1.,fov=10.),
                      ]
                      ),
                      """
                      Time is $(try round(memory[ff][1],2) catch NaN end);          
                      Volume is $(try round(volume(memory[ff][2],memory[ff][3]),4) catch NaN end);
                      Energy is $(try round(memory[ff][4],4) catch NaN end);
                      Normed energy (gamma=1) is $(try round(memory[ff][4]/4/pi/(3/4/pi*volume(memory[ff][2],memory[ff][3]))^(2/3),4) catch NaN end);
                      """,
                      " a,b,c is $(getabc(memory[ff][2]));"
                      )

                  end
    
    )
    
end

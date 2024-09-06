#=

Example of a small fan

This fan is created by relaxation/erosion of a pre-existing plateau

Note how sediments are progressively deposited in the valleys until
the fan is coimpletely full (it has reached steady-state) at which stage
the valleys are progressively emptied

The evolution of the sedimentary flux out of the system and out of the
plateau only is stored in the Fluxes.txt file

=#

using FastScape, WriteVTK
using Random # in order to use the same seed

# Set initial grid size
nx, ny = 101, 201
FastScape_Init()
FastScape_Set_NX_NY(nx,ny)
FastScape_Setup()

xl, yl = 10e3, 20e3 # model dimensions in m
FastScape_Set_XL_YL(xl, yl)
x = range(0, xl, length=nx)
y = range(0, yl, length=ny)



dt = 2e3
FastScape_Set_DT(dt)

# Set erosional parameters
kf = ones(nx,ny)*1e-4
kfsed = 1.5e-4
m = 0.4
n = 1.0
kd = ones(nx,ny)*1.5e-2
kdsed = -1.0
g1 = 1.0
g2 = 1.0
expp = 1.0
FastScape_Set_Erosional_Parameters(kf, kfsed, m, n, kd, kdsed, g1, g2, expp)

# Set BC's - bottom side is fixed only
FastScape_Set_BC(1000)

rng = MersenneTwister(1234);
h = rand!(rng, zeros(nx,ny))    # same random numbers
b = zeros(nx,ny)
FastScape_Init_H(h)
for I in CartesianIndices(h) 
    if y[I[2]] > 0.5*yl
        h[I] += 1000.0  # add a plateau
    end
end
FastScape_Init_H(h)

# set number of time steps and initialize counter istep
nstep = 200

# echo model setup
FastScape_View()

# initializes time step
istep = FastScape_Get_Step()

# To visualize the results in 3D, we need to define 3D arrays:
x2d,y2d = zeros(nx,ny,1), zeros(nx,ny,1)
for I in CartesianIndices(x2d)
    x2d[I], y2d[I] = x[I[1]], y[I[2]]
end
z2d = zeros(nx,ny,1);
basement = zeros(nx,ny,1);

pvd = paraview_collection("Fan")
mkpath("VTK")
while istep<nstep
    global istep, pvd, time, x, y

    # execute step
    FastScape_Execute_Step()
    istep = FastScape_Get_Step()
    
    # outputs h values
    FastScape_Copy_H(h)         # topography
    FastScape_Copy_Basement(b)  # basement
  
    # create VTK file & add it to pvd file
    z2d[:,:,1] = h
    basement[:,:,1] = b
    vtk_grid("VTK/Mountain_$istep", x2d, y2d, z2d) do vtk
        vtk["h [m]"] = z2d
        vtk["Sediment thickness [m]"] = z2d - basement
        vtk["Basement [m]"] = basement
        
        pvd[istep*dt] = vtk
    end
    
    println("step $istep, h-range= $(extrema(h)), nstep=$nstep")    
end

vtk_save(pvd) # save pvd file (open this with Paraview to see an animation)

FastScape_Debug()


FastScape_Destroy()
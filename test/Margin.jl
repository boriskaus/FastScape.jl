#=
Example of the use of the FastScapeInterface
where a square domain (100x100km) is subjected to constant and uniform uplift
of 1 mm/yr while adjacent area (100x100 km) is kept 1000 m below sea level
bottom boundary is at base level, top is nu flux and left and rigth are cyclic
initial random topography

nonlinear erosion law (n=2, m=0.8)
transport coefficient g = 1
marine tranport is activated t

=#
using FastScape, WriteVTK

# Set initial grid size
nx, ny = 101, 151
FastScape_Init()
FastScape_Set_NX_NY(nx,ny)
FastScape_Setup()

xl, yl = 100e3, 150e3 # model dimensions in m
FastScape_Set_XL_YL(xl, yl)
x = range(0, xl, length=nx)
y = range(0, yl, length=ny)

# Timestep
dt = 1e3
FastScape_Set_DT(dt)

# initial topography
h = rand(nx,ny)   # same random numbers
if isfile("h_rand_Margin.txt")
    # use the same random noise for testing
    h = readdlm("h_rand_Margin.txt")
end
b = zeros(nx,ny)
f = zeros(nx,ny)

FastScape_Init_H(h)
for I in CartesianIndices(h) 
    if y[I[2]] < 0.5*yl
        h[I] -= 1000.0  # add a plateau
    end
end
FastScape_Init_H(h)

# Set erosional parameters
kf = ones(nx,ny)*1e-5
kfsed = 1.0e-5
m = 0.8
n = 2.0
kd = ones(nx,ny)*1.0e-2
kdsed = 1e-2
g1 = 1.0
g2 = 1.0
expp = 1.0
FastScape_Set_Erosional_Parameters(kf, kfsed, m, n, kd, kdsed, g1, g2, expp)

# set marine transport parameters
sealevel = 0.0
poro = 0.0
zporo = 1e3
ratio = 0.5
L = 1e2
kds = 3e2
FastScape_Set_Marine_Parameters(sealevel, poro, poro, zporo, zporo, ratio, L, kds, kds/2)

# Set uplift rate
u = ones(nx,ny)*1e-3
for I in CartesianIndices(h) 
    if y[I[2]] < 0.5*yl
        u[I] = 0.0
    end
end

# Set BC's - bottom side is fixed only
FastScape_Set_BC(1000)

# set number of time steps and initialize counter istep
nstep = 200
nfreq = 10

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
Silt = zeros(nx,ny,1);

pvd = paraview_collection("Margin")
mkpath("VTK")
while istep<nstep
    global istep, pvd, time, x, y

    # execute step
    FastScape_Execute_Step()
    istep = FastScape_Get_Step()
    
    # outputs h values
    FastScape_Copy_H(h)         # topography
    FastScape_Copy_Basement(b)  # basement
    FastScape_Copy_F(f)         # silt fraction
    
    if ( (istep/nfreq)*nfreq == istep) 
        # create VTK file & add it to pvd file
        z2d[:,:,1] = h
        basement[:,:,1] = b
        Silt[:,:,1] = f
        vtk_grid("VTK/Mountain_$istep", x2d, y2d, z2d) do vtk
            vtk["h [m]"] = z2d
            vtk["Sediment thickness [m]"] = z2d - basement
            vtk["Basement [m]"] = basement
            vtk["Silt fraction [m]"] = Silt
            
            pvd[istep*dt] = vtk
        end
    
    end
    
    println("step $istep, h-range= $(extrema(h)), nstep=$nstep")    
end

vtk_save(pvd) # save pvd file (open this with Paraview to see an animation)

FastScape_Debug()


FastScape_Destroy()
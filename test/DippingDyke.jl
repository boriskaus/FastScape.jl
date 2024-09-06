#=

problem to test the variability in erodibility (kf)
we assume that a dyke dipping at 30 degrees is buried beneath the landscape
and is progressively exhumed by erosion; for this we use the total erosion
to define the erodibility array kf

=#

using FastScape, WriteVTK
using Random # in order to use the same seed

# run
FastScape_Init()

# Set initial grid size
nx, ny = 201, 201
FastScape_Set_NX_NY(nx,ny)
FastScape_Setup()

xl = 100e3
yl = 100e3
FastScape_Set_XL_YL(xl, yl)
x = range(0, xl, length=nx)
y = range(0, yl, length=ny)

dt = 1e5
FastScape_Set_DT(dt)

rng = MersenneTwister(1234);
h = rand!(rng, zeros(nx,ny))    # same random numbers
FastScape_Init_H(h)

# Set erosional parameters
kf = ones(size(h))*2e-5
kd = ones(size(h))*1e-1
kfsed = -1.0
m = 0.4
n = 1.0
kdsed = -1.0
g = 0.0
p = 1.0
FastScape_Set_Erosional_Parameters(kf, kfsed, m, n, kd, kdsed, g, g, p)

# set uplift rate (uniform while keeping boundaries at base level)
u = ones(size(h))*1e-3
u[:, 1 ] .= 0
u[:, ny] .= 0
u[1,:  ] .= 0
u[nx,: ] .= 0
FastScape_Set_U(u[:])

# Set BC's
FastScape_Set_BC(1111)

# set number of time steps and initialize counter istep
nstep = 50

istep = FastScape_Get_Step()

e = zeros(size(h))
Chi = zeros(nx,ny)

xDyke = xl/10
dxDyke = xl/50
angle = 30
cotana = 1.0/tand(angle)

x2d,y2d = zeros(nx,ny,1), zeros(nx,ny,1)
for I in CartesianIndices(x2d)
    x2d[I], y2d[I] = x[I[1]], y[I[2]]
end
z2d = zeros(nx,ny,1);
Kf = zeros(nx,ny,1);
TotalErosion = zeros(nx,ny,1)


pvd = paraview_collection("DippingDyke")
mkpath("VTK")
while istep<nstep
    global istep, pvd, time, x, y

    FastScape_Copy_Total_Erosion(e)

    # update erodibility where we have the dike
    kf .= 2e-5
    for I in CartesianIndices(h) 
        fac = (x[I[1]]-xDyke-e[I]*cotana-dxDyke)*(x[I[1]]-xDyke-e[I]*cotana+dxDyke) 
        if fac < 0.0
            kf[I] = 1e-5
        end
    end
    FastScape_Set_Erosional_Parameters(kf, kfsed, m, n, kd, kdsed, g, g, p)

    # execute step
    FastScape_Execute_Step()
    
    # get value of time step counter
    istep = FastScape_Get_Step()
    
    # extract solution
    FastScape_Copy_Chi(Chi)
    
    # create VTK file & add it to pvd file
    z2d[:,:,1] = h
    Kf[:,:,1] = kf
    TotalErosion[:,:,1] = e
    vtk_grid("VTK/dike_$istep", x2d, y2d,z2d) do vtk
        vtk["h [m]"] = z2d
        vtk["Erodibility "] = Kf
        vtk["TotalErosion [m]"] = TotalErosion
        pvd[istep*dt] = vtk
    end
   
    #FastScape_VTK(h, 2.0)
    
    x = range(0,xl,length=nx)

    # outputs h values
    FastScape_Copy_H(h)
    
    println("step $istep, h-range= $(extrema(h)), nstep=$nstep")    
end

vtk_save(pvd) # save pvd file (open this with Paraview to see an animation)

FastScape_Debug()

FastScape_Destroy()
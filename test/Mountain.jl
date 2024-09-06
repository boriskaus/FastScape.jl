# Mountain setup
using FastScape 

# run
FastScape_Init()

# Set initial grid size
nx, ny = 101, 101
FastScape_Set_NX_NY(nx,ny)
FastScape_Setup()

xl = 100e3
yl = 100e3
FastScape_Set_XL_YL(xl, yl)

dt = 1e5
FastScape_Set_DT(dt)

h = rand(Float64,nx,ny)
FastScape_Init_H(h)

# Set erosional parameters
kf = ones(size(h))*2e-6
kd = ones(size(h))*1e-1
kfsed = -1.0
m = 0.6
n = 1.5
kdsed = -1.0
g = 0.0
p = -2.0
FastScape_Set_Erosional_Parameters(kf[:], kfsed, m, n, kd[:], kdsed, g, g, p)

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
nstep = 100

istep = FastScape_Get_Step()

Chi = zeros(size(h))

while istep<nstep
    global istep

    # execute step
    FastScape_Execute_Step()
    
    # get value of time step counter
    istep = FastScape_Get_Step()
    
    # extract solution
    FastScape_Copy_Chi(Chi)
    
    # create VTK file
    #FastScape_VTK(h, 2.0)
    
    # outputs h values
    FastScape_Copy_H(h)
    
    println("step $istep, h-range= $(extrema(h)), nstep=$nstep")    
end
FastScape_Debug()


FastScape_Destroy()
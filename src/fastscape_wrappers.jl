
using Libdl
using Fastscapelib_jll

#FSLIB = "../libfastscapelib_fortran.dylib";
#FSLIB = 
function __init__()
    if isfile("libfastscapelib_fortran.dylib")
        global FSLIB = joinpath(pwd(),"libfastscapelib_fortran.dylib")
        println("Using locally compiled version of libfastscapelib_fortran.dylib")
    else
        global FSLIB = Fastscapelib_jll.libfastscapelib_fortran;
        println("Using libfastscapelib_fortran.dylib from Fastscapelib_jll")
    end
end

export FastScape_Init,FastScape_Setup,FastScape_Destroy,FastScape_View,FastScape_Execute_Step
export FastScape_Init_H,FastScape_Init_F,FastScape_Copy_H
export FastScape_Copy_Basement, FastScape_Copy_Total_Erosion, FastScape_Copy_Drainage_Area,FastScape_Copy_Erosion_Rate,FastScape_Copy_Chi, FastScape_Copy_Slope,FastScape_Copy_Curvature,FastScape_Copy_Catchment,FastScape_Copy_F,FastScape_Copy_Lake_Depth
export FastScape_Set_NX_NY,FastScape_Set_XL_YL,FastScape_Set_DT,FastScape_Set_Erosional_Parameters,FastScape_Set_Marine_Parameters
export FastScape_Get_Sizes, FastScape_Get_Step, FastScape_Get_GSSIterations,FastScape_Get_Fluxes
export FastScape_Set_BC,FastScape_Set_U,FastScape_Set_V,FastScape_Set_H,FastScape_Set_All_Layers,FastScape_Set_Basement,FastScape_Set_Precip,FastScape_Set_Tolerance
export FastScape_Debug, FastScape_Strati
export FastScape_Reset_Cumulative_Erosion, FastScape_VTK

"""
    FastScape_Init()
Must be called before any other routine to initialize nx, ny and step
"""
function FastScape_Init()
    ccall((:fastscape_init_, FSLIB), 
        Cvoid, 
        ())
    return
end

"""
    FastScape_Setup()
Must be called to allocate memory for all internal arrays
can only be called once FastScapeSetNXNY has been used to set up nx and ny
"""
function FastScape_Setup()
    ccall((:fastscape_setup_, FSLIB), 
        Cvoid, 
        ())
    return
end

"""
    FastScape_Destroy()
Must be called to deallocate memory
"""
function FastScape_Destroy()
    ccall((:fastscape_destroy_, FSLIB), 
        Cvoid, 
        ())
    return
end

"""
    FastScape_View()
prints to standard output the value of `nx`,`ny`,`nn`,`step`,`xl`,`yl`,`dt`,`K`,`m`,`n`,`kd`,`ibc`
as well as min, mean and max values of `h` and `u`
"""
function FastScape_View()
    ccall((:fastscape_view_, FSLIB), 
        Cvoid, 
        ())
    return
end

"""
    FastScape_Execute_Step()
Executes a single step solving the SPL and diffusion equations
"""
function FastScape_Execute_Step()
    ccall((:fastscape_execute_step_, FSLIB), 
        Cvoid, 
        ())
    return
end

"""
    FastScape_Init_H(h::AbstractArray{Float64})
sets the initial topography (in m) as well as the basement heigh to h
an array of dimension nn(=nx*ny) should be passed
h is double precision of size nn
"""
function FastScape_Init_H(h::AbstractArray{Float64})
    ccall((:fastscape_init_h_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        h)
    return nothing
end

"""
    FastScape_Init_F(F::AbstractArray{Float64})
sets the initial silt fraction to `F`
an array of dimension nn(=nx*ny) should be passed
F is double precision of size nn
"""
function FastScape_Init_F(F::AbstractArray{Float64})
    ccall((:fastscape_init_f_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        F)
    return nothing
end

"""
    FastScape_Copy_H(h::AbstractArray{Float64})
returns the current topographic height (in m)
"""
function FastScape_Copy_H(H::AbstractArray{Float64})
    ccall((:fastscape_copy_h_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        H)
    return 
end

"""
    FastScape_Copy_Basement(b::AbstractArray{Float64})
Returns the current basement height (in m)
"""
function FastScape_Copy_Basement(b::AbstractArray{Float64})
    ccall((:fastscape_copy_basement_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        b)
    return 
end

"""
    FastScape_Copy_Total_Erosion(etot::AbstractArray{Float64})
returns the current cumulative erosion (in m)
"""
function FastScape_Copy_Total_Erosion(etot::AbstractArray{Float64})
    ccall((:fastscape_copy_total_erosion_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        etot)
    return 
end

"""
    FastScape_Copy_Drainage_Area(area::AbstractArray{Float64})
returns the drainage area at each point (in m^2)
"""
function FastScape_Copy_Drainage_Area(area::AbstractArray{Float64})
    ccall((:fastscape_copy_drainage_area_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        area)
    return 
end

"""
    FastScape_Copy_Erosion_Rate(erate::AbstractArray{Float64})
returns the current erosion rate (in m/yr)
"""
function FastScape_Copy_Erosion_Rate(erate::AbstractArray{Float64})
    ccall((:fastscape_copy_erosion_rate_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        erate)
    return 
end

"""
    FastScape_Copy_Chi(chi::AbstractArray{Float64})
returns chi
"""
function FastScape_Copy_Chi(Chi::Array{Float64})
    ccall((:fastscape_copy_chi_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        Chi)
    return 
end

"""
    FastScape_Copy_Slope(slope::AbstractArray{Float64})
returns the current slope
"""
function FastScape_Copy_Slope(slope::AbstractArray{Float64})
    ccall((:fastscape_copy_slope_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        slope)
    return 
end


"""
    FastScape_Copy_Curvature(curvaturep::AbstractArray{Float64})
returns the current curvature
"""
function FastScape_Copy_Curvature(curvaturep::AbstractArray{Float64})
    ccall((:fastscape_copy_curvature_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        curvaturep)
    return 
end

"""
    FastScape_Copy_Catchment(catchment::AbstractArray{Float64})
returns the current catchment
"""
function FastScape_Copy_Catchment(catchment::AbstractArray{Float64})
    ccall((:fastscape_copy_catchment_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        catchment)
    return 
end

"""
    FastScape_Copy_F(F::AbstractArray{Float64})
returns the current surface silt fraction (in m)
"""
function FastScape_Copy_F(F::AbstractArray{Float64})
    ccall((:fastscape_copy_f_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        F)
    return 
end

"""
    FastScape_Copy_Lake_Depth(Lp::AbstractArray{Float64})
returns the lake depth (in m)
"""
function FastScape_Copy_Lake_Depth(Lp::AbstractArray{Float64})
    ccall((:fastscape_copy_lake_depth_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        Lp)
    return 
end

"""
    FastScape_Set_NX_NY(nx::Int64, ny::Int64)

This routine is used to set the resolution of the landscape evolution model. It must be called immediately after `FastScape_Init`.

Arguments:
- `nx` - Resolution or number of grid points in the x-direction (integer)
- `ny` - Resolution or number of grid points in the y-direction (integer)

"""
function FastScape_Set_NX_NY(nnx::Int64,nny::Int64)
    ccall((:fastscape_set_nx_ny_, FSLIB), 
        Cvoid, 
        (Ref{Int32},Ref{Int32}),
        nnx, nny)
    return nothing
end

"""
    FastScape_Set_XL_YL(xl::Float64,yl::Float64)
Sets the value of `xl`,`yl`, the rectangular grid extent (in m)
"""
function FastScape_Set_XL_YL(xl::Float64,yl::Float64)
    ccall((:fastscape_set_xl_yl_, FSLIB), 
        Cvoid, 
        (Ref{Float64},Ref{Float64}),
        xl, yl)
    return nothing
end


"""
    FastScape_Set_DT(dt::Float64)
sets the time step length (in yr)
"""
function FastScape_Set_DT(dt::Float64)
    ccall((:fastscape_set_dt_, FSLIB), 
        Cvoid, 
        (Ref{Float64},),
        dt)
    return nothing
end

"""
    FastScape_Set_Erosional_Parameters (k1::Array{Float64},k2::Float64,m::Float64,n::Float64,kd1::Array{Float64},kd2::Float64,g1::Float64,g2::Float64)

Sets the value of the erosional parameters:
- `k1`,`k2` are rate coefficient in the stream power law (in m^(1-2*m)/yr) for bedrock and sediment respectively
- `m` is the area exponent in the stream power law
- `kd1`, `kd2` are the hillslope transport coefficient or diffusivity (in m^2/yr) for bedrock and sediment respectively
- `g1`, `g2` are the sediment fluvial transport/deposition coefficients (dimensionless) for bedrock and sediment respectively

"""
function FastScape_Set_Erosional_Parameters(kkf::Array{Float64},kkfsed::Float64,mm::Float64,nnn::Float64,kkd::Array{Float64},kkdsed::Float64,gg1::Float64,gg2::Float64,pp::Float64)
    ccall((:fastscape_set_erosional_parameters_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ptr{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64}),
        kkf, kkfsed, mm, nnn, kkd, kkdsed, gg1, gg2, pp)
    return
end

"""
    FastScape_Set_Marine_Parameters(sealevel, poro1, poro2, zporo1, zporo2, ratio, length, kds1, kds2)

Sets the value of the marine transport parameters
- `sl` is sea level (in m)
- `poro1` is surface porosity for silt (dimensionless)
- `poro2` is surface porosity for sand (dimensionless)
- `zporo1` is e-folding porosity depth for silt (in m)
- `zporo2` is e-folding porosity depth for sand (in m)
- `ratio` is the ratio of sand in the incoming flux from the continent (dimensionless)
- `length` is the thickness of the "mixed" surface layer (in m) at the bottom of the ocean
- `kds1` and `kds2` are the marine transport coefficients (diffusivities) for silt and sand respectively (in m^2/yr)
"""
function FastScape_Set_Marine_Parameters(sealevel::Float64, poro1::Float64, poro2::Float64, zporo1::Float64, zporo2::Float64, ratio::Float64, length::Float64, kds1::Float64, kds2::Float64)
    ccall((:fastscape_set_marine_parameters_, FSLIB), 
        Cvoid, 
        (Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64}),
        sealevel, poro1, poro2, zporo1, zporo2, ratio, length, kds1, kds2)
    return
end

"""
    nx, ny = FastScape_Get_Sizes()
returns the value of the grid size
"""
function FastScape_Get_Sizes()
    nx=[Int32(-1)]
    ny=[Int32(-1)]
    
    ccall((:fastscape_get_size_, FSLIB), 
        Cvoid, 
        (Ptr{Int32},Ptr{Int32}),
        nx, ny)
    return nx[1], ny[1]
end

"""
    step = FastScape_Get_Step()
returns the value of the current time step
"""
function FastScape_Get_Step()
    sstep=[Int32(-1)]
    ccall((:fastscape_get_step_, FSLIB), 
        Cvoid, 
        (Ptr{Int32},),
        sstep)
    return sstep[1]
end

"""
    FastScape_Debug()
Display debug info
"""
function FastScape_Debug()
    ccall((:fastscape_debug_, FSLIB), 
        Cvoid, 
        ())
    return
end

"""
    FastScape_Set_BC(ibc::Int64)
sets the boundary conditions
two types are allowed (0 is reflective bc and 1 is fixed base level)
`ibc` should be an integer made of 0 and 1 corresponding to the four boundaries in the
following order: bottom, right, top and left
"""
function FastScape_Set_BC(bc::Int64)
    ccall((:fastscape_set_bc_, FSLIB), 
        Cvoid, 
        (Ref{Int32},),
        bc)
    return nothing
end

"""
    FastScape_Set_U(u::Array{Float64})
sets the uplift velocity/rate (in m/yr)
an array of dimension nn(=nx*ny) should be passed
"""
function FastScape_Set_U(u::Array{Float64})
    ccall((:fastscape_set_u_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        u)
    return nothing
end

"""
    FastScape_Set_V(vx::Array{Float64}, vy:Array{Float64})
sets the x- and y-direction advection velocities/rates (in m/yr)
two array of dimension nn(=nx*ny) should be passed
`vx` and `vy` are double precision of size nn
"""
function FastScape_Set_V(vx::Array{Float64}, vy::Array{Float64})
    ccall((:fastscape_set_v_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},Ptr{Float64}),
        vx,vy)
    return nothing
end

"""
    FastScape_Reset_Cumulative_Erosion()
resets current cumulative erosion
"""
function FastScape_Reset_Cumulative_Erosion()
    ccall((:fastscape_reset_cumulative_erosion_, FSLIB), 
        Cvoid, 
        ())
    return
end

"""
    FastScape_Set_H(h)
resets the surface topography (in m)
"""
function FastScape_Set_H(h::Array{Float64})
    ccall((:fastscape_set_h_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        h)
    return 
end


"""
    FastScape_Set_All_Layers(dhp::Array{Float64})
sets all layers
"""
function FastScape_Set_All_Layers(dhp::Array{Float64})
    ccall((:fastscape_set_all_layers_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        dhp)
    return 
end

"""
    FastScape_Set_Basement(b)
resets the basement topography (in m)
"""
function FastScape_Set_Basement(b::Array{Float64})
    ccall((:fastscape_set_basement_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        b)
    return 
end

"""
    FastScape_Set_Precip(p::Array{Float64})
resets the precipitation rate (in m/yr)
"""
function FastScape_Set_Precip(p::Array{Float64})
    ccall((:fastscape_set_precip_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},),
        p)
    return 
end


"""
     FastScape_Strati(nstepp, nreflectorp, nfreqp, vexp)

"""
function FastScape_Strati(nstepp::Int64, nreflectorp::Int64, nfreqp::Int64, vexp::Float64)
    ccall((:fastscape_strati_, FSLIB), 
        Cvoid, 
        (Ref{Int32},Ref{Int32},Ref{Int32},Ref{Float64}),
        nstepp, nreflectorp, nfreqp, vexp)
    return 
end

"""
    tectonic_flux, erosion_flux, boundary_flux = FastScape_Get_Fluxes()
Gets fluxes
"""
function FastScape_Get_Fluxes()
    ttectonic_flux = [0.0];
    eerosion_flux = [0.0];
    bboundary_flux = [0.0];
    ccall((:fastscape_get_fluxes_, FSLIB), 
        Cvoid, 
        (Ptr{Float64},Ptr{Float64},Ptr{Float64}),
        ttectonic_flux, eerosion_flux, bboundary_flux)
    return ttectonic_flux[1], eerosion_flux[1], bboundary_flux[1]
end

"""
    FastScape_Set_Tolerance(tol_relp, tol_absp, nGSStreamPowerLawMaxp)

Sets tolerances
"""
function FastScape_Set_Tolerance(tol_relp::Float64, tol_absp::Float64, nGSStreamPowerLawMaxp::Int64)
    ccall((:fastscape_set_tolerance_, FSLIB), 
        Cvoid, 
        (Ref{Float64},Ref{Float64},Ref{Int32}),
        tol_relp, tol_absp, nGSStreamPowerLawMaxp)
    return 
end


"""
    nGSSp = FastScape_Get_GSSIterations()
returns the value of the  GSS iterations
"""
function FastScape_Get_GSSIterations()
    nGSSp=[Int32(-1)]
    ccall((:fastscape_get_gssiterations_, FSLIB), 
        Cvoid, 
        (Ptr{Int32},),
        nGSSp)
    return nGSSp[1]
end

function FastScape_VTK(fp::Array{Float64}, vexp::Float64)
    error("This function is not implemented; use the writeVTK.jl package instead")
    return 
end



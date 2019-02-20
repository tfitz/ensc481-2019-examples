module Sphere

export sphere, ellipsoid


function sphere(n=20)

    u = LinRange(0,2*π,n)
    v = LinRange(0,π,n)

    x = cos.(u) * sin.(v)'
    y = sin.(u) * sin.(v)'
    z = ones(n) * cos.(v)'

    return (x,y,z)
end

function ellipsoid(xc,yc,zc,xr,yr,zr; n=20)
    # inputs:
    # xc: center in x
    # yc: center in y
    # zc: center in z
    # xr: radius in x
    # yr: radius in y
    # zr: radius in z
    

    x,y,z = sphere(n)

    x = xr*x .+ xc
    y = yr*y .+ yc
    z = zr*z .+ zc

    return (x,y,z)

end

end

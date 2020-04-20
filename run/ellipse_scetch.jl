# Figure 2 (without labels). This requires working version of Luxor.jl.

import Luxor.ellipse

function ellipse(cpoint::Point, a,b, action=:none;
                 stepvalue=pi/100,
                 vertices=false,
                 tiltangle=0,
                 reversepath=false)
#     a = k/2  # a = ellipse's major axis, the widest part
#     cpoint = midpoint(focus1, focus2)
#     dc = distance(focus1, cpoint)
#     b = sqrt(abs(a^2 - dc^2)) # minor axis, hopefuly not 0
    phi = tiltangle # angle between the major axis and the x-axis
    points = Point[]
    drawing = false
    for t in 0:stepvalue:2pi
        xt = cpoint.x + a * cos(t) * cos(phi) - b * sin(t) * sin(phi)
        yt = cpoint.y + a * cos(t) * sin(phi) + b * sin(t) * cos(phi)
        push!(points, Point(xt, yt))
    end
    vertices ? points : poly(points, action, close=true, reversepath=reversepath)
end

function drawcolumn(a,b_e,h,color)
    setcolor(color...)

    a_gc = a*√(1-h^2/b^2)
    b_gc = b_e/a*a_gc
    perspective=(b-h/2)/b

    #top
    epoly = ellipse(Point(0,-h),a_gc,b_gc*perspective, vertices=true)
    setdash("solid")
    poly(epoly[1:end÷2], :stroke,  close=false)
    setdash("dash")
    poly(epoly[end÷2+1:end], :stroke,  close=false)

    #center
    epoly = ellipse(cp,a_gc,b_gc, vertices=true)
    setdash("solid")
    poly(epoly[1:end÷2], :stroke, close=false)
    setdash("dash")
    poly(epoly[end÷2+1:end], :stroke,  close=false)


    #bottom
    epoly = ellipse(Point(0,h),a_gc,b_gc*perspective, vertices=true)
    setdash("solid")
    poly(epoly[1:end÷2], :stroke,  close=false)
    setdash("dash")
    poly(epoly[end÷2+1:end], :stroke,  close=false)


    setdash("solid")
    line(Point(a_gc,-h),Point(a_gc,h),:stroke)
    line(Point(-a_gc,-h),Point(-a_gc,h),:stroke)
end

a,b = 200,120
resa = Int(2.1*a)
resb = Int(3*b)
@show resa, resb


C2=(0.17254901960784313, 0.6274509803921569, 0.17254901960784313)
forestgreen=(0.13333333333333333, 0.5450980392156862, 0.13333333333333333)
C0=(31/255,119/255,180/255) #color






@pdf begin
    fontsize(16)
    fontface("CMU Serif")
    f1 = Point(-150, 0)
    f2 = Point(150, 0)
    cp = Point(0,0)


    #rotation axis:
    Luxor.arrow(Point(0,b+b/2),Point(0,-b-b/2))

    Luxor.text("u", Point(10,-b*1.3))

    #main ellipse
    epoly1 = ellipse(cp,a,b, vertices=true)

    #equatorial ellipse
    b_e=b/4
    epoly2 = ellipse(cp,a,b_e, vertices=true)

    #meridional ellipse
    a_m = a/8
    epoly3 = ellipse(cp,a_m,b, vertices=true)


    #draw background dashed half ellipses

    setdash("dash")
    setcolor(forestgreen)
    poly(epoly2[end÷2+1:end], :stroke,  close=false)
    setcolor("black")
    poly(epoly3[[3*end÷4+1:end;1:end÷4]], :stroke,  close=false)

    #fill equatorial plane
    setcolor(forestgreen)
    setopacity(0.3)
    poly(epoly2, :fill, close=true)
    setopacity(1)
    setcolor("black")
    ##semi axes
    setdash("solid")
    Luxor.arrow(O,Point(-a_m,b_e))
    Luxor.arrow(O, Point(0,-b))
    Luxor.arrow(O, Point(a,0))


    ##geostrophic column

    drawcolumn(a,b_e,0.85*b,C0)




    #draw foreground solid half ellipses
    setdash("solid")
    setcolor("black")
    poly(epoly1, :stroke,  close=true)
    setcolor(forestgreen)
    poly(epoly2[1:end÷2], :stroke,  close=false)
    setcolor("black")
    poly(epoly3[[end÷4+1:end÷2;end÷2+1:3*end÷4]], :stroke,  close=false)

end 420 360 "../figs/aligned_ellipsoid_trunc.pdf"

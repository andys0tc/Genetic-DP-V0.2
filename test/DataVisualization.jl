using DataFrames, Gadfly, CSV

df = DataFrame(CSV.File("/home/andrea/Julia/GenDP-main v02/test/raw_results_VC.csv"))



sizemap(p::Float64; min=1mm, max=3mm)= min + p/15*(max-min)
#sizemap(p::Float64,min=1mm,max=3mm)=(p*max)/10

#palettef = Scale.lab_gradient("red", "orange", "blue")
coord = Coord.cartesian(xmin=0.5, ymin=0, xmax=1.0, ymax=200)
plot(df,coord,Scale.size_area(sizemap, minvalue=1, maxvalue=5),
    layer(df,x =:NaiveGAAppRatio, y=:TreeSize,Theme(highlight_width=0pt), 
        Geom.point,color=[colorant"grey"],alpha=[0.25],size=:NumPopu,shape =[Shape.circle]),

    layer(df,x =:HalfGenDPAppRatio, y=:TreeSize,Theme(highlight_width=0pt),
        Geom.point,color=[colorant"purple"],alpha=[0.25],size=:NumPopu,shape =[Shape.circle]),

     layer(df,x =:GenDPAppRatio, y=:TreeSize,Theme(highlight_width=0pt),
        Geom.point,color=[colorant"blue"],alpha=[0.25],size=:NumPopu,shape =[Shape.circle]))



# Heat Transfer Equation (2D Finite-Difference Time-Domain approximation)
# www.overfitting.net
# https://www.overfitting.net/2021/02/transferencia-de-calor-por-elementos_10.html

library(png)
library(stringr)

object="walls"


# Read simulation parameters
heatparams=read.table("heatparams.csv", header=TRUE, sep=",")
dx=heatparams$value[heatparams$desc=='cell size (m)']
dt=heatparams$value[heatparams$desc=='time step (s)']
N=as.integer(heatparams$value[heatparams$desc=='number of iterations'])
NSNAPSHOTS=as.integer(heatparams$value[heatparams$desc=='number of snapshots'])

# Read objects
# Temp (K or ºC), alpha (m2/s), k (w/(m*K)), q (w/m3), type, desc
# Vales for type: 'boundary', 'solid', 'fluid', 'source'
heatobjparams=read.table("heatobjects.csv", header=TRUE, sep=",")
heatobjects=readPNG("heatobjects.png")

plot(as.raster(heatobjects), interpolate=F)  # display objects
heatobjectsunique=heatobjects[,,3]+heatobjects[,,1]*2+heatobjects[,,2]*4
writePNG(heatobjectsunique/max(heatobjectsunique),
         paste0(object,"_gray.png"))


# Stability condition (if UNSTABLE reduce dt or increase dx)
# Threshold is 1/2 for 1D, 1/4 for 2D, 1/6 for 3D
r=max(heatobjparams$alpha)*dt/dx^2
print(paste0("r=", r, " -> ", ifelse(r<=1/4, "STABLE", "UNSTABLE"),
             " (",round(r/(1/4)*100),"%)"))


# Obtain unique colours and round to material position
INITMARK=heatobjectsunique[1,1]
i=2
while (heatobjectsunique[1,i] != INITMARK) i=i+1
NMAT=i-2  # number of materials (colours) defined
print(paste0(NMAT, " materials defined"))

colours=heatobjectsunique[1,2:(NMAT+1)]  # just colours (in first row)
heatobjectsunique=heatobjectsunique[-1,]  # drop first row
NROW=nrow(heatobjectsunique)
NCOL=ncol(heatobjectsunique)

# Row/Col precalculation for 'insulate' rectangular objects
MINROW=array(0,NMAT)
MAXROW=MINROW
MINCOL=MINROW
MAXCOL=MINROW

lst=list()
for (i in 1:NMAT) {
    FLAG_INSULATE=(heatobjparams$type[i]=='insulate')
    indices=which(heatobjectsunique==colours[i], arr.ind=FLAG_INSULATE)
    heatobjectsunique[indices]=i
    lst[[i]]=indices  # create indexing list for each object
    
    # Corners of 'insulate' rectangular object
    if (FLAG_INSULATE) {
        MINROW[i]=min(lst[[i]][,1])
        MAXROW[i]=max(lst[[i]][,1])
        
        MINCOL[i]=min(lst[[i]][,2])
        MAXCOL[i]=max(lst[[i]][,2])        
    }
}
plot(as.raster(heatobjectsunique/NMAT), interpolate=F)


# Plot scene contour
heatobjectscontour=heatobjectsunique*0
heatobjectscontour[2:(NROW-1),2:(NCOL-1)]=
    abs(heatobjectsunique[1:(NROW-2),2:(NCOL-1)] -
        heatobjectsunique[2:(NROW-1),2:(NCOL-1)]) +
    abs(heatobjectsunique[2:(NROW-1),1:(NCOL-2)] -
        heatobjectsunique[2:(NROW-1),2:(NCOL-1)])
heatobjectscontour[heatobjectscontour != 0]=1
writePNG(heatobjectscontour, paste0(object,"_contour.png"))


# Precalculate working arrays
# alpha=k/(rho*cp), rho*cp=k/alpha
Temp=heatobjectsunique*0
for (i in 1:NMAT) {
    Temp[lst[[i]]]=Temp[lst[[i]]]+1
}
if ((min(Temp) != 1) || (max(Temp) != 1)) print("ERROR!") else print("OK")

k=Temp
alpha=Temp
for (i in 1:NMAT) {
    Temp[lst[[i]]] =heatobjparams$Temp[i]
    k[lst[[i]]]    =heatobjparams$k[i]
    alpha[lst[[i]]]=heatobjparams$alpha[i]
}
rhocp=k/alpha  # we won't use alpha, just k and rhocp (=rho*cp)



# Time domain T iteration using an explicit FD scheme
# valid for heterogeneous conductivity media
MINT=min(heatobjparams$Temp)
MAXT=max(heatobjparams$Temp)

# Material (object) of special interest
MAINMAT=5
tempe=c()

SKIP=round(N/NSNAPSHOTS)
for (j in 0:N) {
    MINT=min(Temp)
    MAXT=max(Temp)
    
    # Snapshot T distribution
    if (j %% SKIP==0) {
        # Save PNG
        nombre=paste0("heattransfer_",
                      str_pad(j, nchar(N), pad='0'), ".png")
        writePNG(((Temp-MINT)/(MAXT-MINT))^0.2, nombre)
        
        # Print AVG T per material
        txt=paste0("Iter ", j, "/", N, ": ")
        for (i in 1:NMAT) {
            txt=paste0(txt, ifelse(i==1,""," - "), heatobjparams$desc[i], ": ",
                       round(min(Temp[lst[[i]]]), 1), "/",
                       round(mean(Temp[lst[[i]]]), 1), "/",
                       round(max(Temp[lst[[i]]]), 1))
        }
        print(txt)
        tempe=c(tempe, mean(Temp[lst[[MAINMAT]]]))
    }

    # Iterate T for the whole grid using the standard formula
    Temp[2:(NROW-1),2:(NCOL-1)] = Temp[2:(NROW-1),2:(NCOL-1)] +
        dt/(rhocp[2:(NROW-1),2:(NCOL-1)] * dx^2) *
        (
            (   k[3: NROW   ,2:(NCOL-1)] +    k[2:(NROW-1),2:(NCOL-1)])/2 *
            (Temp[3: NROW   ,2:(NCOL-1)] - Temp[2:(NROW-1),2:(NCOL-1)]) +
            (   k[1:(NROW-2),2:(NCOL-1)] +    k[2:(NROW-1),2:(NCOL-1)])/2 *
            (Temp[1:(NROW-2),2:(NCOL-1)] - Temp[2:(NROW-1),2:(NCOL-1)]) +
                
            (   k[2:(NROW-1),3: NCOL   ] +    k[2:(NROW-1),2:(NCOL-1)])/2 *
            (Temp[2:(NROW-1),3: NCOL   ] - Temp[2:(NROW-1),2:(NCOL-1)]) +
            (   k[2:(NROW-1),1:(NCOL-2)] +    k[2:(NROW-1),2:(NCOL-1)])/2 *
            (Temp[2:(NROW-1),1:(NCOL-2)] - Temp[2:(NROW-1),2:(NCOL-1)])             
        )
    
    # Special materials: boundaries, fluids, internal sources and insulates
    for (i in 1:NMAT) {
        if (heatobjparams$type[i]=='boundary') {
            Temp[lst[[i]]] = heatobjparams$Temp[i]  # constant T
        } else if (heatobjparams$type[i]=='fluid') {
            Temp[lst[[i]]] = mean(Temp[lst[[i]]])  # instantaneous convection
        } else if (heatobjparams$type[i]=='source') {
            Temp[lst[[i]]] = Temp[lst[[i]]] +
                dt/rhocp[lst[[i]]]*heatobjparams$q[i]  # heat source
        } else if (heatobjparams$type[i]=='insulate') {  # copy T
            Temp[lst[[i]]]=heatobjparams$Temp[i]  # reset T
            
            # Copy T along top and bottom
            Temp[MINROW[i], MINCOL[i]:MAXCOL[i]]=
                Temp[MINROW[i]-1, MINCOL[i]:MAXCOL[i]]
            Temp[MAXROW[i], MINCOL[i]:MAXCOL[i]]=
                Temp[MAXROW[i]+1, MINCOL[i]:MAXCOL[i]]
            
            # Copy T along left and right
            Temp[MINROW[i]:MAXROW[i], MINCOL[i]]=
                Temp[MINROW[i]:MAXROW[i], MINCOL[i]-1]
            Temp[MINROW[i]:MAXROW[i], MAXCOL[i]]=
                Temp[MINROW[i]:MAXROW[i], MAXCOL[i]+1]
            
            # Refine 4 corners
            Temp[MINROW[i], MINCOL[i]]=
                (Temp[MINROW[i]-1, MINCOL[i]]+Temp[MINROW[i], MINCOL[i]-1])/2
            Temp[MAXROW[i], MINCOL[i]]=
                (Temp[MAXROW[i]+1, MINCOL[i]]+Temp[MAXROW[i], MINCOL[i]-1])/2
            Temp[MINROW[i], MAXCOL[i]]=
                (Temp[MINROW[i]-1, MAXCOL[i]]+Temp[MINROW[i], MAXCOL[i]+1])/2
            Temp[MAXROW[i], MAXCOL[i]]=
                (Temp[MAXROW[i]+1, MAXCOL[i]]+Temp[MAXROW[i], MAXCOL[i]+1])/2
        }
    }
}

# Evolution of T in object of special interest
plot(seq(0,dt*N/60,length.out=length(tempe)), tempe, type='l', col='red',
     xlab='Time (min)', ylab='T (ºC/K)', main=heatobjparams$desc[MAINMAT])



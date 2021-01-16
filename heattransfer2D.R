# Heat Transfer Equation (2D Finite-Difference Time-Domain approximation)
# www.overfitting.net

library(png)
library(stringr)

object="cube"


# Read simulation parameters
heatparams=read.table("heatparams.csv", header=TRUE, sep=",")
dx=heatparams$value[heatparams$desc=='cell size (m)']
dt=heatparams$value[heatparams$desc=='time step (s)']
N=as.integer(heatparams$value[heatparams$desc=='number of iterations'])
NSNAPSHOTS=as.integer(heatparams$value[heatparams$desc=='number of snapshots'])

# Read objects
heatobjparams=read.table("heatobjects.csv", header=TRUE, sep=",")
heatobjects=readPNG("heatobjects.png")

plot(as.raster(heatobjects), interpolate=F)  # display objects
heatobjectsunique=heatobjects[,,3]+heatobjects[,,1]*2+heatobjects[,,2]*4
writePNG(heatobjectsunique/max(heatobjectsunique),
         paste0(object,"_gray.png"))


# Stability condition (if UNSTABLE reduce dt or increase dx)
# Threshold is 1/2 for 1D, 1/4 for 2D, 1/6 for 3D
r=dt*max(heatobjparams$alpha)/dx^2
print(paste0("r=", r, " -> ", ifelse(r<=1/4, "STABLE", "UNSTABLE"),
             " (",round(r/(1/4)*100),"%)"))


# Obtain unique colours and round to object position
INITMARK=heatobjectsunique[1,1]
i=2
while (heatobjectsunique[1,i] != INITMARK) i=i+1
NOBJECTS=i-2  # number of objects (colours) defined
print(paste0(NOBJECTS, " objects defined"))

colours=heatobjectsunique[1,2:(NOBJECTS+1)]  # just colours (in first row)
heatobjectsunique=heatobjectsunique[-1,]  # drop first row
NROW=nrow(heatobjectsunique)
NCOL=ncol(heatobjectsunique)

lst=list()
for (i in 1:NOBJECTS) {
    indices=which(heatobjectsunique==colours[i])
    heatobjectsunique[indices]=i
    lst[[i]]=indices  # create indexing list for each object
}
plot(as.raster(heatobjectsunique/NOBJECTS), interpolate=F)


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
k=Temp
alpha=Temp
for (i in 1:NOBJECTS) {
    Temp[lst[[i]]]=heatobjparams$Temp[i]
    k[lst[[i]]]=heatobjparams$k[i]
    alpha[lst[[i]]]=heatobjparams$alpha[i]
}
rhocp=k/alpha  # won't use alpha, just k and rhocp=rho*cp


# Time domain T iteration using an explicit FD scheme
# valid for heterogeneous conductivity media
MINT=min(heatobjparams$Temp)
MAXT=max(heatobjparams$Temp)

SKIP=round(N/NSNAPSHOTS)
for (j in 0:N) {
    # Snapshot T distribution
    if (j %% SKIP==0) {
        nombre=paste0("heattransfer_",
                      str_pad(j, nchar(N), pad='0'), ".png")
        writePNG((Temp-MINT)/(MAXT-MINT), nombre)
    }

    # Iterate T for the whole grid
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
    
    # Reset T on boundaries and assume instantaneous convection on fluids
    for (i in 1:NOBJECTS) {
        if (heatobjparams$type[i]=='boundary') {
            Temp[lst[[i]]]=heatobjparams$Temp[i]  # reset T
        } else if (heatobjparams$type[i]=='fluid') {
            Temp[lst[[i]]]=mean(Temp[lst[[i]]])  # average T (=E conservation)           
        }
    }
}

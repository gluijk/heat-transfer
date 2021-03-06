# Heat Transfer Equation (2D Finite-Difference Time-Domain approximation)
# www.overfitting.net
# https://www.overfitting.net/2021/02/transferencia-de-calor-por-elementos_10.html

library(png)
library(stringr)
library(rgl)
library(viridis)

OBJPARAMS="heatobjects.csv"  # objects thermal definition
OBJCOLOURS="heatobjects.png"  # simulation geometry definition
SIMPARAMS="heatparams.csv"  # simulation parameters
object="flux"  # simulation name


# Read simulation parameters
heatparams=read.table(SIMPARAMS, header=TRUE, sep=",")
dx=heatparams$value[heatparams$desc=='cell size (m)']
dt=heatparams$value[heatparams$desc=='time step (s)']
N=as.integer(heatparams$value[heatparams$desc=='number of iterations'])
NSNAPSHOTS=as.integer(heatparams$value[heatparams$desc=='number of snapshots'])
gamma=heatparams$value[heatparams$desc=='gamma']

# Read objects
# Temp (K or �C), alpha (m2/s), k (w/(m*K)), q (w/m3), type, desc
# Valid values for type:
# 'solid', 'source', 'fluid', 'boundary', 'insulate', 'isoflux'
heatobjparams=read.table(OBJPARAMS, header=TRUE, sep=",")
heatobjects=readPNG(OBJCOLOURS)

plot(as.raster(heatobjects), interpolate=F)  # display objects
heatobjectsunique=heatobjects[,,3]+heatobjects[,,1]*2+heatobjects[,,2]*4


# Stability condition (if UNSTABLE reduce dt and/or increase dx)
# Threshold is 1/2 for 1D, 1/4 for 2D, 1/6 for 3D
r=max(heatobjparams$alpha)*dt/dx^2
print(paste0(ifelse(r<=1/4, "STABILITY", "ERROR"),
            ": r=", r, " (",round(r/(1/4)*100),"%)",
            " -> ", ifelse(r<=1/4, "STABLE",
            "UNSTABLE (reduce dt and/or increase dx)")))


# Obtain unique colours and round to material position
INITMARK=heatobjectsunique[1,1]
i=2
while (heatobjectsunique[1,i] != INITMARK) i=i+1
NMAT=i-2  # number of materials (colours) defined
if (NMAT != nrow(heatobjparams)) {
    print(paste0("ERROR: ", nrow(heatobjparams),
                 " materials ('", OBJPARAMS, "') and ", NMAT,
                 " colours ('", OBJCOLOURS, "') mismatch"))
} else print(paste0("OK: ", NMAT, " materials defined"))


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
    FLAG_IDEAL=(heatobjparams$type[i]=='insulate' |
                heatobjparams$type[i]=='isoflux')
    indices=which(heatobjectsunique==colours[i], arr.ind=FLAG_IDEAL)
    if (length(indices)) {  # there is at least one pixel
        heatobjectsunique[indices]=i
        lst[[i]]=indices  # create indexing list for each object
        
        # Precalculate corners of 'insulate' rectangular object
        if (FLAG_IDEAL) {
            MINROW[i]=min(lst[[i]][,1])
            MAXROW[i]=max(lst[[i]][,1])
            MINCOL[i]=min(lst[[i]][,2])
            MAXCOL[i]=max(lst[[i]][,2])
            
            # Check if insulate material is made of a single rectangle
            if ((MAXROW[i]-MINROW[i]+1)*(MAXCOL[i]-MINCOL[i]+1)
                !=nrow(lst[[i]])) print(paste0(
                    "ERROR: insulate/isoflux material ", i, " ('",
                    heatobjparams$desc[i],"') is not a single rectangle"))
        }        
    } else {
        print(paste0("WARNING: material ", i, " ('",
                     heatobjparams$desc[i],"') doesn't appear in the scene"))
    }
}


# Plot scene contour
heatobjectscontour=heatobjectsunique*0
heatobjectscontour[2:(NROW-1),2:(NCOL-1)]=
    abs(heatobjectsunique[1:(NROW-2),2:(NCOL-1)] -
        heatobjectsunique[2:(NROW-1),2:(NCOL-1)]) +
    abs(heatobjectsunique[2:(NROW-1),1:(NCOL-2)] -
        heatobjectsunique[2:(NROW-1),2:(NCOL-1)])
heatobjectscontour[heatobjectscontour != 0]=1
writePNG(heatobjectscontour, paste0(object,"_contour.png"))


# Precalculate arrays
# alpha=k/(rho*cp), rho*cp=k/alpha
Temp=heatobjectsunique*0
for (i in 1:NMAT) {
    Temp[lst[[i]]]=Temp[lst[[i]]]+1
}
if (min(Temp) != 1) {
    print(paste0("ERROR: wrong colours in 'heatobjects.png', ",
        "some part of the scene has an unknown colour"))
    writePNG(Temp, paste0(object,"heatobjects_ERROR.png"))
} else print("OK: all parts in the scene have a defined colour")

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

SKIP=round(N/NSNAPSHOTS)
T_Evol=array(0, c(NMAT, NSNAPSHOTS+1))
N_OUT=0
for (j in 0:N) {
    MINT=min(Temp)
    MAXT=max(Temp)
    
    # Snapshot T distribution
    if (j %% SKIP==0) {
        # Save PNG
        nombre=paste0("heattransfer_",
                      str_pad(N_OUT, nchar(NSNAPSHOTS), pad='0'), "_",
                      str_pad(j, nchar(N), pad='0'), ".png")
        writePNG(((Temp-MINT)/(MAXT-MINT))^(1/gamma), nombre)
        N_OUT=N_OUT+1
        
        # Print and save AVG T per material
        txt=paste0("Iter ", j, "/", N, ": ")
        for (i in 1:NMAT) {
            txt=paste0(txt, ifelse(i==1,""," - "), heatobjparams$desc[i], ": ",
                       round(min(Temp[lst[[i]]]), 1), "/",
                       round(mean(Temp[lst[[i]]]), 1), "/",
                       round(max(Temp[lst[[i]]]), 1))
            T_Evol[i,N_OUT]=mean(Temp[lst[[i]]])
        }
        print(txt)
    }

    # Iterate T for the whole grid using the standard equation for 'solid'
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
    
    # Refine special materials
    for (i in 1:NMAT) {
        if (heatobjparams$type[i]=='source') {
            Temp[lst[[i]]] = Temp[lst[[i]]] +
                dt/rhocp[lst[[i]]]*heatobjparams$q[i]  # heat source
            
        } else if (heatobjparams$type[i]=='fluid') {
            Temp[lst[[i]]] = mean(Temp[lst[[i]]])  # instantaneous convection
            
        } else if (heatobjparams$type[i]=='boundary') {
            Temp[lst[[i]]] = heatobjparams$Temp[i]  # constant T

        } else if (heatobjparams$type[i]=='insulate') {  # copy T
            Temp[lst[[i]]] = heatobjparams$Temp[i]  # reset T
            
            # Copy T along bottom and top
            if (MINROW[i]>1) Temp[MINROW[i], MINCOL[i]:MAXCOL[i]] =
                    Temp[MINROW[i]-1, MINCOL[i]:MAXCOL[i]]
            if (MAXROW[i]<NROW) Temp[MAXROW[i], MINCOL[i]:MAXCOL[i]] =
                    Temp[MAXROW[i]+1, MINCOL[i]:MAXCOL[i]]
            
            # Copy T along left and right
            # (left/right prevails over bottom/top)
            if (MINCOL[i]>1) Temp[MINROW[i]:MAXROW[i], MINCOL[i]] =
                    Temp[MINROW[i]:MAXROW[i], MINCOL[i]-1]
            if (MAXCOL[i]<NCOL) Temp[MINROW[i]:MAXROW[i], MAXCOL[i]] =
                    Temp[MINROW[i]:MAXROW[i], MAXCOL[i]+1]

        } else if (heatobjparams$type[i]=='isoflux') {  # transfer heat flux
            Temp[lst[[i]]] = heatobjparams$Temp[i]  # reset T
            
            # Transfer heat flux along bottom and top
            if (MINROW[i]>2) Temp[MINROW[i], MINCOL[i]:MAXCOL[i]] =
                    2*Temp[MINROW[i]-1, MINCOL[i]:MAXCOL[i]] -
                      Temp[MINROW[i]-2, MINCOL[i]:MAXCOL[i]]
            if (MAXROW[i]<NROW-1) Temp[MAXROW[i], MINCOL[i]:MAXCOL[i]] =
                    2*Temp[MAXROW[i]+1, MINCOL[i]:MAXCOL[i]] -
                      Temp[MAXROW[i]+2, MINCOL[i]:MAXCOL[i]]
            
            # Transfer heat flux along left and right
            # (left/right prevails over bottom/top)
            if (MINCOL[i]>2) Temp[MINROW[i]:MAXROW[i], MINCOL[i]] =
                    2*Temp[MINROW[i]:MAXROW[i], MINCOL[i]-1] -
                      Temp[MINROW[i]:MAXROW[i], MINCOL[i]-2]
            if (MAXCOL[i]<NCOL-1) Temp[MINROW[i]:MAXROW[i], MAXCOL[i]] =
                    2*Temp[MINROW[i]:MAXROW[i], MAXCOL[i]+1] -
                      Temp[MINROW[i]:MAXROW[i], MAXCOL[i]+2]
        }
    }
}


# Plot mean(T) evolution for every object
for (i in 1:NMAT) {
    png(paste0("tprofile_",heatobjparams$desc[i],".png"),width=512, height=400)
    plot(seq(0,dt*N/60,length.out=ncol(T_Evol)), T_Evol[i,],
         type='l', col='red',
         xlab='Time (min)', ylab='T (�C/K)', main=heatobjparams$desc[i],
         ylim=c(min(T_Evol[i,]),max(T_Evol[i,])))
    dev.off()   
}



# HEAT FLUX
MAINMAT=1
Temp[-lst[[MAINMAT]]]=NaN  # Evaluate q only on MAINMAT

# Calculate heat flux q=-k*grad(T): (NROW-2)x(NCOL-2) matrix
qx=-k[2:(NROW-1),2:(NCOL-1)] *
    (Temp[2:(NROW-1),3:NCOL] - Temp[2:(NROW-1),1:(NCOL-2)])/(2*dx)
qy=-k[2:(NROW-1),2:(NCOL-1)] *
    (Temp[1:(NROW-2),2:(NCOL-1)] - Temp[3:NROW,2:(NCOL-1)])/(2*dx)


# 3D plot of q modulus
qmod=(qx^2+qy^2)^0.5
z=qmod
color=magma(256)
persp3d(seq(0,1,length.out=nrow(z)), seq(0,1,length.out=ncol(z)), z,
        col=color[cut(z, 256)], axes=TRUE, box=TRUE,
        aspect=c(1,1*NCOL/NROW,0.5*NCOL/NROW))



# Plot T countour map
png("flux_simulation.png", width=640, height = 598)
tempe=t(Temp)[,nrow(Temp):1]
x=10*(1:nrow(tempe))
y=10*(1:ncol(tempe))
image(x, y, tempe, col = magma(256), axes = FALSE, asp=1)
contour(x, y, tempe, levels = seq(-70, 70, by = 5), labcex = 1.2,
        add = TRUE, col = "brown")
dev.off()



# FOLLOWING CODE NEEDS REVISION

# Plot T vector field

library(rasterVis)
library(raster)

# Plot 2D flux vectors
# flux.png is the last simulation iteration
img=readPNG("flux.png")
img[img==1]=NaN
image(img, col = magma(256), asp=1)
proj <- CRS('+proj=longlat +datum=WGS84')
df <- expand.grid(x = seq(-1, 1, length.out=nrow(img)),
                  y = seq(-222/200, -222/200, length.out=ncol(img)))
dim(img)=c(length(img), 1)
df$z <- as.data.frame(img)
r <- rasterFromXYZ(df, crs=proj)  # ERROR ???
pdf("fluxvector.pdf", lwd.arrows=3)
vectorplot(r, par.settings=RdBuTheme())
dev.off()

pdf("streamplot.pdf")
streamplot(r, maxpixels=300, aspX=0.02, aspY=0.02)
dev.off()


# THIS CODE WORKS
proj <- CRS('+proj=longlat +datum=WGS84')
df <- expand.grid(x = seq(-2, 2, .01), y = seq(-2, 2, .01))
df$z <- with(df, (3*x^2 + y)*exp(-x^2-y^2))
r <- rasterFromXYZ(df, crs=proj)
vectorplot(r, par.settings=RdBuTheme())  # plot vectorplot (arrows)
streamplot(r)  # plot stream

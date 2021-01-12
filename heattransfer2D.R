# Heat Transfer Equation (2D Finite-Difference Time-Domain approximation)
# www.overfitting.net

# iPad 9309

library(png)
library(stringr)
setwd("D:/R/43_HeatTransfer")


# Read simulation parameters
heatparams=read.table("heatparams.csv", header=T, sep=",")  # Excel CSV
dx=heatparams$value[heatparams$desc=='cell size (m)']
dt=heatparams$value[heatparams$desc=='time step (s)']
N=as.integer(heatparams$value[heatparams$desc=='number of iterations'])

# Read objects
heatobjobjparams=read.table("heatobjects.csv", header=T, sep=",")  # Excel CSV
heatobjects=readPNG("heatobjects.png")

plot(as.raster(heatobjects), interpolate=F)  # display objects
heatobjectsunique=heatobjects[,,3]+heatobjects[,,1]*2+heatobjects[,,2]*4
writePNG(heatobjectsunique/max(heatobjectsunique), "heatobjects_gray.png")

# Stability condition (if UNSTABLE reduce dt or increase dx)
# Threshold is 1/2 for 1D, 1/4 for 2D, 1/6 for 3D
r=dt*max(heatobjobjparams$alpha)/dx^2
print(paste0("r=", r, " -> ", ifelse(r<=1/4, "STABLE", "UNSTABLE"),
             " (",round(r/(1/4)*100),"%)"))


# Obtain unique colours and round to object position
INITMARK=heatobjectsunique[1,1]
i=2
while (heatobjectsunique[1,i] != INITMARK) i=i+1
NOBJECTS=i-2  # number of objects (colours) defined
print(paste0(NOBJECTS, " objects defined"))

colours=heatobjectsunique[1,2:(NOBJECTS+1)]  # keep just colours
heatobjectsunique=heatobjectsunique[-1,]  # drop first row (colour definition)
NROW=nrow(heatobjectsunique)
NCOL=ncol(heatobjectsunique)

lst=list()
for (i in 1:NOBJECTS) {
    indices=which(heatobjectsunique==colours[i])
    heatobjectsunique[indices]=i
    lst[[i]]=indices  # create list of objects
}
plot(as.raster(heatobjectsunique/NOBJECTS), interpolate=F)

# Plot scene in lines
heatobjectslines=heatobjectsunique*0
heatobjectslines[2:(NROW-1),2:(NCOL-1)]=
    abs(heatobjectsunique[3:NROW,2:(NCOL-1)] - heatobjectsunique[2:(NROW-1),2:(NCOL-1)])+
    abs(heatobjectsunique[2:(NROW-1),3:NCOL] - heatobjectsunique[2:(NROW-1),2:(NCOL-1)])+
    abs(heatobjectsunique[1:(NROW-2),2:(NCOL-1)] - heatobjectsunique[2:(NROW-1),2:(NCOL-1)])*0+
    abs(heatobjectsunique[2:(NROW-1),1:(NCOL-2)] - heatobjectsunique[2:(NROW-1),2:(NCOL-1)])*0
heatobjectslines[heatobjectslines !=0 ]=1

writePNG(heatobjectslines, "heatobjects_lines.png")

# Define working arrays
# alpha = k/(rho*cp), rho*cp=k/alpha
T=heatobjectsunique*0
k=T
alpha=T
for (i in 1:NOBJECTS) {
    T[lst[[i]]]=heatobjobjparams$T[i]
    k[lst[[i]]]=heatobjobjparams$k[i]
    alpha[lst[[i]]]=heatobjobjparams$alpha[i]
}
rhocp=k/alpha  # won't use alpha, just k and rhocp=rho * cp


# Time domain T iteration (vector notation)
MINT=min(heatobjobjparams$T)
MAXT=max(heatobjobjparams$T)

NSNAPSHOTS=200
SKIP=round(N/NSNAPSHOTS)
for (j in 0:N) {
    if (j %% SKIP==0) {
        nombre=paste0("heattransfer_",
                      str_pad(j, nchar(N), pad='0'), ".png")
        writePNG((T-MINT)/(MAXT-MINT), target=nombre)
    }

    # Formulation valid for multilayer walls
    T[2:(NROW-1),2:(NCOL-1)] = T[2:(NROW-1),2:(NCOL-1)] +
        dt/(rhocp[2:(NROW-1),2:(NCOL-1)] * dx^2) *
        (
            (k[3:NROW    ,2:(NCOL-1)] + k[2:(NROW-1),2:(NCOL-1)])/2 * (T[3:NROW    ,2:(NCOL-1)] - T[2:(NROW-1),2:(NCOL-1)]) -
            (k[1:(NROW-2),2:(NCOL-1)] + k[2:(NROW-1),2:(NCOL-1)])/2 * (T[2:(NROW-1),2:(NCOL-1)] - T[1:(NROW-2),2:(NCOL-1)]) +
                
            (k[2:(NROW-1),3:NCOL    ] + k[2:(NROW-1),2:(NCOL-1)])/2 * (T[2:(NROW-1),3:NCOL    ] - T[2:(NROW-1),2:(NCOL-1)]) -
            (k[2:(NROW-1),1:(NCOL-2)] + k[2:(NROW-1),2:(NCOL-1)])/2 * (T[2:(NROW-1),2:(NCOL-1)] - T[2:(NROW-1),1:(NCOL-2)])             
        )
    
    # Reset T in the boundaries and apply instantaneous convection on fluids
    for (i in 1:NOBJECTS) {
        if (heatobjobjparams$type[i]=='boundary') {
            T[lst[[i]]]=heatobjobjparams$T[i]  # boundary T conditions
        } else if (heatobjobjparams$type[i]=='fluid') {
            T[lst[[i]]]=mean(T[lst[[i]]])  # average T = quick convection           
        }
    }
}



# Algoritmo:
# NOTA: es obligatorio que todo el perímetro sea 'boundary'
#       como aprox. el perímetro también puede ser 'fluid'
# 1. Iterar T en todos los puntos salvo periferia (vectorización)
# 2. Sobreescribir T en los objetos 'boundary' (which() precalculados)
# 3. Promediar T en los objetos 'fluid' (which() precalculados)
#
# La periferia debería haberse reescrito en 2 o 3 (preferiblemente 2)

# Mejoras a introducir:
# 1. Fuentes de calor internas
# 2. Convección en fluidos
# 3. Cálculos de flujo de calor (vectorial)








NROW=60
NCOL=80
a=array(1:48,c(NROW,NCOL))

# Method 1: which() indexing
indices=which(col(a)>1 & col(a)<ncol(a) & row(a)>1 & row(a)<NROW)
# a[indices+1] baja una fila
# a[indices-1] sube una fila
# a[indices+NROW] avanza una columna
# a[indices-NROW] retrocede una columna
start_time=Sys.time()
for (i in 1:100000) a[indices]=a[indices]+
    (a[indices+1] + a[indices-1] +
     a[indices+NROW] + a[indices-NROW])/1000
end_time=Sys.time()
print(end_time - start_time)

# Method 2: array indexing
start_time=Sys.time()
for (i in 1:10000) a[2:(NROW-1),2:(NCOL-1)]=a[2:(NROW-1),2:(NCOL-1)]+
    (a[3:NROW,2:(NCOL-1)] + a[1:(NROW-2),2:(NCOL-1)] +
     a[2:(NROW-1),3:NCOL] + a[2:(NROW-1),1:(NCOL-2)])/1000
end_time=Sys.time()
print(end_time - start_time)



indices: 21.06793 secs
matricial: 5.608142 secs -> ELEGIMOS ESTA FORMA
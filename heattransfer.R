# Heat Transfer Equation (FDM approximation)
# www.overfitting.net

library(png)

# 1D HEAT TRANSFER
n=25  # number of cells
N=500  # number of iterations
L=0.1  # 0.1 m wall thickness
dx=L/n  # 0.01 m
dt=0.05  # 0.1 s

alpha=0.0001  # 0.0001 m2/s
T0=-5  # initial T in the wall
Tleft=-5  # left end boundary T
Tright=20  # right end boundary T

T=matrix(T0, nrow=1, ncol=n)
T[1]=Tleft
T[n]=Tright
indices=which(col(T)!=1 & col(T)!=n)
x=seq(from=dx/2, to=L-dx/2, length.out=n)  # position array

val=alpha*dt/dx^2
for (j in 0:N) {
    nombre=paste0("heat_transfer_",
        ifelse(j<10, "000", ifelse(j<100, "00", ifelse(j<1000, "0", ""))), j,
        ".png")
    png(nombre)
    plot(x, T, type='l', col='red', lwd=2,
         xlim=c(0,L), ylim=c(min(Tleft,Tright), max(Tleft,Tright)),
         ylab='T (ºC)',
         main=paste0('Iteration: ', j, '/',N, ', t=', round(j*dt), 's'))
    abline(v=c(x[1], x[n]), lty=2)
    abline(h=0)
    dev.off()
    
    # Iterate T in vector notation
    T[indices]=T[indices]+val*(T[indices+1] - 2*T[indices] + T[indices-1])
}

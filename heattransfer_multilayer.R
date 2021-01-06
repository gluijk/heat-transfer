# Heat Transfer Equation (1D Finite-Difference Time-Domain approximation)
# www.overfitting.net

library(png)
library(stringr)


# Parameters
N=500000L  # number of iterations
n=100L  # number of cells
L=0.2  # wall thickness (m)
dx=L/n  # cell size (m)
dt=0.005*3  # time step (s)

alpha0=0.0001  # layer A
alpha1=0.000001  # layer B
alpha2=0.00001  # layer C
T0=-5  # initial T in the wall
Tleft=-5  # left end boundary condition for T
Tright=20  # right end boundary condition for T


# Build vectors
T=matrix(T0, nrow=1, ncol=n)
T[1]=Tleft
T[n]=Tright
indices=which(col(T)!=1 & col(T)!=n)

SPLIT1=33  # layer A | layer B boundary
SPLIT2=66  # layer B | layer C boundary
alpha=matrix(alpha0, nrow=1, ncol=n)
alpha[SPLIT1:SPLIT2]=alpha1
alpha[SPLIT2:n]=alpha2

x=seq(from=dx/2, to=L-dx/2, length.out=n)  # position array


# Stability condition (if UNSTABLE reduce dt or increase dx)
r=dt*max(alpha)/dx^2
print(paste0("r=", r, " -> ", ifelse(r<=0.5, "STABLE", "UNSTABLE"),
             " (",round(r/0.5*100),"%)"))


# Time domain T iteration (vector notation)
val=alpha*dt/dx^2
for (j in 0:N) {
    if (j %% 1000 ==0) {
        nombre=paste0("heattransfer_",
                      str_pad(j, nchar(N), pad='0'), ".png")
        png(nombre)
        plot(x, T, type='l', col='red', lwd=2,
             xlim=c(0,L), ylim=c(min(Tleft,Tright), max(Tleft,Tright)),
             ylab='T (ºC)',
             main=paste0('Iteration: ', j, '/',N, ', t=', round(j*dt), 's'))
        abline(v=c(x[1], x[n],
                (x[SPLIT1-1]+x[SPLIT1])/2,
                (x[SPLIT2-1]+x[SPLIT2])/2), lty=2)
        abline(h=0)
        dev.off()
    }

    T[indices]=T[indices]+
        val[indices]*(T[indices+1] - 2*T[indices] + T[indices-1])  # delta T
}

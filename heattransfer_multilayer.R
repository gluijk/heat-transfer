# Heat Transfer Equation (1D Finite-Difference Time-Domain approximation)
# www.overfitting.net

library(png)
library(stringr)


# Parameters
N=100000L  # number of iterations
n=100L  # number of cells (2 are the outer ends)
n=101L  # number of cells (2 are the outer ends)

L=0.2  # wall thickness (m)
dx=L/n  # cell size (m)
dt=4  # time step (s)

alpha0=9.19e-8  # layer A: madera de pino
alpha1=0.23e-6  # layer B: lana de vidrio
alpha2=5.19e-8  # layer C: yeso (pladur)

T0=-10  # initial T in the wall
Tleft=-10  # left end boundary condition for T
Tright=20  # right end boundary condition for T


# Build vectors
T=matrix(T0, nrow=1, ncol=n)
T[1]=Tleft
T[n]=Tright
indices=which(col(T)!=1 & col(T)!=n)

SPLIT1=80  # first cell of layer B
SPLIT2=94  # first cell of layer C

alpha=matrix(alpha0, nrow=1, ncol=n)
alpha[SPLIT1:(SPLIT2-1)]=alpha1
alpha[SPLIT2:n]=alpha2

x=seq(from=dx/2, to=L-dx/2, length.out=n)  # position array


# Stability condition (if UNSTABLE reduce dt or increase dx)
r=dt*max(alpha)/dx^2
print(paste0("r=", r, " -> ", ifelse(r<=0.5, "STABLE", "UNSTABLE"),
             " (",round(r/0.5*100),"%)"))


# Time domain T iteration (vector notation)
val=alpha*dt/dx^2  # precalculated vector
for (j in 0:N) {
    if (j %% 1000==0) {
        nombre=paste0("heattransfer_",
                      str_pad(j, nchar(N), pad='0'), ".png")
        png(nombre)
        plot(x, T, type='l', col='red', lwd=2,
             xlim=c(0,L), ylim=c(min(Tleft,Tright), max(Tleft,Tright)),
             ylab='T (ºC)',
             main=paste0('Iteration: ', j, '/',N,
                         ', t=', as.integer(j*dt), 's'))
        abline(v=c((x[1]+x[2])/2, (x[n-1]+x[n])/2,
                (x[SPLIT1-1]+x[SPLIT1])/2,
                (x[SPLIT2-1]+x[SPLIT2])/2), lty=2)
        abline(h=0)
        dev.off()
    }

    # Formulation valid for a single layer
    # T[indices]=T[indices]+
    #     val[indices]*(T[indices+1] - 2*T[indices] + T[indices-1])  # delta T

    # Formulation valid for multilayer walls
    T[indices] = T[indices] +
        (val[indices]+val[indices+1])/2*(T[indices+1] - T[indices]) -
        (val[indices-1]+val[indices])/2*(T[indices] - T[indices-1])
    
    # Formulation suggested by Carlos Gil Bellosta
    # T_{t+1)i = T_ti + dt / dx^2 * (alpha[i+1] (T_t(i+1) - T_ti) - alpha[i-1] (T_ti - T_t(i-1)))
    #   T[indices] = T[indices] +
    #       val[indices+1]*(T[indices+1] - T[indices]) -
    #       val[indices-1]*(T[indices] - T[indices-1])

}

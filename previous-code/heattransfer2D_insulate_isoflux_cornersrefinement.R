# Refinement of corners on 'insulate' and 'isoflux' materials
# is avoided since it can lead to unstability

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
        if (MINCOL[i]>1) Temp[MINROW[i]:MAXROW[i], MINCOL[i]] =
                Temp[MINROW[i]:MAXROW[i], MINCOL[i]-1]
        if (MAXCOL[i]<NCOL) Temp[MINROW[i]:MAXROW[i], MAXCOL[i]] =
                Temp[MINROW[i]:MAXROW[i], MAXCOL[i]+1]
        
        # Refine T on 4 corners:
        # Bottom-Left
        if (MINROW[i]>1 & MINCOL[i]>1) {
            Temp[MINROW[i], MINCOL[i]] =
              (Temp[MINROW[i]-1, MINCOL[i]]+Temp[MINROW[i], MINCOL[i]-1])/2
        } else if (MINROW[i]>1) {  # MINCOL[i]=1
            Temp[MINROW[i], MINCOL[i]] = Temp[MINROW[i]-1, MINCOL[i]]              
        } else if (MINCOL[i]>1) {  # MINROW[i]=1
            Temp[MINROW[i], MINCOL[i]] = Temp[MINROW[i], MINCOL[i]-1]                 
        }
        
        # Top-Left
        if (MAXROW[i]<NROW & MINCOL[i]>1) {
            Temp[MAXROW[i], MINCOL[i]] =
              (Temp[MAXROW[i]+1, MINCOL[i]]+Temp[MAXROW[i], MINCOL[i]-1])/2
        } else if (MAXROW[i]<NROW) {  # MINCOL[i]=1
            Temp[MAXROW[i], MINCOL[i]] = Temp[MAXROW[i]+1, MINCOL[i]]        
        } else if (MINCOL[i]>1) {  # MAXROW[i]=NROW
            Temp[MAXROW[i], MINCOL[i]] = Temp[MAXROW[i], MINCOL[i]-1]               
        }
        
        # Bottom-Right
        if (MINROW[i]>1 & MAXCOL[i]<NCOL) {
            Temp[MINROW[i], MAXCOL[i]] =
              (Temp[MINROW[i]-1, MAXCOL[i]]+Temp[MINROW[i], MAXCOL[i]+1])/2
        } else if (MINROW[i]>1) {  # MAXCOL[i]=NCOL
            Temp[MINROW[i], MAXCOL[i]] = Temp[MINROW[i]-1, MAXCOL[i]]              
        } else if (MAXCOL[i]<NCOL) {  # MINROW[i]=1
            Temp[MINROW[i], MAXCOL[i]] = Temp[MINROW[i], MAXCOL[i]+1]                 
        }

        # Top-Right
        if (MAXROW[i]<NROW & MAXCOL[i]<NCOL) {
            Temp[MAXROW[i], MAXCOL[i]] =
              (Temp[MAXROW[i]+1, MAXCOL[i]]+Temp[MAXROW[i], MAXCOL[i]+1])/2
        } else if (MAXROW[i]<NROW) {  # MAXCOL[i]=NCOL
            Temp[MAXROW[i], MAXCOL[i]] = Temp[MAXROW[i]+1, MAXCOL[i]]              
        } else if (MAXCOL[i]<NCOL) {  # MAXROW[i]=NROW
            Temp[MAXROW[i], MAXCOL[i]] = Temp[MAXROW[i], MAXCOL[i]+1]                 
        }
        
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
        if (MINCOL[i]>2) Temp[MINROW[i]:MAXROW[i], MINCOL[i]] =
                2*Temp[MINROW[i]:MAXROW[i], MINCOL[i]-1] -
                  Temp[MINROW[i]:MAXROW[i], MINCOL[i]-2]
        if (MAXCOL[i]<NCOL-1) Temp[MINROW[i]:MAXROW[i], MAXCOL[i]] =
                2*Temp[MINROW[i]:MAXROW[i], MAXCOL[i]+1] -
                  Temp[MINROW[i]:MAXROW[i], MAXCOL[i]+2]
        
        # Refine heat flux on 4 corners:
        # Bottom-Left
        if (MINROW[i]>2 & MINCOL[i]>2) {
            Temp[MINROW[i], MINCOL[i]] =
              (Temp[MINROW[i]-1, MINCOL[i]]+Temp[MINROW[i], MINCOL[i]-1]) -
              (Temp[MINROW[i]-2, MINCOL[i]]+Temp[MINROW[i], MINCOL[i]-2])/2  
        } else if (MINROW[i]>2) {  # MINCOL[i]=1,2
            Temp[MINROW[i], MINCOL[i]] =
                2*Temp[MINROW[i]-1, MINCOL[i]] -
                  Temp[MINROW[i]-2, MINCOL[i]]              
        } else if (MINCOL[i]>2) {  # MINROW[i]=1,2
            Temp[MINROW[i], MINCOL[i]] =
                2*Temp[MINROW[i], MINCOL[i]-1] -
                  Temp[MINROW[i], MINCOL[i]-2]             
        }
        
        # Top-Left
        if (MAXROW[i]<NROW-1 & MINCOL[i]>2) {
            Temp[MAXROW[i], MINCOL[i]] =
              (Temp[MAXROW[i]+1, MINCOL[i]]+Temp[MAXROW[i], MINCOL[i]-1]) -
              (Temp[MAXROW[i]+2, MINCOL[i]]+Temp[MAXROW[i], MINCOL[i]-2])/2 
        } else if (MAXROW[i]<NROW-1) {  # MINCOL[i]=1,2
            Temp[MAXROW[i], MINCOL[i]] =
                2*Temp[MAXROW[i]+1, MINCOL[i]] -
                  Temp[MAXROW[i]+2, MINCOL[i]]
        } else if (MINCOL[i]>2) {  # MAXROW[i]=NROW
            Temp[MAXROW[i], MINCOL[i]] =
                2*Temp[MAXROW[i], MINCOL[i]-1] -
                  Temp[MAXROW[i], MINCOL[i]-2]
        }
        
        # Bottom-Right
        if (MINROW[i]>2 & MAXCOL[i]<NCOL-1) {
            Temp[MINROW[i], MAXCOL[i]] =
              (Temp[MINROW[i]-1, MAXCOL[i]]+Temp[MINROW[i], MAXCOL[i]+1]) -
              (Temp[MINROW[i]-2, MAXCOL[i]]+Temp[MINROW[i], MAXCOL[i]+2])/2
        } else if (MINROW[i]>2) {  # MAXCOL[i]=NCOL,NCOL-1
            Temp[MINROW[i], MAXCOL[i]] =
                2*Temp[MINROW[i]-1, MAXCOL[i]] -
                  Temp[MINROW[i]-2, MAXCOL[i]]
        } else if (MAXCOL[i]<NCOL-1) {  # MINROW[i]=1,2
            Temp[MINROW[i], MAXCOL[i]] =
                2*Temp[MINROW[i], MAXCOL[i]+1] -
                  Temp[MINROW[i], MAXCOL[i]+2]
        }
        
        # Top-Right
        if (MAXROW[i]<NROW-1 & MAXCOL[i]<NCOL-1) {
            Temp[MAXROW[i], MAXCOL[i]] =
              (Temp[MAXROW[i]+1, MAXCOL[i]]+Temp[MAXROW[i], MAXCOL[i]+1]) -
              (Temp[MAXROW[i]+2, MAXCOL[i]]+Temp[MAXROW[i], MAXCOL[i]+2])/2
        } else if (MAXROW[i]<NROW-1) {  # MAXCOL[i]=NCOL,NCOL-1
            Temp[MAXROW[i], MAXCOL[i]] =
                2*Temp[MAXROW[i]+1, MAXCOL[i]] -
                  Temp[MAXROW[i]+2, MAXCOL[i]]
        } else if (MAXCOL[i]<NCOL-1) {  # MAXROW[i]=NROW,NROW-1
            Temp[MAXROW[i], MAXCOL[i]] =
                2*Temp[MAXROW[i], MAXCOL[i]+1] -
                  Temp[MAXROW[i], MAXCOL[i]+2]
        }        
    }
}


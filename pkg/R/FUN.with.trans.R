FUN.with.trans <- function(z, N, T, is.intercept,
                           effect = c("none", "individual", "time", "twoways")) 
  {
    with.trans <- match.arg(effect)
    if(!is.intercept & mean(z)!=0){
      stop("You want to estimate a model without intercept,
           but the arithmetic mean of one or more of your variables is not zero.")
    }
    
    mean.z     <- ifelse(is.intercept, mean(z), 0);

    switch(with.trans	
           , none = {
             if(is.intercept){Z  = z - mean.z}
             else{            Z  = z};
             liste <- list("Tr"  = "none", # Name of *Tr*ansformation
                           "I"   = ifelse(is.intercept, TRUE, FALSE), 
                           "ODM" = z,        # *O*rig. *D*ata *M*atrix
                           "TDM" = Z,        # *T*ransformed *D*ata *M*atrix
                           "TDV" = c(Z),     # *T*ransformed *D*ata *V*ector
                           "OVm" = mean.z,   # *OV*erall *m*ean
                           "TRm" = 0);       # *TR*ansformation *m*ean. 
                           liste
           }
           , individual = {
             if(is.intercept){Z  = z - matrix(colMeans(z), T, N, byrow = TRUE)- mean.z}
             else{            Z  = z - matrix(colMeans(z), T, N, byrow = TRUE)};
             liste <- list("Tr"  = "individual", # Name of *Tr*ansformation
                           "I"   = ifelse(is.intercept, TRUE, FALSE), 
                           "ODM" = z,            # *O*rig. *D*ata *M*atrix
                           "TDM" = Z,            # *T*ransformed *D*ata *M*atrix
                           "TDV" = c(Z),         # *T*ransformed *D*ata *V*ector
                           "OVm" = mean.z,       # *OV*erall *m*ean
                           "TRm" = colMeans(z)); # *TR*ansformation *m*ean
             liste								
                         }
           , time = {
             if(is.intercept){Z  = z - rowMeans(z) - mean.z}
             else{            Z  = z - rowMeans(z)};
             liste <- list("Tr"  = "time",       # Name of *Tr*ansformation
                           "I"   = ifelse(is.intercept, TRUE, FALSE), 
                           "ODM" = z,            # *O*rig. *D*ata *M*atrix
                           "TDM" = Z,            # *T*ransformed *D*ata *M*atrix
                           "TDV" = c(Z),         # *T*ransformed *D*ata *V*ector
                           "OVm" = mean.z,       # *OV*erall *m*ean
                           "TRm" = rowMeans(z)); # *TR*ansformation *m*ean
             liste                                  
           }
           , twoways = {
             if(is.intercept){Z  = z - matrix(colMeans(z), T, N, byrow = TRUE) - rowMeans(z) + mean.z}
             else{            Z  = z - matrix(colMeans(z), T, N, byrow = TRUE) - rowMeans(z)};
             liste <- list("Tr"  = "twoway",   # Name of *Tr*ansformation
                           "I"   = ifelse(is.intercept, TRUE, FALSE), 
                           "ODM" = z,          # *O*rig. *D*ata *M*atrix
                           "TDM" = Z,          # *T*ransformed *D*ata *M*atrix
                           "TDV" = c(Z),       # *T*ransformed *D*ata *V*ector
                           "OVm" = mean.z,     # *OV*erall *M*ean
                           "TRm" = list("individual" = colMeans(z),
                                        "time"       = rowMeans(z)));# *TR*ansformation *m*eans
             liste
           }
           )		 		
  }

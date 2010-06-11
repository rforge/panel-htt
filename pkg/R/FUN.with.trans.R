FUN.with.trans <- function(z, N, T, is.intercept, effect = c("pooled", "individual", "time", "twoways")) #, "none"
		{
                  with.trans <- match.arg(effect)
                  mean.z <- mean(z);
                  switch(with.trans	
                         , pooled = {
                           if(is.intercept){Z  = z - mean.z}
                           else{            Z  = z};
                           liste <- list("T"   = "pooled", # Name of *T*ransformation
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
                           liste <- list("T"   = "individual", # Name of *T*ransformation
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
                           liste <- list("T"   = "time",       # Name of *T*ransformation
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
                           liste <- list("T"   = "twoway",   # Name of *T*ransformation
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

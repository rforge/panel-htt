## Änderungen:
## 1: eine cover-Funktion geschrieben 'WithTrans'
## 2: logische Umfrage 'if(!is.intercept & mean(z)!=0)' eliminiert "intercep kann in einem 
##    dummy design matrix eingegeben. Sonst werden die Konstante Effecte in der Factorstructure 
##    erscheinen" 
## 3: OVm ist jetzt unter TRm
## 4: Trm ist immer eine Liste für alle Effekte
## 5: für 'Twoways' ein Intercept muss immer im model sein sonst braucht man zusätiche
##    Restriktionen für die Identifikation.
##    Also für 'Twoways' das argument 'is.intercept will be ingonred'

## To Do: Evtl sollte man Einige listen-Componenten mit "drop=FALSE" ausstatten,
## da z.B. ODM bei P=1 ein vector ist und keine matrix! 


FUN.with.trans <- function(z, N, T, is.intercept, effect = c("none", "individual", "time", "twoways")){
	with.trans <- match.arg(effect)
        const <- ifelse(is.intercept, mean(z), 0)
        
	switch(with.trans	
           , none = {
			Z  	<- z - const
             	liste <- list(
			   "Tr"  = "none", # Name of *Tr*ansformation
                           "I"   = ifelse(is.intercept, TRUE, FALSE), 
                           "ODM" = z,        # *O*rig. *D*ata *M*atrix
                           "TDM" = Z,        # *T*ransformed *D*ata *M*atrix
                           "TDV" = c(Z),     # *T*ransformed *D*ata *V*ector
                           "TRm" = list(     # *TR*ansformation *m*eans. 
						"OVc"  = const,	# *OV*erall *c*onstant
						"InC"  = rep.int(const, N),	# *In*dividual *C*onstants
						"TiVC" = rep.int(const, T))	# *Ti*m V*arying *C*onstants
				   )    
                  return(liste)
           }
           , individual = {
			InC	<- colMeans(z)
			Z  	<- z - matrix(InC, T, N, byrow = TRUE) 
             	liste <- list(
			   "Tr"  = "individual", # Name of *Tr*ansformation
                           "I"   = ifelse(is.intercept, TRUE, FALSE), 
                           "ODM" = z,            # *O*rig. *D*ata *M*atrix
                           "TDM" = Z,            # *T*ransformed *D*ata *M*atrix
                           "TDV" = c(Z),         # *T*ransformed *D*ata *V*ector
                           "TRm" = list(	  	 # *TR*ansformation *m*eans. 
						"OVc"  = const,	# *OV*erall *c*onstant
						"InC"  = InC,	# *In*dividual *C*onstants
						"TiVC" = rep.int(const, T))	# *Ti*m V*arying *C*onstants
				   )  
                  return(liste)								
           }
           , time = {
			TiVC <- rowMeans(z)
	           	Z  	<- z - TiVC
             	liste <- list(
				   "Tr"  = "time",       # Name of *Tr*ansformation
                           "I"   = ifelse(is.intercept, TRUE, FALSE), 
                           "ODM" = z,            # *O*rig. *D*ata *M*atrix
                           "TDM" = Z,            # *T*ransformed *D*ata *M*atrix
                           "TDV" = c(Z),         # *T*ransformed *D*ata *V*ector
                           "TRm" = list(	  	 # *TR*ansformation *m*eans. 
						"OVc"  = const,	# *OV*erall *c*onstant
						"InC"  = rep.int(const, N), 	# *In*dividual *C*onstants
						"TiVC" = TiVC)	# *Ti*m V*arying *C*onstants
				   )  
                  return(liste)								
           }
           , twoways = {
			InC   <- colMeans(z)
			TiVC  <- rowMeans(z)
			Cons  <- mean(z)
             	Z  	<- z - matrix(InC, T, N, byrow = TRUE) - TiVC + Cons
             	liste <- list(
                              "Tr"  = "twoway",   # Name of *Tr*ansformation
                              "I"   = ifelse(is.intercept, TRUE, FALSE), 
                              "ODM" = z,          # *O*rig. *D*ata *M*atrix
                              "TDM" = Z,          # *T*ransformed *D*ata *M*atrix
                              "TDV" = c(Z),       # *T*ransformed *D*ata *V*ector
                              "TRm" = list(	  	 # *TR*ansformation *m*eans. 
						"OVc"  = Cons,	# *OV*erall *c*onstant
						"InC"  = InC, 	# *In*dividual *C*onstants
						"TiVC" = TiVC)	# *Ti*m V*arying *C*onstants
				   )  
                  return(liste)								
           }
    )		 		
  }

## cover =================================================================================================
WithTrans <- function(z, intercept= TRUE, heto.effect = c("none", "individual", "time", "twoways"))
  {
	is.regular.panel(z, stopper = TRUE)
	nr <- nrow(z)
	nc <- ncol(z)
	with.trans <- match.arg(heto.effect)
	FUN.with.trans(z, N = nc, T = nr, is.intercept = intercept, effect = with.trans)
  }
##========================================================================================================


############################ test

## Y <- matrix(rnorm(100), 10 , 10)
## WithTrans(Y, intercept= T, heto.effect = "individual")

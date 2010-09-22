FUN.Pformula <- function(formula, effect)
  {
    data.fra <- model.frame(formula)
    dat.term <- attr(data.fra, "terms")	

    ## Construct data from formula
    ## dim(response)=TxN == dim(y.matrix)=TxN:
    y.matrix <- model.response(data.fra, "numeric")  

    ## 1)Extract Regressors
    ## 2)Check the presence of a 'intercept' in the formula
    ## 3)And built the x.all.matrix (TxN*P)

    regressors.mat <- model.matrix(dat.term, data.fra)
    is.intercept   <- ifelse(colnames(regressors.mat)[1] == "(Intercept)", TRUE, FALSE)
    if(is.intercept){
      x.all.matrix       <- regressors.mat[,-1]
    }else{x.all.matrix   <- regressors.mat}

    if(!is.intercept & effect=="twoways"){stop("Effects >> twoways << need an Intercept!")}
    
    ## Dimension parameters

    N  <- ncol(y.matrix)
    T  <- nrow(y.matrix)
    NT <- N*T
    P  <- as.integer(ncol(x.all.matrix)/N)

    ## Write the response, Y, and the 'p' regressors, X, in a list,
    ## where each component contains one of p+1 TxN-Matrices

    data.all.mat  <- cbind(y.matrix, x.all.matrix)
    
    ## New of Object: From Matrix data.all.mat a List model.in.list 

    model.in.list <- lapply(1:(P+1), function(z, i) z[,seq((i-1)*N+1,i*N)], z = data.all.mat)
    data.in.list  <- sapply(model.in.list,
                            function(z) FUN.with.trans(z,
                                                       N            = N,
                                                       T            = T,
                                                       is.intercept = is.intercept,
                                                       effect       = effect),
                            simplify=FALSE)
    data.in.list
  }

FUN.Pformula <- function(formula, effect)
	{
        data.fra <- model.frame(formula)
        dat.term <- attr(data.fra, "terms")	

# construct data from formula
	y.matrix <- model.response(data.fra, "numeric")  # dim(response)=TxN and dim(y.matrix)=TxN:

  # check the presence of a 'intercept' in the formula

	regressors.mat <- model.matrix(dat.term, data.fra)
	is.intercept   <- ifelse(colnames(regressors.mat)[1] == "(Intercept)", TRUE, FALSE)
	if(is.intercept){
          x.all.matrix <- regressors.mat[,-1]
        }else{x.all.matrix   <- regressors.mat}

  # dimention parameters

	N  <- ncol(y.matrix)
	T  <- nrow(y.matrix)
        NT <- N*T
	P  <- as.integer(ncol(x.all.matrix)/N)

  # write the 'P' regressors in a list where each component contains a T x N matix

	data.all.mat <- cbind(y.matrix, x.all.matrix)
	model.in.list <- lapply(1:(P+1), function(z, i) z[,seq((i-1)*N+1,i*N)], z = data.all.mat)
        data.in.list <- sapply(model.in.list,
                            function(z) FUN.with.trans(z
                                                     , N=N
                                                     , T=T
                                                     , is.intercept=is.intercept
                                                     , effect=effect), simplify=FALSE)


	data.in.list
	}

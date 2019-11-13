augment <- function(corm,t=0.001){

    n <- ncol(corm)
    names <- colnames(as.data.frame(corm))
    index <- as.character(1:n)
    colnames(corm) <- rownames(corm) <- index
    perm <- sample(1:n)
    corm <- (corm[perm,])[,perm]

    B <- as.data.frame(t(chol(corm)))

    B[n+1,] <- B[,n+1] <- 0
    B[n+1,1] <- runif(1,-1,1)
    for (j in 2:n){
        B_up <- B_lo <- as.matrix(B)
        B_up[n+1,j] <-  sqrt(1-sum(B_up[n+1,1:(j-1)]^2))
        B_lo[n+1,j] <- -sqrt(1-sum(B_lo[n+1,1:(j-1)]^2))

        upper <- tcrossprod(B_up)[n+1,j]
        lower <- tcrossprod(B_lo)[n+1,j]

        if (upper - lower <= t){ cor <- (upper + lower)/2
        }else{cor <- runif(1,lower,upper)}

        B[n+1,j] <- 1/B[j,j] * (cor - sum(B[n+1,1:(j-1)] * B[j,1:(j-1)]))
    }
    B[n+1,n+1] <- sqrt(1 - sum(B[(n+1),1:n]^2))

    out <- tcrossprod(as.matrix(B))
    out <- out[match(c(index,n+1),c(rownames(out),n+1)),match(c(index,n+1),c(colnames(out),n+1))]
    colnames(out)[1:n] <- rownames(out)[1:n] <- names

    out
}

augmentcpp <- function(corm,buff){

    n <- ncol(corm)
    names <- colnames(as.data.frame(corm))
    index <- as.character(1:n)
    colnames(corm) <- rownames(corm) <- index
    perm <- sample(1:n)
    corm <- (corm[perm,])[,perm]

    B <- as.data.frame(t(chol(corm)))

    B[n+1,] <- B[,n+1] <- 0
    B[n+1,1] <- runif(1,-1,1)
    B[n+1,n+1] <- sqrt(1-sum(B[n+1,]^2))
    B <- as.matrix(B)

    out <- matprod(augment_loop(B,n,buff))

    colnames(out) <- rownames(out) <- c(perm,n+1)
    out <- out[match(c(index,n+1),c(rownames(out),n+1)),match(c(index,n+1),c(colnames(out),n+1))]
    colnames(out)[1:n] <- rownames(out)[1:n] <- names

    out
}

factory <- function(fun){
  warn <- err <- NULL
  res <- withCallingHandlers(
    tryCatch(fun, error=function(e) {
      err <<- conditionMessage(e)
      NULL
    }), warning=function(w) {
      warn <<- append(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(res, warn=warn, err=err)
}

RandomCorm <- function(nvars){
	cor <- runif(1,-1,1)
	base <- matrix(c(1,cor,cor,1),nrow=2)

	while(ncol(base) < nvars){
		base <- augment(base)
		}
	base
}

RandomCormCPP <- function(nvars,buff=0.01){
	cor <- runif(1,-1,1)
	base <- matrix(c(1,cor,cor,1),nrow=2)

	if(nvars >= 3){
	base <- augmentcpp(base,buff)

	while(ncol(base) < nvars){
		base <- augmentcpp(base,buff)
		}
	}
	base
}

DOPE <- function(mod,nsims=10000,language="cpp",cl=NULL){
              output <- list()
              mod_mat <- model.frame(m)
              names <- c(colnames(mod_mat)[-1],"ControlFunction","R_Squared")
              vcvm <- cov(mod_mat)

              if(language == "cpp"){
                out <- as.data.frame(t(pbapply::pbsapply(1:nsims,function(x)simfuncpp(vcvm),cl=cl)))
              }

              if(language == "R"){
                out <- as.data.frame(t(pbapply::pbsapply(1:nsims,function(x)simfun(vcvm),cl=cl)))
              }
      colnames(out) <- names
      out
}

pctpos <- function(DOPE_OUTPUT){t(apply(DOPE_OUTPUT,2,function(x)(length(x[which(x>0)])/length(x))))}

simfuncpp <- function(vcvm){
                aug <- augmentcpp(vcvm,buff=.Machine$double.eps)
                zz <- aug[-1,-1]
                zy <- as.matrix(aug[1,2:ncol(aug)])
                betas <- solve(zz) %*% zy
                rsq <- (t(zy) %*% solve(zz) %*% zy)/aug[1,1]
                output <- c(t(betas),rsq)
                output
}

simfun <- function(vcvm){
                aug <- augment(vcvm,t=.Machine$double.eps)
                zz <- aug[-1,-1]
                zy <- as.matrix(aug[1,2:ncol(aug)])
                betas <- solve(zz) %*% zy
                rsq <- (t(zy) %*% solve(zz) %*% zy)/aug[1,1]
                output <- c(t(betas),rsq)
                output
}

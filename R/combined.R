augment <- function(corm,buff=sqrt(.Machine$double.eps)){

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

        if (upper - lower <= buff){ cor <- (upper + lower)/2
        }else{cor <- runif(1,lower,upper)}

        B[n+1,j] <- 1/B[j,j] * (cor - sum(B[n+1,1:(j-1)] * B[j,1:(j-1)]))
    }
    B[n+1,n+1] <- sqrt(1 - sum(B[(n+1),1:n]^2))

    out <- tcrossprod(as.matrix(B))
    out <- out[match(c(index,n+1),c(rownames(out),n+1)),match(c(index,n+1),c(colnames(out),n+1))]
    colnames(out)[1:n] <- rownames(out)[1:n] <- names
    colnames(out)[n+1] <- "Augment"
    
    out
}

augmentcpp <- function(corm,buff=sqrt(.Machine$double.eps)){

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

RandomCorm <- function(nvars,buff=sqrt(.Machine$double.eps)){
	cor <- runif(1,-1,1)
	base <- matrix(c(1,cor,cor,1),nrow=2)

	while(ncol(base) < nvars){
		base <- augment(base,buff)
	}
	rownames(base) <- colnames(base) <- paste0("V",1:nvars)
	base
}

RandomCormCPP <- function(nvars,buff=sqrt(.Machine$double.eps)){
	cor <- runif(1,-1,1)
	base <- matrix(c(1,cor,cor,1),nrow=2)

	if(nvars >= 3){
	base <- augmentcpp(base,buff)

	while(ncol(base) < nvars){
		base <- augmentcpp(base,buff)
		}
	}
	rownames(base) <- colnames(base) <- paste0("V",1:nvars)
	base
}

DOPE <- function(mod,nsims=10000,language="cpp",n.cores=1){
  output <- list()
  mm <- model.matrix(mod)
  mod_mat <- as.matrix(data.frame(y=model.frame(mod)[,1],mm[,-1]))
  names <- c(colnames(mm)[-1],"ControlFunction","R_Squared")
  vcvm <- cov(mod_mat)
  
  if(n.cores==1){
    cl <- NULL
  }else{
    cl <- parallel::makeCluster(n.cores)
    parallel::clusterEvalQ(cl,library(DOPE))
  }
  
  if(language == "cpp"){
    out <- as.data.frame(t(pbapply::pbsapply(1:nsims,function(x)simfuncpp(vcvm),cl=cl)))
    parallel::stopCluster(cl)
  }
  
  if(language == "R"){
    out <- as.data.frame(t(pbapply::pbsapply(1:nsims,function(x)simfun(vcvm),cl=cl)))
    parallel::stopCluster(cl)
  }
  colnames(out) <- names
  out
}

pctpos <- function(DOPE_OUTPUT){
    t(apply(DOPE_OUTPUT,2,function(x)(length(x[which(x>0)])/length(x))))
}

simfuncpp <- function(vcvm){
                psd <- FALSE
                attempts <- 0
                while(psd == FALSE){
                  attempts <- attempts + 1
                  if(attempts > 10){
                    print("Ill Conditioned System")
                    break
                  }
                  aug <- augmentcpp(vcvm,buff=0.01)
                  psd <- !any(eigen(aug)$values < 1e-8)
                }
                
                zz <- aug[-1,-1]
                zy <- as.matrix(aug[1,2:ncol(aug)])
                betas <- solve(zz) %*% zy
                rsq <- (t(zy) %*% solve(zz) %*% zy)/aug[1,1]
                output <- c(t(betas),rsq)
                output
}

simfun <- function(vcvm){
                psd <- FALSE
                attempts <- 0
                while(psd == FALSE){
                  attempts <- attempts + 1
                  if(attempts > 10){
                    print("Ill Conditioned System")
                    break
                  }
                  aug <- augment(vcvm,buff=0.01)
                  psd <- !any(eigen(aug)$values < 1e-8)
                }
                zz <- aug[-1,-1]
                zy <- as.matrix(aug[1,2:ncol(aug)])
                betas <- solve(zz) %*% zy
                rsq <- (t(zy) %*% solve(zz) %*% zy)/aug[1,1]
                output <- c(t(betas),rsq)
                output
}

plot_DOPE <- function(DOPE_OUTPUT,which_var=1,type="histogram",width=2,bindiv=2){
  
  tmp <- DOPE_OUTPUT[,which_var]
  pos <- pctpos(DOPE_OUTPUT)[which_var]
  d <- density(tmp)
  mode <- d$x[which.max(d$y)]
  
  if(type=="density"){
    plot(d,xlim=c(median(tmp)-width*sd(tmp),median(tmp)+width*sd(tmp)),
         main=paste0("DOPE: ", names(DOPE_OUTPUT)[which_var]),lwd=3)
    legend("topleft",legend = c(paste0("Percent Positive: ",round(pos,3)),
                                paste0("Median: ",round(median(tmp),3)),
                                #paste0("Mode: ", round(mode,3)),
                                paste0("Lower95: ",round(quantile(tmp,0.025),3)),
                                paste0("Upper95: ",round(quantile(tmp,0.975),3)),
                                paste0("Draws: ", length(tmp))
                                )
           )
  }
  
  if(type=="histogram"){
    hist(tmp,breaks=nrow(DOPE_OUTPUT)/bindiv,freq=F,main=paste0("DOPE: ", names(DOPE_OUTPUT)[which_var]),
                       xlab="Parameter Value",
                       xlim=c(median(tmp)-width*sd(tmp),median(tmp)+width*sd(tmp)))
    legend("topleft",legend = c(paste0("Percent Positive: ",round(pos,3)),
                                paste0("Median: ",round(median(tmp),3)),
                                #paste0("Mode: ", round(mode,3)),
                                paste0("Lower95: ",round(quantile(tmp,0.025),3)),
                                paste0("Upper95: ",round(quantile(tmp,0.975),3)),
                                paste0("Draws: ", length(tmp))
                                )
          )
  }

}


add_DOPE <- function(DOPE_OUTPUT,restriction,which_var=1,type="histogram",bindiv=2){

  sub <- with(DOPE_OUTPUT,eval(parse(text=restriction)))
  tmp <- DOPE_OUTPUT[sub,which_var]
  pos <- pctpos(DOPE_OUTPUT[sub,])[which_var]
  d <- density(tmp)
  mode <- d$x[which.max(d$y)]
  
  if(type=="density"){
    lines(d,col="blue",lwd=3)
    legend("topright",legend = c(paste0("Percent Positive: ",round(pos,3)),
                                paste0("Median: ",round(median(tmp),3)),
                                #paste0("Mode: ", round(mode,3)),
                                paste0("Lower95: ",round(quantile(tmp,0.025),3)),
                                paste0("Upper95: ",round(quantile(tmp,0.975),3)),
                                paste0("Draws: ", length(tmp))
    )
    )
  }
  
  if(type=="histogram"){
    hist(tmp,breaks=nrow(DOPE_OUTPUT)/bindiv,freq=F,add=T,col="red")
    legend("topright",legend = c(paste0("Percent Positive: ",round(pos,3)),
                                paste0("Median: ",round(median(tmp),3)),
                                #paste0("Mode: ", round(mode,3)),
                                paste0("Lower95: ",round(quantile(tmp,0.025),3)),
                                paste0("Upper95: ",round(quantile(tmp,0.975),3)),
                                paste0("Draws: ", length(tmp))
    )
    )
  }
  
}

mode_DOPE <- function(DOPE_OUTPUT,which_var=1){
  x <- DOPE_OUTPUT[,which_var]
  d <- density(x)
  d$x[which.max(d$y)]
}

noise_plot <- function(DOPE_OUTPUT, which_var = 1){
  
  DOPE_OUTPUT <- selectorate
  which_var <- 2
  
  rsqs <- seq(min(DOPE_OUTPUT$R_Squared)+0.0001,1,length=50)
  ppos <- unlist(lapply(rsqs,function(x)pctpos(DOPE_OUTPUT[which(DOPE_OUTPUT$R_Squared < x),])[which_var]))
  sgn <- sign(median(DOPE_OUTPUT[,which_var]))
  if(sgn < 0){
    ppos <- 1-ppos
  }
  infoloss <- unlist(lapply(rsqs,function(x)nrow(DOPE_OUTPUT[which(DOPE_OUTPUT$R_Squared < x),])))
  infoloss <- 1- infoloss/nrow(DOPE_OUTPUT)
  
  par(mar=c(5,5,2,5))
  plot(rsqs,ppos,xlab="Unmodeled Systematic Variation",ylab="Percent Results Same Sign as Naive",type="p",
       ylim=c(0.3,1))
  legend("bottomleft",legend=c("Coefficient Sign Certainty","Assumption Restrictiveness","Maximal Uncertainty")
         ,pch=c(1,19,NA),lty=c(0,0,2))
  mpps <- pctpos(DOPE_OUTPUT)[which_var]
  if(sgn < 0){
    mpps <- 1-mpps
    abline(h=mpps,lty=2,lwd=2)
  }else{
    abline(h=mpps,lty=2,lwd=2)
  }
  par(new=T)
  plot(rsqs,infoloss,axes=F,xlab=NA,ylab=NA,pch=19)
  axis(side=4)
  mtext(side = 4, line =3, "Percent DOPE Draws Rejected")
  
  #auc <- sum(diff(rsqs)[1]*ppos)/((1-min(rsqs))*(1-mpps))
  #legend("topright",legend=paste0("AUC: ",round(auc,4)),lty=0,pch=NA)
}

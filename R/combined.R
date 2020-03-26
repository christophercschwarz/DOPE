augment <- function(covm, buff = sqrt(.Machine$double.eps)){
  Is <- sqrt(1/diag(covm))
  corm <- diag(Is) %*% covm %*% t(diag(Is))
  
  n <- ncol(corm)
  names <- colnames(as.data.frame(corm))
  index <- as.character(1:n)
  colnames(corm) <- rownames(corm) <- index
  perm <- sample(1:n)
  corm <- (corm[perm, ])[, perm]
  
  B <- as.data.frame(t(chol(corm)))
  
  B[n + 1, ] <- B[, n + 1] <- 0
  B[n + 1, 1] <- runif(1, -1, 1)
  for (j in 2:n) {
    B_up <- B_lo <- as.matrix(B)
    B_up[n + 1, j] <- sqrt(1 - sum(B_up[n + 1, 1:(j - 1)]^2))
    B_lo[n + 1, j] <- -sqrt(1 - sum(B_lo[n + 1, 1:(j - 1)]^2))
    upper <- tcrossprod(B_up)[n + 1, j]
    lower <- tcrossprod(B_lo)[n + 1, j]
    if (upper - lower <= buff) {
      cor <- (upper + lower)/2
    }
    else {
      cor <- runif(1, lower, upper)
    }
    B[n + 1, j] <- 1/B[j, j] * (cor - sum(B[n + 1, 1:(j - 
                                                        1)] * B[j, 1:(j - 1)]))
  }
  B[n + 1, n + 1] <- sqrt(1 - sum(B[(n + 1), 1:n]^2))
  
  out <- tcrossprod(as.matrix(B))
  
  out <- out[match(c(index, n + 1), c(rownames(out), n + 1)), 
             match(c(index, n + 1), c(colnames(out), n + 1))]
  colnames(out)[1:n] <- rownames(out)[1:n] <- names
  colnames(out)[n + 1] <- "Augment"
  
  Is <- sqrt(1/Is^2)
  diag(Is) %*% out %*% t(diag(Is))
}

augmentcpp <- function(covm,buff=sqrt(.Machine$double.eps)){
    Is <- sqrt(1/diag(covm))
    corm <- diag(Is) %*% covm %*% t(diag(Is))
    
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

    Is <- sqrt(1/Is^2)
    diag(Is) %*% out %*% t(diag(Is))
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

stats <- function(coefs){
  draws <- length(coefs)
  med <- median(coefs)
  l95 <- quantile(coefs,0.025)
  u95 <- quantile(coefs,0.975)
  ppos <- sum(coefs>0)/draws
  out <- data.frame(Summary = c("Percent Positive",
                                "Median",
                                "Lower 95",
                                "Upper 95",
                                "Draws"),
                    Value = c(ppos*100,med,l95,u95,draws))
  out$Value <- round(out$Value,2)
  out
}


plot_DOPE <- function(output,vname,xmin=NULL,xmax=NULL,bw=0.2,shade=FALSE){
  
  if(shade){
    g <- ggplot(output,aes(x=output[,vname])) + 
      geom_histogram(binwidth=bw, color="black",na.rm=T)
    ncuts <- nrow(ggplot_build(g)$data[[1]])
    output$fillz <- cut(output$R_Squared,ncuts)
    col = NA
  }else{
    fillz <- "grey35"
    col <- "black"
  }
  
  lims <- c(ifelse(is.null(xmin),min(output[,vname]),xmin),
            ifelse(is.null(xmax),max(output[,vname]),xmax))
  
  ggplot(output,aes(x=output[,vname],fill=fillz)) + 
    geom_histogram(binwidth=bw,color=col,show.legend = F,na.rm=T) +
    scale_fill_grey() +
    theme_bw() +
    xlim(lims) +
    xlab("Coefficient Value") +
    ylab("Frequency") + 
    annotate("table",-Inf,Inf,
             label=list(stats(output[,vname])),
             hjust=0,vjust=1)
}



noise_plot <- function(output,vname,adj=0.3){
  rsqs <- seq(min(output$R_Squared) + 0.0001,1,length=50)
  
  ppos <- unlist(lapply(rsqs,function(x)nrow(
    output[which(output$R_Squared < x & output[,vname] > 0),])/nrow(
      output[which(output$R_Squared < x),])))
  
  sgn <- sign(median(output[,vname]))
  if(sgn < 0){
    ppos <- 1-ppos
  }
  infoloss <- unlist(lapply(rsqs,function(x)nrow(output[which(output$R_Squared < x),])))
  infoloss <- 1- infoloss/nrow(output)
  
  tmp <- data.frame(rsqs,ppos,infoloss)
  tmp %>% ggplot(aes(x=rsqs,y=ppos)) +
    theme_bw() +
    geom_hline(aes(yintercept = min(ppos),color="Maximal Uncertainty"),linetype="dashed") +
    scale_y_continuous(name="Proportion Same Sign as Naive",limits = c(adj,1),
                       sec.axis=sec_axis(~ (.-adj)/(1-adj),
                                         name="Proportion Draws Rejected")) +
    xlab("Maximum R Squared") +
    geom_point(aes(color="Proportion Same Sign")) +
    geom_point(aes(y=infoloss/(1/(1-adj))+adj,color="Draws Rejected"),shape=1) +
    scale_color_manual(name="",values=c("black","black","black"),
                       labels=c("Maximal Uncertainty","Proportion Same Sign","Draws Rejected"),
                       guide="legend") +
    guides(colour = guide_legend(override.aes = list(linetype=c(2,0,0),
                                                     shape=c(NA,19,1)))) +
    theme(legend.justification = c(0,0),legend.position = c(0,0),
          legend.background = element_blank(), legend.title = element_blank(),
          legend.box.background = element_rect(colour = "black"))
}


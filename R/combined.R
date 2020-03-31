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
  
  Is <- c(sqrt(1/Is^2),1)
  diag(Is) %*% out %*% t(diag(Is))
}

augmentcpp <- function(covm,buff=sqrt(.Machine$double.eps)){
  
    Is <- sqrt(1/diag(covm))
    corm <- diag(Is) %*% covm %*% t(diag(Is))
    
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

    Is <- c(sqrt(1/Is^2),1)
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

DOPE <- function(mod,nsims=10000,language="cpp",n.cores=1,buff=sqrt(.Machine$double.eps)){
  
  output <- list()
  mm <- model.matrix(mod)
  mmn <- colnames(mm)
  g <- grep("(Intercept)",mmn,fixed=T)
  
  mod_mat <- as.matrix(data.frame(y=model.frame(mod)[,1],if(length(g)==0){mm}else{mm[,-g]}))
  names <- c(if(length(g)==0){mmn}else{mmn[-g]},"ControlFunction","R_Squared")
  vcvm <- cov(mod_mat)
  
  require(foreach)
  require(doSNOW)
  require(doParallel)
  pb <- txtProgressBar(max = nsims, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cl <- parallel::makeCluster(n.cores)
  parallel::clusterEvalQ(cl,library(DOPE))
  parallel::clusterExport(cl,"vcvm",envir=environment())
  registerDoSNOW(cl)
  
  if(language=="cpp"){
    out <- foreach(i=1:nsims,
                   .combine=rbind,
                   .options.snow = opts) %dopar% try(simfuncpp(vcvm,buff))
  }
  if(language=="R"){
    out <- foreach(i=1:nsims,
                   .combine=rbind,
                   .options.snow = opts) %dopar% try(simfun(vcvm,buff))
  }
  suppressWarnings(try(stopCluster(cl)))
  
  colnames(out) <- names
  
  if(length(g)!=0){
    tmp <- apply(mod_mat,2,mean)
    int <- sapply(1:nrow(out),function(x)
      tmp[1] - tmp[2:ncol(mod_mat)] %*% as.matrix(out[x,1:(ncol(out)-2)])) 
    out <- data.frame(Intercept = int,out)
  }else{
    out <- data.frame(out)
  }
  
  old <- c(coef(mod)[1],coef(mod)[-1],NA,summary(mod)$r.squared)
  
  out[(nrow(out)+1),] <- old
  rownames(out)[nrow(out)] <- "Naive"
  
  out
}

simfuncpp <- function(vcvm,buff=sqrt(.Machine$double.eps)){
                psd <- FALSE
                attempts <- 0
                while(psd == FALSE){
                  attempts <- attempts + 1
                  if(attempts > 10){
                    print("Ill Conditioned System")
                    break
                  }
                  aug <- augmentcpp(vcvm,buff)
                  psd <- !any(eigen(aug)$values < 1e-8)
                }
                
                zz <- aug[-1,-1]
                zy <- as.matrix(aug[1,2:ncol(aug)])
                betas <- solve(zz) %*% zy
                rsq <- (t(zy) %*% solve(zz) %*% zy)/aug[1,1]
                output <- c(t(betas),rsq)
                output
}

simfun <- function(vcvm,buff=sqrt(.Machine$double.eps)){
                psd <- FALSE
                attempts <- 0
                while(psd == FALSE){
                  attempts <- attempts + 1
                  if(attempts > 10){
                    print("Ill Conditioned System")
                    break
                  }
                  aug <- augment(vcvm,buff)
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


plot_DOPE <- function(output,vname,xmin=NULL,xmax=NULL,bw=NULL,shade=FALSE,include_naive = TRUE){
  
  cond <- which(is.na(output$ControlFunction))
  if(length(cond)==0){
    include_naive <- FALSE
  }
  if(include_naive){
    old <- output[cond,]
    output <- output[-cond,]  
  }
  
  lims <- c(ifelse(is.null(xmin),quantile(output[,vname],probs=0.02),xmin),
            ifelse(is.null(xmax),quantile(output[,vname],probs=0.98),xmax))
  
  bw <- ifelse(is.null(bw),abs(range(lims))/50,bw)
    
  tmp <- output[which(output[,vname] >= lims[1] & output[,vname] <= lims[2]),]
  
  if(shade){
    g <- ggplot(tmp,aes(x=tmp[,vname])) + 
      geom_histogram(binwidth=bw, color="black",na.rm=T)
    ncuts <- nrow(ggplot_build(g)$data[[1]])
    tmp$fillz <- cut(tmp$R_Squared,ncuts)
    col = NA
  }else{
    fillz <- "grey35"
    col <- "black"
  }
  
  ggplot(tmp,aes(x=tmp[,vname],fill=fillz)) + 
    geom_histogram(binwidth=bw,color=col,show.legend = F,na.rm=T) +
    scale_fill_grey() +
    theme_bw() +
    xlim(lims) +
    xlab("Coefficient Value") +
    ylab("Frequency") + 
    annotate("table",-Inf,Inf,
             label=list(stats(output[,vname])),
             hjust=0,vjust=1) -> p1
  if(include_naive){
    p1 + geom_vline(xintercept = old[,vname],color="red",size=1.25)
  }else{
    p1
  }
}

sensitivity_plot <- function(output,vname,adj=NULL){
  
  cond <- which(is.na(output$ControlFunction))
  if(length(cond)!=0){
    old <- output[cond,]
    output <- output[-cond,]  
  }else{
    old <- output[0,]
  }
  
  rsqs <- seq(min(output$R_Squared) + 0.0001,1,length=50)
  
  ppos <- unlist(lapply(rsqs,function(x)nrow(
    output[which(output$R_Squared < x & output[,vname] > 0),])/nrow(
      output[which(output$R_Squared < x),])))
  
  ppos2 <- unlist(lapply(rsqs,function(x)nrow(
    output[which(output$R_Squared >= x & output[,vname] > 0),])/nrow(
      output[which(output$R_Squared >= x),])))
  
  sgn <- ifelse(nrow(old)==0,sign(median(output[,vname])),sign(old[,vname]))
  if(sgn < 0){
    ppos <- 1-ppos
    ppos2 <- 1-ppos2
  }
  
  adj <- ifelse(is.null(adj),min(ppos2,na.rm=T)-0.75*(1-min(ppos2,na.rm=T)),adj)
  
  infoloss <- unlist(lapply(rsqs,function(x)nrow(output[which(output$R_Squared < x),])))
  infoloss <- 1- infoloss/nrow(output)
  ila1 <- infoloss*(1-adj)
  ila1 <- ila1 / max(ila1,na.rm=T) * (min(ppos2,na.rm=T) - adj) + adj
  
  infoloss2 <- unlist(lapply(rsqs,function(x)nrow(output[which(output$R_Squared >= x),])))
  infoloss2 <- 1- infoloss2/nrow(output)
  ila2 <- infoloss2*(1-adj)
  ila2 <- ila2 / max(ila2,na.rm=T) * (min(ppos2,na.rm=T) - adj) + adj
  
  tmp <- data.frame(rsqs,ppos,ila1,ila2)
  
  tmp %>% ggplot(aes(x=rsqs,y=ppos)) +
    annotate(geom = "rect", xmin = min(rsqs,na.rm=T), xmax = max(rsqs,na.rm=T), ymin = adj, ymax = min(ppos2,na.rm=T),
             fill = "grey", colour = "white", alpha = 0.5) +
    theme_bw() +
    geom_hline(aes(yintercept = min(ppos,na.rm=T),color="Ignorance"),linetype="dashed") +
    geom_hline(aes(yintercept = min(ppos2,na.rm=T),color="Pessimism"),linetype="dotdash") +
    scale_y_continuous(name="Proportion Same Sign as Naive",limits = c(adj,max(ppos,na.rm=T)),
                       sec.axis=sec_axis(~ (.-adj)/(min(ppos2,na.rm=T)-adj),
                                         name="Draws Rejected",
                                         breaks=c(0,0.25,0.5,0.75,1))) +
    theme(axis.title.y.right = element_text(hjust = 1)) +
    xlab("Maximum R Squared") +
    geom_point(aes(color="Proportion Same Sign 1")) +
    geom_point(aes(y=ppos2,color="Proportion Same Sign 2"),shape=10,na.rm=T) +
    geom_point(aes(y=ila1,color="Draws Rejected (Ig)"),shape=1) +
    geom_point(aes(y=ila2,color="Draws Rejected (Ps)"),shape=18) +
    scale_color_manual(name="",values=c("black","black","black","black","black","black"),
                       labels=c(paste0("Ignorance: ",round(min(ppos,na.rm=T),3)),
                                paste0("Pessimism: ",round(min(ppos2,na.rm=T),3)),
                                "Upper Thresholding",
                                "Lower thresholding",
                                "Draws Rejected (UT)",
                                "Draws Rejected (LT)"),
                       guide="legend") +
    guides(colour = guide_legend(override.aes = list(linetype=c(2,5,0,0,0,0),
                                                     shape=c(NA,NA,19,10,1,18)))) +
    theme(legend.position = "bottom")
}

DOPE_irls <- function(X, y, family=binomial(link="logit"), maxit=25, tol=1e-08){
  X <- as.matrix(data.frame(Intercept = 1,X))
  b = rep(0,ncol(X))
  for(j in 1:maxit){
    eta    = X %*% b
    g      = family$linkinv(eta)
    gprime = family$mu.eta(eta)
    z      = eta + (y - g) / gprime
    W      = as.vector(gprime^2 / family$variance(g))
    bold   = b
    b      = solve(crossprod(X,W*X), crossprod(X,W*z), tol=2*.Machine$double.eps)
    if(sqrt(crossprod(b-bold)) < tol) break
  }
  tz <- solve(chol(solve(diag(W)))) %*% z
  tx <- solve(chol(solve(diag(W)))) %*% X
  mod <- lm(tz~.-1,data=as.data.frame(tx))
  
  list(model = mod, iterations = j, working_variable = z, weights = W)
}
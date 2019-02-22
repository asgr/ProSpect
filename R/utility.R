interp_param=function(x, params, log=FALSE, method='linear'){
  flag=rep(2,length(x))
  flag[x<min(params)]=1
  flag[x>max(params)]=3
  x[x<min(params)]=min(params)
  x[x>max(params)]=max(params)
  if(log){
    temp=interp1(log10(params), 1:length(params), xi=log10(x), method=method)
  }else{
    temp=interp1(params, 1:length(params), xi=x, method=method)
  }
  flag[temp%%1==0 & flag==2]=0
  ID_lo=floor(temp)
  ID_hi=ceiling(temp)
  ID_mode=ID_lo
  ID_mode[temp%%1>0.5]=ID_hi[temp%%1>0.5]
  return(data.frame(x=x, param_lo=params[floor(temp)], param_hi=params[ceiling(temp)], ID_lo=ID_lo, ID_hi=ID_hi, ID_mode=ID_mode, weight_lo=1-temp%%1, weight_hi=temp%%1, flag=flag))
}

#This is a direct copy of the interval function from LaplacesDemon. Since I only use this one function I didn't want to add the whole LD dependency

.interval=function (x, a = -Inf, b = Inf, reflect = TRUE) 
{
    if (missing(x)) 
        stop("The x argument is required.")
    if (a > b) 
        stop("a > b.")
    if (reflect & is.finite(a) & is.finite(b) & any(!is.finite(x))) {
        if (is.array(x)) {
            d <- dim(x)
            x <- as.vector(x)
        }
        x.inf.pos <- !is.finite(x)
        x[x.inf.pos] <- .interval(x[x.inf.pos], a, b, reflect = FALSE)
        if (is.array(x)) 
            x <- array(x, dim = d)
    }
    if (is.vector(x) & {
        length(x) == 1
    }) {
        if (reflect == FALSE) 
            x <- max(a, min(b, x))
        else if (x < a | x > b) {
            out <- TRUE
            while (out) {
                if (x < a) 
                  x <- a + a - x
                if (x > b) 
                  x <- b + b - x
                if (x >= a & x <= b) 
                  out <- FALSE
            }
        }
    }
    else if (is.vector(x) & {
        length(x) > 1
    }) {
        if (reflect == FALSE) {
            x.num <- which(x < a)
            x[x.num] <- a
            x.num <- which(x > b)
            x[x.num] <- b
        }
        else if (any(x < a) | any(x > b)) {
            out <- TRUE
            while (out) {
                x.num <- which(x < a)
                x[x.num] <- a + a - x[x.num]
                x.num <- which(x > b)
                x[x.num] <- b + b - x[x.num]
                if (all(x >= a) & all(x <= b)) 
                  out <- FALSE
            }
        }
    }
    else if (is.array(x)) {
        d <- dim(x)
        x <- as.vector(x)
        if (reflect == FALSE) {
            x.num <- which(x < a)
            x[x.num] <- a
            x.num <- which(x > b)
            x[x.num] <- b
        }
        else if (any(x < a) | any(x > b)) {
            out <- TRUE
            while (out) {
                x.num <- which(x < a)
                x[x.num] <- a + a - x[x.num]
                x.num <- which(x > b)
                x[x.num] <- b + b - x[x.num]
                if (all(x >= a) & all(x <= b)) 
                  out <- FALSE
            }
        }
        x <- array(x, dim = d)
    }
    return(x)
}

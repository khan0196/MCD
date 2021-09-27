
######################################
# MCD: Mutiple Change Points Detection
# W: window size
# M: bootstrap sample size
# lag.max:  max lag
# 2*r: number of blocks
# s: start
# e: end
# output: locations of change points
######################################
MCD = function(Y,W,M,lag.max,r,tau,s,e){
  
  if(!is.matrix(Y)) print('Data is not a matrix')
  
  Y = Y[s:e,] # truncate data set
  T = nrow(Y) # length
  d = ncol(Y) # dimension
  
  cp = vector(mode = 'list',length = lag.max+1)
  name = vector(length = lag.max+1 )
  
  if(e-s+1<2*W) print('Invalid Window Size')
  
  for(l in 0:lag.max) {
    
    H = array(0,c(d,d,T-l)) # outer product. Here domain of H is [1, T-l]  while in paper it is [1+l, T]
    
    for (i in (1+l):T) {
      H[,,i-l] = Y[i-l,]%*%t(Y[i,])
    }
    
    tau.mat = matrix(tau,nrow = d,ncol = d,byrow = T)
    
    psi.H =  array(NA,c(d,d,T-l)) # truncated outer product by tau , initialized as an array with same dimensions as H
    
    for (i in 1:(T-l)) {
      for (k in 1:d) {
        x = H[,k,i]
        psi.H[,k,i] = sign(x)*ifelse(abs(x)>tau.mat[,k],tau.mat[,k],abs(x))
      }
    }
    
    
    # calculate values of statistic w/ max norm through out points in [W,T-W], corresponding to  [s+W-1, e-W] when e = 1 & s = T
    
    S.max = vector(length=length((W):(T-W))) # Here S.max[1] = S.max(W) 
    
    for (x in 1:length(S.max)){
  
      S = MOSUM( (x+W-1), psi.H, W,l) 
      S.max[x] = max(abs(S))
      
    }
    

    lb = 1 # left boundary
    rb = T # right boundary
    CP = NULL # change point
    FLAG = 0 # denote whether to stop segamenting the current interval
    # FLAG = 1 means either no change points found within the interval 
    #          or new interval < 4W
    
    while (sum(FLAG)<length(FLAG)) {
      # loop over each level of segamentation
      # ex. ideally, in ith loop, we look at 2^(i-1) intervals
      # FLAG is a vector with same length as # of intervals in the current loop 
      # sum(FLAG) = length(FLAG) means no more segamentation needed on any intervals and break the loop
      
      FLAG = rep(0, times = length(lb))
      new.lb = NULL
      new.rb = NULL
      lb.ind = NULL # indices of left boundaries to be elimated before next seg
      rb.ind = NULL 
      

      for (i in 1:length(lb)){
        
        c = seq(from = lb[i]+W-1, to = rb[i]-W) 
        T.l = S.max[c-W+1] # Here T.l[1] = T_l(lb[i]) in range [s,e]
        t.hat = which.max(T.l)+lb[i]+W-1 # index of t.hat in the orginal range [s,e]

        
        T.star.that = BGMB(t.hat,d,W,l,psi.H,M,r) # doing Block-wise Gaussian Multiplier Bootstrap
        a = max(0,max(T.star.that))
        
        if(max(T.l)<a){
          FLAG[i] = 1
        }
        else if( (t.hat-1-lb[i])<2*W & (rb[i]-t.hat)<2*W ) { # if both left and right intervals are too short
          # then stop segamenting the two
          FLAG[i] = 1
          CP = append(CP,t.hat)
          lb.ind = c(lb.ind,i)
          rb.ind = c(rb.ind,i)
        }
        else if( (t.hat-1-lb[i])<2*W & (rb[i]-t.hat)>2*W ){ # if left interval is short then keep right one
          # and eliminate old lb of left interval
          # note that rb of left interval is always a CP (same for other round)
          CP = append(CP,t.hat)
          new.lb = append(new.lb,t.hat)
          lb.ind = c(lb.ind,i)
        }
        else if( (t.hat-1-lb[i])>2*W & (rb[i]-t.hat)<2*W ){ # if right interval is short then keep left one
          # and eliminate old rb of right interval
          CP = append(CP,t.hat)
          new.rb = append(new.rb,t.hat-1)
          rb.ind = c(rb.ind,i)
        }
        else{  #if both are fine, then save both
          CP = append(CP,t.hat)
          new.lb = append(new.lb,t.hat)
          new.rb = append(new.rb,t.hat-1)
        }
      }
      # removing redundant boundaries 
      if(!is.null(lb.ind)) lb = lb[-lb.ind] 
      if(!is.null(rb.ind)) rb = rb[-rb.ind]
      # updating new boundaries for next loop
      lb = sort(union(lb,new.lb))
      rb = sort(union(rb,new.rb))
      
    }
    
    name[l+1] =paste('lag =',l)
    if(!is.null(CP)){cp[[l+1]] = CP+s-1} #returning to original range of the series
  }
  names(cp) = name
  return(cp)
}

######################################
# sigma.hat: autocovariance estimate 
# psi.H: truncated outer product
# l: lag
# output: autocovariance estimate
######################################
sigma.hat = function(psi.H,start,end,l){
  
  gamma = 0
  
  for (t in start:(end-l)) {
    
    gamma = gamma + (end - start - l + 1)^(-1)*psi.H[,,t] #H[,,start] = H_(start+l)
    
  } 
  
  return(gamma)
}
######################################
# MOSUM
# c: check point in range[1,T]
# W: window size
# l :lag
# output: value of S
######################################
MOSUM = function(c,psi.H,W,l){
  
  B.s = c-W+1
  B.e = c
  
  A.s = c+1
  A.e = c+W
  stat = ((W-l)/2)^(0.5)*(sigma.hat(psi.H,B.s,B.e,l)  - sigma.hat(psi.H,A.s,A.e,l))
  
  return(stat)
  
}
################################################
# BGMB: Block-wise Gaussian Multiplier Bootstrap
# v: check point
# d: dimension
# W: window size
# l: lag
# H: outer product
# M: bootstrap sample size
# r: half # of blocks
# output: T.star at chk point v
################################################
BGMB = function(v,d,W,l,psi.H,M,r){
  
  size=floor((W-l)/(2*r))
  if(size==0) print('Invalid Block Size Found')
  
  O = array(0,c(d,d,2*r))
  E = array(0,c(d,d,2*r))
  
  for(q in 1:r){ # within each block,there are # = size check points
    
    index.left.O = (v-W+2*(q-1)*size):(v-W+(2*q-1)*size-1)
    index.left.E = (v-W+(2*q-1)*size):(v-W+2*q*size-1)
    
    u = q + r
    
    index.right.O = (v-W+l+1+2*(u-1)*size):(v-W+l+1+(2*u-1)*size-1)            
    index.right.E = (v-W+l+1+(2*u-1)*size):(v-W+l+1+2*u*size-1)                             
    
    O[,,q] =  apply(psi.H[,,index.left.O],1:2,sum)
    O[,,u] =  apply(psi.H[,,index.right.O],1:2,sum)
    
    E[,,q] =  apply(psi.H[,,index.left.E],1:2,sum)
    E[,,u] =  apply(psi.H[,,index.right.E],1:2,sum)
    
  }
  
  D = (O - E)
  
  
  T.star = vector(length = M)
  
  for (n in 1:M) {
    
    set.seed(124+n)
    e = rnorm(2*r)
    S.star = matrix(0,nrow = d,ncol = d)
    
    for (i in 1:(2*r)) {
      S.star = S.star + e[i]*D[,,i]
    }
    
    T.star[n] = max(abs((2*(W-l))^(-0.5)*S.star))
  }
  
  return(T.star)
  
}



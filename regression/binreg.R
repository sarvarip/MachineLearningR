sigm = function(x){
  return (1/(1+exp(-x)))
}

evaluate = function(y,mu){
  norows = length(y)
  mse = 1/norows * (sum((y-mu)^2))
  return (mse)
}

updateweight = function(X,n,y,lam,weight){
 #X should have the column of ones prepended
 nofeat = ncol(X)
 theta = X %*% weight
 mu = n * sigm(theta)
 sigmoid = sigm(theta)
 diagS = as.vector(n * sigmoid * (1 - sigmoid))
 S = diag(diagS)
 hessian = - t(X) %*% S %*% X - 2 * lam * diag(nofeat)
 gradient = t(X) %*% (y - mu) - 2 * lam * weight
 update = solve(hessian) %*% gradient
 #update = normalize(update)
 return (weight - update)
}

binreg = function(X,n,y,lam,iter){
  i = 0
  norows = length(y)
  std = apply(X, 2, sd)
  avg = apply(X, 2, mean)
  Xscaled = scale(X)
  Xscaled = cbind(rep(1, norows), Xscaled)
  nofeat = ncol(Xscaled)
  weight = rnorm(nofeat)
  ypred = binpredict(Xscaled,n,weight)
  mse.baseline = evaluate(y,ypred)
  mse.baseline = mse.baseline / 3
  #We expect to mse to be reduced at least 
  #by 67% after correct Newton's iteration
  #print(mse.baseline)
  while(i<iter){
    #print(weight)
    weight = updateweight(Xscaled,n,y,lam,weight)
    ypred = binpredict(Xscaled,n,weight)
    mse = evaluate(y,ypred)
    #print(mse)
    i = i + 1
  }
  #transform weights back as if was unscaled
  weight[2:nofeat] = weight[2:nofeat] / std
  weightbias = weight[1] - sum(avg*weight[2:nofeat]) 
  weight[1] = weightbias
  if (mse>mse.baseline){
    #In case Newton's does not converge
    #just run again
    print("Did not converge, retrying...")
    weight = binreg(X,n,y,lam,iter)
  }
  return (weight)
}

binpredict = function(X,n,weight){
  norows = length(n)
  nofeat = length(weight)
  if (ncol(X) != nofeat){
    X = cbind(rep(1, norows), X)
  }
  mu = n * sigm(X %*% weight)
  return (mu)
}

cv.binpredict = function(X,n,y){
  #5-fold cross-validation
  norows = length(n)
  iter = 10*norows
  parts <- sample(rep(1:5, length.out=norows))
  lam.totry = exp(c(-20,seq(-10,1,1)))
  #exp(-20) essentially represents 
  #zero regularization but helps with
  #numerical issues when inverting hessian
  mse.lambda = rep(0,length(lam.totry))
  for (lambda in seq(lam.totry)){
    lam = lam.totry[lambda]
    for (j in 1:5){
      trainind = which(parts != j)
      testind = which(parts != j)
      train = X[trainind,]
      test = X[testind,]
      n.train = n[trainind]
      y.train = y[trainind]
      n.test = n[testind]
      y.test = y[testind]
      mse.fold = rep(0,5)
      
      weight = binreg(train,n.train,y.train,lam,iter)
      ypred = binpredict(test,n.test,weight)
      mse = evaluate(y.test,ypred)
      mse.fold[j] = mse
    }
    mse.lambda[lambda] = mean(mse.fold)
  }
  best.lambda.idx = which.min(mse.lambda)
  best.lambda = lam.totry[best.lambda.idx]
  return (best.lambda)
}

data.train = matrix(c(0.3, 0.8, 0.2, 0.1, 0.1, 0.2, 0.4, 0.9, 0.1, 0.8, 0.9, 0.1, 0.5, 0.1, 0.7), 5, 3)
y.train = c(0, 4, 10, 5, 1)
n.train = c(10, 14, 12, 6, 8)
iter = 10
lam = 0.001

#weight = binreg(data.train,n.train,y.train,lam,iter)
#ypred = binpredict(data.train,n.train,weight)
#mse = evaluate(y.train,ypred)

best.lambda = cv.binpredict(data.train, n.train, y.train)
iter = nrow(data.train)
weight = binreg(data.train,n.train,y.train,best.lambda,iter)
ypred = binpredict(data.train,n.train,weight)
mse = evaluate(y.train,ypred)
print(mse)

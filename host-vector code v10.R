######################Meta Data######################
#v8 updates:
#Both host and Vector have SI infection model
#Infection starts in a single species, simulations are run with an infection start in each host species in a network
#Infection simulations record the prevalence of infection at each time step, up to 400 steps.

#v10 updates:
#checkLambda variable in functions simulateData, simHostTrait, and simVectorTrait
#checkLambda resimulates data if the difference between the estimated lambda of the simulated data and the imput lambda parameter is greater than checkLambda
#checkLambda is default Null, should be a number between 0 and 1
#no contraints limiting the value of checkLambda yet
######################Set Up#########################

library(geiger)
library(PHYLOGR)
library(igraph)
library(bipartite)
library(phytools)
library(ggplot2)
library(lhs)
library(MuMIn)


###################PrefMatrix Funtions################################
createTree = function(birthrate, deathrate, numSpecies) {
	tree = sim.bdtree(b = birthrate, d = deathrate, stop = "taxa", n = numSpecies, t = 100, extinct = FALSE) #create tree
	tree = drop.extinct(tree) 		#prune extinct lineages
	for (y in 1: length(tree$edge.length)) {							#remove polytomies
		if (tree$edge.length[y] == 0.00000000) tree$edgelength[y] = 0.01
	}
	tree
}

createVCV = function(var1, var2, CoV) {
#Create a variance-covariance matrix for two traits
	matrix(c(var1, CoV, CoV, var2), nrow = 2, ncol = 2)
}

simCoVaryTraits = function(treeNum, h1Var, h2Var, hCoV, lambda.host, allTrees) {
	VCV = createVCV(h1Var, h2Var, hCoV)
	treeAdj = rescale(allTrees[[treeNum]], "lambda", lambda.host) 
	treeData = sim.char(treeAdj, VCV, nsim = 1, model = "BM", root = 0)
	treeData = treeData[,,1]
	rownames(treeData) = gsub("s", "h", rownames(treeData))
	treeData
}

simHostTrait = function(treeNum, hVar, lambda.host, allTrees, checkLambda = NULL) {  ##v10: added checkLambda
  tree = allTrees[[treeNum]]
  treeAdj = rescale(tree, "lambda", lambda.host)
  treeData = sim.char(treeAdj, hVar, nsim = 1, model = "BM", root = 0)
  treeData = treeData[,,1]
  if (! is.null(checkLambda)) {  ##V10: following lines check estimate lambda from the generated data, and regenerate data while the estimated lambda is further from the input lambda than checkLambda
    estLambda = phylosig(tree, treeData, method = "lambda")$lambda
    ##print(c(lambda.host, estLambda)) ##debug for checkLambda
    while (abs(estLambda- lambda.host) > checkLambda) {
      ##print(1) ##debug for checkLambda
      treeData = sim.char(treeAdj, hVar, nsim = 1, model = "BM", root = 0)
      treeData = treeData[,,1]
      estLambda = phylosig(tree, treeData, method = "lambda")$lambda
      ##print(c(lambda.host, estLambda)) ##debug for checkLambda
    }
  }
  names(treeData) = gsub("s", "h", names(treeData))
  treeData
}


simVectorTrait = function(treeNum, v1Var, lambda.vector, allTrees, checkLambda = NULL) { ##v10:see comments in simHostTrait
	tree = allTrees[[treeNum]]
  treeAdj = rescale(tree, "lambda", lambda.vector)
	treeData = sim.char(treeAdj, v1Var, nsim = 1, model = "BM", root = 0)
	treeData = treeData[,,1]
	if (! is.null(checkLambda)) {  ##V10: following lines check estimate lambda from the generated data, and regenerate data while the estimated lambda is further from the input lambda than checkLambda
	  estLambda = phylosig(tree, treeData, method = "lambda")$lambda
	  ##print(c(lambda.vector, estLambda)) ##debug for checkLambda
	  while (abs(estLambda- lambda.vector) > checkLambda) {
	    ##print(2)  ##debug for checkLambda
	    treeData = sim.char(treeAdj, v1Var, nsim = 1, model = "BM", root = 0)
	    treeData = treeData[,,1]
	    estLambda = phylosig(tree, treeData, method = "lambda")$lambda
	    ##print(c(lambda.vector, estLambda)) debug for checkLambda
	  }
	}
	names(treeData) = gsub("s", "v", names(treeData))
	treeData
}

createPrefMatrix = function(hTrait, vTrait, g) {
#input: vTrait is a vector of simulated vector trait, and g is the vector specialization index
	hTraitScaled = ((hTrait-min(hTrait))/max(hTrait-min(hTrait)))  #scales all values in hTrait vector from 0 to 1
	vTraitScaled = ((vTrait-min(vTrait))/max(vTrait-min(vTrait)))  #scales all values in vTrait vector from 0 to 1
	diffMatrix = matrix(nrow = length(hTraitScaled), ncol = length(vTraitScaled), dimnames = list(names(hTraitScaled), names(vTraitScaled)))
	for (i in 1:length(hTraitScaled)) {
		for (j in 1:length(vTraitScaled)) {
			diffMatrix[i, j] = abs(hTraitScaled[i] - vTraitScaled[j]) #fills each cell of the matrix with the difference between scaled host value and the scaled vector value
		}
	}
	prefMatrix = apply(diffMatrix, MARGIN = c(1,2), FUN = function(X) (1-X)^exp(g))
	prefMatrix
}

getModularity = function(prefMatrix) {
	modularity = computeModules(prefMatrix)
	modularity@likelihood
}


extinction = function(h1data, h2data, prop.extinction) { #input:preferred trait data, extinction trait data, and proportion of species to go extinct; output: new preferred trait data
	numExtinct = floor(length(h2data) * prop.extinction)
	hExtant = names(h2data)[order(h2data)][1:(length(h2data) - numExtinct)]
	h1data[hExtant]
}

createHypercube = function(numSim = 100, h = c(20), b.h = c(0.9), d.h = c(0.5), v = c(20), b.v = c(.9), d.v = c(.5), lambda.vector = c(1), lambda.host = c(0,1), lambda.extinction = c(0,1), h1Var = c(1), h2Var = c(1), hCoV = c(0), v1Var = c(1), g = c(0, 2), prop.extinct = c(0.5)) { #lambda.vector = lambda for vector preference, lambda.host = lambda for host preferred trait
	
	#Parameter values are desired range for hypercube. To exclude a value from the hypercube, enter a single fixed value for that paramter
	
	P = list(h = h, b.h = b.h, d.h = d.h, v = v, b.v = b.v, d.v = d.v, h1Var = h1Var, h2Var = h2Var, hCoV = hCoV, v1Var = v1Var, lambda.vector = lambda.vector, lambda.host = lambda.host, lambda.extinction = lambda.extinction, g = g, prop.extinct = prop.extinct)
	
	whichP = rep(0, times = length(P))
	for (i in 1:length(P)) {
		if (length(P[[i]]) == 2) {whichP[i] = 1
		} else P[[i]] = rep(P[[i]], numSim)
	}
	if (numSim>200)
		lhs = randomLHS(numSim, sum(whichP))
	else 
		lhs = optimumLHS(numSim, sum(whichP))

	colnames(lhs) = names(P)[whichP]
	for(i in 1:ncol(lhs)) {
		P[which(whichP==1)][[i]] = lhs[,i] * (P[which(whichP == 1)][[i]][2] - P[which(whichP == 1)][[i]][1]) + P[which(whichP ==1)][[i]][1]
	}
	P$h = sapply(P$h, ceiling)
	P$v = sapply(P$v, ceiling)
	P = as.data.frame(P)
	 P[order(P$d.h),]
}

simulateTrees = function(ntree, P, hypercube = FALSE) {
	if (hypercube) {
		treesMatrix = matrix(nrow = ntree * nrow(P), ncol = 8) #8 columns represent host tree, vector tree, and the 7 parameters used to simulate
		colnames(treesMatrix) = c("hTree", "vTree", "h", "b.h", "d.h", "v", "b.v", "d.v")
		trees = data.frame(treesMatrix)
		treesP = matrix(nrow = 0, ncol = 8)
		colnames(treesP) = c("hTree", "vTree", "h", "b.h", "d.h", "v", "b.v", "d.v")
		for (x in 1:ntree) treesP = rbind(treesP, P)
		trees[,which(colnames(trees) %in% colnames(treesP))] = treesP[,which(colnames(treesP) %in% colnames(trees))]
		trees = trees[order(trees$h),]
		rownames(trees) = 1:nrow(trees)
	} else {
		treesMatrix = matrix(nrow = 0, ncol = 8) #8 columns represent host tree, vector tree, and the 7 parameters used to simulate
		colnames(treesMatrix) = c("hTree", "vTree", "h", "b.h", "d.h", "v", "b.v", "d.v")
		trees = data.frame(treesMatrix)
		for(h in 1:length(P$h)) {
			for(b.h in 1:length(P$b.h)) {
			for(d.h in 1:length(P$d.h)) {
			for(v in 1:length(P$v)) {
			for(b.v in 1:length(P$b.v)) {
			for(d.v in 1:length(P$d.v)) {
				addTrees = data.frame(hTree = rep(NA, times = ntree), vTree = NA, h = P$h[h], b.h = P$b.h[b.h], d.h = P$d.h[d.h], v = P$v[v], b.v = P$b.v[b.v], d.v = P$d.v[d.v])
				trees = rbind(trees, addTrees)
	}}}}}}}
	trees$hTree = vector("list", length = length(trees$hTree))
	trees$vTree	= vector("list", length = length(trees$vTree))

	trees$hTree = mapply(createTree, birthrate = trees$b.h, deathrate = trees$d.h, numSpecies = trees$h, SIMPLIFY = FALSE) 
	trees$vTree = mapply(createTree, birthrate = trees$b.v, deathrate = trees$d.v, numSpecies = trees$v, SIMPLIFY = FALSE)
	return(trees)
}

simulateData = function(nsim, P, trees, hypercube = FALSE, checkLambda = NULL) {
	if (hypercube) {
		oneSim = data.frame(treeNum = 1:nrow(trees), h1data = NA, h2data = NA, vdata = NA, prefMatrixAll = NA, prefMatrixExt = NA, h1Var = NA, h2Var = NA, hCoV = NA, v1Var = NA, lambda.vector = NA, lambda.host = NA, lambda.extinction = NA, g = NA, prop.extinct = NA)
		for (i in 1:nrow(oneSim)) {
			oneSim[i, -(1:6)] = P[ceiling(i/(nrow(trees)/nrow(P))), -(1:6)]
		}
		allData = data.frame(treeNum = NA, h1data = NA, h2data = NA, vdata = NA, prefMatrixAll = NA, prefMatrixExt = NA, h1Var = NA, h2Var = NA, hCoV = NA, v1Var = NA, lambda.vector = NA, lambda.host = NA, lambda.extinction = NA, g = NA, prop.extinct = NA)
		for (rep in 1:nsim) {
			allData = rbind(oneSim, allData)
		}
		allData = allData[-nrow(allData),]
	} else {		

		allData = matrix(NA, nrow = 0, ncol = 15) #1 col each for tree number, pref matrix with and without extinction, each data simulation parameter, plus 2 for host traits and 1 for vector data,
		colnames(allData) = c("treeNum", "h1data", "h2data", "hCoV", "vdata", "prefMatrixAll", "prefMatrixExt", "h1Var", "h2Var", "lambda.vector", "v1Var", "lambda.host", "lambda.extinction", "g", "prop.extinction")		
		allData = data.frame(allData)


		for(h1Var in P$h1Var) {
			for (h2Var in P$h2Var) {
			for (hCoV in P$hCoV) {
			for (v1Var in P$v1Var) {
			for (lambda.vector in P$lambda.vector){
			for(lambda.host in P$lambda.host) {
			for(lambda.extinction in P$lambda.extinction){
			for (g in P$g) {
			for (prop.extinct in P$prop.extinct) {
			for (treeNum in 1:nrow(trees)) {
				addData = data.frame(treeNum = treeNum, h1data = rep(NA, times = nsim), h2data = NA, vdata = NA, prefMatrixAll = NA, prefMatrixExt = NA, h1Var = h1Var, h2Var = h2Var, hCoV = hCoV, v1Var = v1Var, lambda.vector = lambda.vector, lambda.host = lambda.host, lambda.extinction = lambda.extinction, g = g, prop.extinct = prop.extinct)
				allData = rbind(allData, addData)
		}}}}}}}}}}
	}

	#hostTraits = mapply(simHostTraits, treeNum = allData$treeNum, h1Var = allData$h1Var, h2Var = allData$h2Var, hCoV = allData$hCoV, lambda.host = allData$lambda.host, MoreArgs = list(allTrees = trees$hTree), SIMPLIFY = FALSE) for covarying traits 
	#allData$h1data = mapply(function(x) x[,1], hostTraits, SIMPLIFY = FALSE)
	#allData$h2data = mapply(function(x) x[,2], hostTraits, SIMPLIFY = FALSE)
  allData$h1data = mapply(simHostTrait, treeNum = allData$treeNum, hVar = allData$h1Var, lambda.host = allData$lambda.host, MoreArgs = list(allTrees = trees$hTree, checkLambda = checkLambda), SIMPLIFY = FALSE) ##v10: added checkLambda
  allData$h2data = mapply(simHostTrait, treeNum = allData$treeNum, hVar = allData$h2Var, lambda.host = allData$lambda.extinction, MoreArgs = list(allTrees = trees$hTree, checkLambda = checkLambda), SIMPLIFY = FALSE) #v10: added checkLambda
  allData$vdata = mapply(simVectorTrait, treeNum = allData$treeNum, v1Var = allData$v1Var, lambda.vector = allData$lambda.vector, MoreArgs = list(allTrees = trees$vTree, checkLambda = checkLambda), SIMPLIFY = FALSE) 	 #v10: added checkLambda			
	allData$prefMatrixAll = mapply(createPrefMatrix, hTrait = allData$h1data, vTrait = allData$vdata, g = allData$g, SIMPLIFY = FALSE)
	h1Extinct = mapply(extinction, h1data = allData$h1data, h2data = allData$h2data, prop.extinction = allData$prop.extinct, SIMPLIFY = FALSE)
	allData$prefMatrixExt = mapply(createPrefMatrix, hTrait = h1Extinct, vTrait = allData$vdata, g = allData$g, SIMPLIFY = FALSE)
	return(allData)
}

calculateZ = function(prefMatrix, N = 100) {
	nulls <- nullmodel(prefMatrix, N=N, method="r2d")
	modules.nulls <- sapply(nulls, computeModules)
	like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
	z <- (getModularity(prefMatrix) - mean(like.nulls))/sd(like.nulls)
	z
}

getQscores = function(traitData, trees, N = 100) {
	allData = list(traitData = traitData, trees = trees)
	#allData$traitData$Q_All = sapply(traitData$prefMatrixAll, calculateZ, N = N)
	#allData$traitData$Q_Ext = sapply(traitData$prefMatrixExt, calculateZ, N = N)
	allData$traitData$Q_All = sapply(traitData$prefMatrixAll, getModularity)
	allData$traitData$Q_Ext = sapply(traitData$prefMatrixExt, getModularity)
	allData
}
	
	
##############################################################################

######################Infection Model Functions###############################
allInfections = function(allData, infect_0 = 0.01, vPop = 200, hPop = 100, hInfRate = 0.05, vInfRate = 0.05, steps = 300) {  #set infection simulation parameters here
  allData$traitData$prev_All = lapply(allData$traitData$prefMatrixAll, runInfection, infect_0 = infect_0, vPop = vPop, hPop = hPop, hInfRate = hInfRate, vInfRate = vInfRate, steps = steps)
  allData$traitData$prev_Ext = lapply(allData$traitData$prefMatrixExt, runInfection, infect_0 = infect_0, vPop = vPop, hPop = hPop, hInfRate = hInfRate, vInfRate = vInfRate, steps = steps)
  allData
}

runInfection = function(prefMatrix, infect_0 = 0.01, vPop = 200, hPop = 100, hInfRate = 0.05, vInfRate = 0.05, steps = 300) {
  allPrev = matrix(1, nrow = steps, ncol = nrow(prefMatrix))
  colnames(allPrev) = rownames(prefMatrix)
  allPrev = data.frame(allPrev)
  propMatrix = getPropMatrix(prefMatrix)
  for(host_0 in rownames(prefMatrix)) {
    hosts = data.frame(name = rownames(prefMatrix), S = hPop , I = 0)
    hosts$S[which(hosts$name == host_0)] = hosts$S[which(hosts$name == host_0)] - ceiling(hPop * infect_0)
    hosts$I[which(hosts$name == host_0)] = hosts$I[which(hosts$name == host_0)] + ceiling(hPop * infect_0)
    # print(hosts)
    vectors = data.frame(name = colnames(prefMatrix), S = vPop, I = 0)
    t = 0
    #allHostData = vector("list", 0) #stores the infection levels at every times step. Sucks system memory, only use if necessary
    #allVectorData = vector("list", 0) #see above
    while(sum(hosts$S) > 0 & t <= steps) {
      t = t + 1
      deltaHosts = apply(hosts, MARGIN = 1, FUN = getDeltaHosts, vectors = vectors, propMatrix = propMatrix, hInfRate = hInfRate)
      deltaHosts = data.frame(t(deltaHosts))
      vectorInfections = apply(vectors, MARGIN = 1, FUN = getVectorInf, hosts = hosts, propMatrix = propMatrix, vInfRate = vInfRate)
      hosts[, -1] = hosts[,-1] + deltaHosts
      vectors$S = vectors$S - as.numeric(vectorInfections)
      vectors$I = vectors$I + as.numeric(vectorInfections)
      prev = sum(hosts$I)/(sum(hosts$S) + sum(hosts$I))
      allPrev[t, host_0] = prev
      # plot(allPrev[,host_0])
      #allHostData[[t]] = hosts #see AllHostData above
      #allVectorData[[t]] = vectors #see AllHostData above
    }
  }
  print(allPrev[300,])
  return(allPrev)
}

getVectorInf = function(vector, hosts, propMatrix, vInfRate) {
  infections = sum(apply(hosts, MARGIN = 1, perHost, vector = vector, propMatrix = propMatrix, vInfRate = vInfRate))
  if (infections > as.numeric(vector["S"])) infections = vector["S"]
  infections
}

perHost = function(host, vector, propMatrix, vInfRate) {
  hostName = host["name"]
  vectorName = vector["name"]
  host = as.numeric(host[-1])
  vector = as.numeric(vector[-1])
  bites = runif(vector[1], 0 , 1)
  infections = sum(bites < propMatrix[hostName, vectorName] * (host[2]/sum(host)) * vInfRate)
  infections
}

getPropMatrix = function(prefMatrix) {
  apply(prefMatrix, MARGIN = 2, FUN = getPropCol)
}

getPropCol = function(column) {
  sapply(column, function(X) X/sum(column))
}

getDeltaHosts = function(host, vectors, propMatrix, hInfRate) {
  hostInfections = sum(apply(vectors, MARGIN = 1, FUN = getHostInfections, host = host, propMatrix = propMatrix, hInfRate = hInfRate))
  if (hostInfections > as.numeric(host["S"])) hostInfections = as.numeric(host["S"])
  return(c(0 - hostInfections, hostInfections))
}

getHostInfections = function(vector, host, propMatrix, hInfRate) {
  hostName = host["name"]
  vectorName = vector["name"]
  host = as.numeric(host[-1])
  vector = as.numeric(vector[-1])
  bites = runif(vector[2], 0, 1) 
  infections = sum(bites < propMatrix[hostName, vectorName] * (host[1]/sum(host)) * hInfRate)
  infections
} 
##############################################################################

#########################Analysis Functions###################################
getPrev = function(allPrev, t) {
  #input: dataframe of prevelence where each column represents a species and each row a time step, the time step of interest
  #returns the mean and sd of disease prevelance for all species of a network at a given time (t)
  prev = allPrev[t,]
  out = list(mean = mean(as.numeric(prev)), sd = sd(prev))
  out
}

getPrevAllSim = function(prevs, t) {
  #Gives a vector with average prevalence at a given time step for each network in all Data
  PrevAllSim = sapply(prevs, FUN = function(X, t) getPrev(X, t)[[1]], t = t)
  PrevAllSim
}

plotAllPrev = function(allPrev, interval = 1) {
  #scatterplot of prevelance data with time step on x and prevelance of each species on y. Interval determines distance between each plotted time step
  allPrev = ggplotDataTransform((allPrev))
  allPrev = allPrev[which(allPrev$t/interval == round(allPrev$t/interval)),]
  g = ggplot(allPrev, aes(t, prev)) + geom_point()
  g
}

ggplotDataTransform = function(raw) {
  #transforms a marix or dataframe with one variable on x and one on y into dataframe with one row for every observation and a column for each variable
  ggplotFormat = data.frame(matrix(nrow = 0, ncol = 3))
  colnames(ggplotFormat) = c("t", "host", "prev")
  for(t in 1:nrow(raw)) {
    for(host in colnames(raw)) {
      ggplotFormat = rbind(ggplotFormat, data.frame("t" = t, "host" = host, "prev" = raw[t, host]))
    }
  }
  ggplotFormat
}

time2saturation = function(host, allPrev, satLevel = 1) {
  #input: prevelance data for all hosts in a network, host of interest
  #output: time steps until disease saturation when starting in given host
  t = which(allPrev[,host] >= satLevel)[1]
  if (is.na(t)) t = 300
  t
}

getMeanSaturation = function(allPrev, satLevel = 1) {
  #gets the mean of saturation times for disease starting in each host species
  all_t = sapply(colnames(allPrev), time2saturation, allPrev = allPrev, satLevel = satLevel)
  return(mean(all_t))
}

getSDSaturation = function(allPrev, satLevel = 1) {
  all_t = sapply(colnames(allPrev), time2saturation, allPrev = allPrev, satLevel = satLevel)
  return(sd(all_t))
}

modelAveraging = function(allData, formula) {
  linear = lm(formula, data = allData$traitData, na.action = 'na.fail')
  d = dredge(linear, beta = 'sd')
  avg = model.avg(d)
  summary(avg)
}

summaryTable = function(allData) {
  traitData = cbind(allData$traitData, allData$trees)
  sumT = ddply(traitData, c("h", "v", "lambda.host", "lambda.vector"), summarise, N = length(Q_All), mean = mean(Q_All), sd = sd(Q_All), SE = sd/sqrt(N))
  sumT
}
  
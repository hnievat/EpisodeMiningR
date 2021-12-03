#WINEPI algorithm for R by hnievat
#Original algorithm by Mannila, et al. 1997
#Based on Python implementation by duytri and nkmrtty

library(gtools)
library(rapportools)

data <- read.csv(file.choose(), header=TRUE, row.names=1) #column 1=time unit, needs header row

#### User parameters ####

width <- 5 #window width in time units
step <- 1 #time units for window scanning
minFrequent <- 0.1 #minimum support for sequences
minConfidence <- 0.7 #minimum confidence for rules
episodeType <- "parallel" #"serial" or "parallel"

#### Functions ####

scanWindows_parallel <- function (windows, Ck) {
  ssCnt <- data.frame(val=double())
  for (i in 1:length(windows)) {
    if (length(Ck) == 0) {
      return("done")
    }
    if (typeof(Ck) == "character") { #for length 1 combos
      cklength=length(Ck)
      for (j in 1:cklength) {
        if (all(Ck[j] %in% windows[i][[1]]) == TRUE) { 
          if (any(row.names(ssCnt) == Ck[j]) == FALSE) { 
            newSubset <- data.frame(row.names=c(Ck[j]), val=c(1)) 
            ssCnt <- rbind(ssCnt, newSubset)
          } else {
            ssCnt[Ck[j],] <- ssCnt[Ck[j],] + 1 
          }
        }
      }
    }
    if (typeof(Ck) == "list") { #for length >1 combos
      cklength=length(Ck)
      for (j in 1:cklength) { 
        if (all(Ck[[j]] %in% windows[i][[1]]) == TRUE) { 
          listString <- paste(c(Ck[[j]]), collapse=" ") 
          if (any(row.names(ssCnt) == c(listString)) == FALSE) { 
            newSubset <- data.frame(row.names=c(listString), val=c(1)) 
            ssCnt <- rbind(ssCnt, newSubset)
          } else {
            ssCnt[listString,] <- ssCnt[listString,] + 1 
          }
        }
      }
    }
  }
  numItems <- length(windows)
  retList <- list()
  supportData <- matrix(ncol=2, nrow=0)
  for (key in 1:nrow(ssCnt)) {
    support <- ssCnt[key,]/numItems
    if (support >= minFrequent) {
      temp <- c(rownames(ssCnt)[key])
      retList <- rbind(temp, retList)
    }
    temp <- data.frame(row.names=c(rownames(ssCnt)[key]), val=c(support))
    supportData <- rbind(supportData, temp)
  }
  return <- list(retList,supportData)
  return(return)
}

isSubsetInOrderWithGap <- function (sub, lst) { 
  ln <- length(sub) 
  j <- 0 
  for (elem in 1:length(lst[[1]])) { 
    if (lst[[1]][elem] == sub[j+1]) { 
      j <- j + 1
    }
    if (isTRUE(j == ln)) { 
      return(TRUE) 
    }
  }
  return(FALSE)
}

scanWindows_serial <- function (windows, Ck) { 
  ssCnt <- data.frame(val=double()) 
  for (i in 1:length(windows)) { 
    if (length(Ck) == 0) {
      return("done")
    }
    if (typeof(Ck) == "character") { #for length 1 combos
      cklength=length(Ck)
      for (j in 1:cklength) { 
        if (isSubsetInOrderWithGap(Ck[j],windows[i]) == TRUE) { 
          if (any(row.names(ssCnt) == Ck[j]) == FALSE) { 
            newSubset <- data.frame(row.names=c(Ck[j]), val=c(1)) 
            ssCnt <- rbind(ssCnt, newSubset)
          } else {
            ssCnt[Ck[j],] <- ssCnt[Ck[j],] + 1 
          }
        }
      }
    }
    if (typeof(Ck) == "list") { #for length >1 combos
      cklength=length(Ck)
      for (j in 1:cklength) { 
        if (isSubsetInOrderWithGap(Ck[[j]],windows[i]) == TRUE) { 
          listString <- paste(c(Ck[[j]]), collapse=" ") 
          if (any(row.names(ssCnt) == c(listString)) == FALSE) { 
            newSubset <- data.frame(row.names=c(listString), val=c(1)) 
            ssCnt <- rbind(ssCnt, newSubset)
          } else {
            ssCnt[listString,] <- ssCnt[listString,] + 1 
          }
        }
      }
    }
  }
  numItems <- length(windows)
  retList <- list()
  supportData <- matrix(ncol=2, nrow=0)
  for (key in 1:nrow(ssCnt)) {
    support <- ssCnt[key,]/numItems
    if (support >= minFrequent) {
      temp <- c(rownames(ssCnt)[key])
      retList <- rbind(temp, retList)
    }
    temp <- data.frame(row.names=c(rownames(ssCnt)[key]), val=c(support))
    supportData <- rbind(supportData, temp)
  }
  return <- list(retList,supportData)
  return(return)
}

checkSubsetFrequency <- function (candidate, Lk, k) {
  if (k>1) {
    subsets <- combinations(length(candidate),k,candidate)
  } else {
    return(TRUE)
  }
  subsetsList <- as.list(as.data.frame(t(subsets)))
  subsetsPaste <- vector()
  for (j in 1:length(subsetsList)) {
    subsetsPaste <- c(subsetsPaste,paste(subsetsList[[j]], collapse=" "))
  }
  for (i in 1:length(subsetsPaste)) {
    if (any(Lk == subsetsPaste[i]) == FALSE) {
      return(FALSE)
    }
  }
  return(TRUE)
}

aprioriGen_parallel <- function (Lk, k) {
  resList <- list()
  candidatesK <- list()
  lk <- list()
  for (t in 1:length(Lk)) {
    for (item in 1:length(Lk[[t]])) {
      temp <- c(Lk[[t]][item])
      lk <- rbind(lk, temp)
    }
  }
  lksplit <- vector()
  for (row in 1:length(lk)) {
    lksplit <- c(lksplit, unlist(strsplit(lk[[row]], " ")))
  }
  lksplit <- unique(lksplit)
  lk <- sort(lksplit)
  if (length(lk) >= k) {
    candidatesK <- list(combinations(length(lk),k,lk))
  } else {
    candidatesK <- list()
  }
  if (length(candidatesK) > 1) {
    for (i in 1:length(candidatesK)) {
      if (checkSubsetFrequency(candidatesK[[1]][i,], Lk, k-1) == TRUE) {
        resList <- append(resList, candidatesK[i])
      }
    }
  }
  if (tryCatch({nrow(candidatesK[[1]])},error=function(e){print("0")}) > 1) {
    for (i in 1:nrow(candidatesK[[1]])) {
      if (checkSubsetFrequency(candidatesK[[1]][i,], Lk, k-1) == TRUE) {
        resList <- append(resList, list(candidatesK[[1]][i,]))
      }
    }
  }
  return(resList)
}

aprioriGen_serial <- function (Lk, k) {
  resList <- list()
  candidatesK <- list()
  lk <- list()
  for (t in 1:length(Lk)) {
    for (item in 1:length(Lk[[t]])) {
      temp <- c(Lk[[t]][item])
      lk <- rbind(lk, temp)
    }
  }
  lksplit <- vector()
  for (row in 1:length(lk)) {
    lksplit <- c(lksplit, unlist(strsplit(lk[[row]], " ")))
  }
  lksplit <- unique(lksplit)
  lk <- sort(lksplit)
  if (length(lk) >= k) {
    candidatesK <- list(permutations(length(lk),k,lk))
  } else {
    candidatesK <- list()
  }
  if (length(candidatesK) > 1) {
    for (i in 1:length(candidatesK)) {
      if (checkSubsetFrequency(candidatesK[[1]][i,], Lk, k-1) == TRUE) {
        resList <- append(resList, candidatesK[i])
      }
    }
  }
  if (tryCatch({nrow(candidatesK[[1]])},error=function(e){print("0")}) > 1) {
    for (i in 1:nrow(candidatesK[[1]])) {
      if (checkSubsetFrequency(candidatesK[[1]][i,], Lk, k-1) == TRUE) {
        resList <- append(resList, list(candidatesK[[1]][i,]))
      }
    }
  }
  return(resList)
}

#### Main program flow ####
#### Window building ####
windows <- list(data[1,1])
t_end <- as.numeric(rownames(data)[1]) + step
t_start <- t_end - width
noWins <- as.integer((as.numeric(rownames(data)[nrow(data)]) - as.numeric(rownames(data)[1]) + width)/step)
for (i in 1:noWins-1) { #iterating through possible windows
  row <- vector()
  t_start <- t_start + step
  t_end <- t_end + step
  for (j in 1:nrow(data)) {
    if (t_start <= as.numeric(rownames(data)[j]) && as.numeric(rownames(data)[j]) < t_end) {
      row <- c(row, data[j,1])
    }
  }
  windows <- append(windows,list(c(row)))
}
windows[[length(windows)]] <- NULL #delete empty row at the end
#### Combination making ####
C1 <- vector() #first combination
for (i in 1:length(windows)) {
  for (j in 1:length(windows[[i]])) {
    if (is.element(windows[[i]][j], C1) == FALSE) {
      C1 <- c(C1, windows[[i]][j])
    }
  }
}
C1 <- sort(C1)

if (episodeType == "serial") {
  return <- scanWindows_serial(windows, C1) #store C1 data
}
if (episodeType == "parallel") {
  return <- scanWindows_parallel(windows, C1) #store C1 data
}
L1 <- return[[1]]
supportData <- return[[2]]
L <- list(L1)
k <- 2

while (length(L[[k-1]]) > 0) { #next combinations
  if (episodeType == "serial") {
    Ck <- aprioriGen_serial(L[[k-1]], k) #generate combinations and store at Ck
    return <- scanWindows_serial(windows, Ck)
  }
  if (episodeType == "parallel") {
    Ck <- aprioriGen_parallel(L[[k-1]], k) #generate combinations and store at Ck
    return <- scanWindows_parallel(windows, Ck)
  }
  if (length(return) == 1 && return == "done") { break }
  Lk <- return[[1]]
  supK <- return[[2]]
  supportData <- rbind(supportData, supK)
  L <- append(L, list(Lk))
  k <- k + 1
}
#### Rule generation ####
largeItemSet <- L #alias for symmetry with original python code
bigRuleList <- list()
lhsList <- list()
for (i in 2:length(largeItemSet)) { #for each list with 2-element sets or more; i=list number
  for (item in 1:length(largeItemSet[[i]])) { #for each set in a list; item=set number in list
    for (j in 1:(i-1)) { #for each combination of set elements; j=combination size
      itemList <- strsplit(largeItemSet[[i]][[item]], " ")[[1]] #turn set string into list()
      itemString <- largeItemSet[[i]][[item]] #set string alias for easy access
      lhsList <- list(combn(itemList,j)) #create combinations
      if (length(lhsList) > 0) {#if there's at least 1 combination
        if ((nrow(lhsList[[1]]) > 1 && ncol(lhsList[[1]]) > 1)) {
          lhsListLength <- ncol(lhsList[[1]])
        } else {
          lhsListLength <- length(lhsList[[1]])
        }
        for (k in 1:lhsListLength) { #for every individual combination; k=combination number
          if ((nrow(lhsList[[1]]) > 1 && ncol(lhsList[[1]]) > 1)) { #if the combination in itself is length>1
            lhs <- paste(lhsList[[1]][,k], collapse=" ") #parse first
          } else {
            lhs <- lhsList[[1]][k] #get combination directly
          }
          if (!is.na(supportData[itemString,])) { #if the combination exists in the support data
            conf <- (supportData[itemString,])/(supportData[lhs,]) #confidence=set conf./combination conf.
            if (conf >= minConfidence) { #if conf >= minimum confidence
                rule <- paste("[",lhs,"] ==> [",itemString,"] [",width,"] [",supportData[itemString,],", ",conf,"]",sep="")
                bigRuleList <- c(bigRuleList, rule) #parse rule and add to list
            }
          }
        }
      }
    }
  }
}
#### Rule printing ####
print(bigRuleList)

#program by hnievat
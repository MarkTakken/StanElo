args = list(
    utilities_path = "/path/to/DBDA2E-utilities.R",  #Optional, for plotting
    root_dir = "/path/to/root/directory",
    data_file = "example.csv",                    #Must be csv format
    ratings_write_file = "exampleratings.txt",    #Saved as text file if extension is txt, as csv file if extension is csv
    meanElo = 2000,
    maxElo = 5000,                                #The posterior for each rating is uniform from meanElo-maxElo to meanElo+maxElo
    burnin = 300,
    chain.iters = 2500
)

if (!dir.exists(args$root_dir)) {
    stop("Root directory does not exist.  Please edit it in the args list at the top of the source code.")
}

library(comprehenr)
library(hash)
library(rstan)                               #This code assumes the latest version 2.26.1
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#Loading/compiling model
modelString = "
functions {
    real logbig(real power) {
        if (power < 5) return log(1+10^power);
        else return power*log(10);
    }
    real likelihood(int y,int N,real dif) {
        return lchoose(N,y) - N*logbig(dif/400) + (N-y)*dif/400*log(10);
    }
}
data {
    int<lower=1> Npla;
    int<lower=1> Npair;
    real meanElo;
    real sumElo;
    real maxElo;
    array[Npair] int<lower=0> y;
    array[Npair] int<lower=0> N;
    array[Npair,2] int<lower=1,upper=Npla> p;
}
parameters {
    array[Npla-1] real<lower=meanElo-maxElo,upper=meanElo+maxElo> raw;
}
transformed parameters {
    array[Npla] real ratings;
    ratings[1:(Npla-1)] = raw;
    ratings[Npla] = sumElo - sum(raw);
}
model {
    for (i in 1:Npair) {
        target += likelihood(y[i],N[i],ratings[p[i,2]]-ratings[p[i,1]]);
    }
}
"
compiled_model_path = paste(args$root_dir,"/stanDso.RData",sep="")
if (file.exists(compiled_model_path)) {
    load(compiled_model_path)
} else {
    cat("\nCompiling model...\n\n")
    stanDso = stan_model(model_code=modelString)
    save("stanDso",file=compiled_model_path)
}

#Reading data from csv
read_data = function(path) {
    df = read.csv(path)
    df[,1] = trimws(as.character(df[,1]))
    df[,2] = trimws(as.character(df[,2]))
    Npair = nrow(df)
    ind = 1
    h = hash()
    p = matrix(nrow=Npair,ncol=2)
    for (r in 1:Npair) {
        for (c in 1:2) {
            if (is.null(h[[df[r,c]]])) {
                h[[df[r,c]]] = ind
                p[r,c] = ind
                ind = ind + 1
            }
            else {
                p[r,c] = h[[df[r,c]]]
            }
        }
    }
    itostr = c()
    for (n in keys(h)) itostr[h[[n]]] = n
    Npla = length(itostr)
    y = df[,3]
    N = df[,4]
    return (list(y,N,p,Npla,Npair,itostr))
}
data = read_data(paste(args$root_dir,"/",args$data_file,sep=""))
itostr = data[[6]]
Npla = data[[4]]
sumElo = args$meanElo * data[[4]]
dataList = list(y=data[[1]],N=data[[2]],p=data[[3]],Npla=data[[4]],Npair=data[[5]],meanElo=args$meanElo,sumElo=sumElo,maxElo=args$maxElo)

#Running model
stanFit = sampling(object=stanDso,data=dataList,chains=4,iter=args$chain.iters,warmup=args$burnin,thin=1)
print(stanFit)
sampleMat = as.matrix(stanFit)[,Npla:(2*Npla-1)]

#String formatting, writing the ratings file
ext = substring(args$ratings_write_file,gregexpr(pattern="\\.",args$ratings_write_file)[[1]][1]+1)
repstr = function(string,n) {
    if (n == 0) return ("")
    newstr = ""
    for (i in 1:n) newstr = paste(newstr,string,sep="")
    return (newstr)
}
ratings = matrix(nrow=Npla,ncol=5)
for (i in 1:Npla) {
    m = mean(sampleMat[,i])
    low = quantile(sampleMat[,i],0.025)
    high = quantile(sampleMat[,i],0.975)
    ratings[i,1] = round(m)
    ratings[i,2] = round(m - low)
    ratings[i,3] = round(high - m)
    ratings[i,4] = round(low)
    ratings[i,5] = round(high)
}
ord = order(ratings[,1],decreasing=T)
ratings = ratings[ord,]
sampleMat = sampleMat[,ord]
itostr = itostr[ord]
if (ext == "csv") {
    dfo = data.frame(Pla=to_vec(for (ind in 1:Npla) itostr[ind]),mean=ratings[,1],
        minus=ratings[,2],plus=ratings[,3],low=ratings[,4],high=ratings[,5])
    write.csv(dfo,file=paste(args$root_dir,"/",args$ratings_write_file,sep=""))
} else if (ext == "txt") {
    mlcol = c()
    mlcol[1] = max(nchar(itostr))
    for (i in 1:5) {
        mlcol[i+1] = max(to_vec(for (n in ratings[,i]) nchar(n)))
    }
    mlcol[1] = max(mlcol[1],3)
    mlcol[2] = max(mlcol[2],4)
    mlcol[5] = max(mlcol[5],3)
    mlcol[6] = max(mlcol[6],3)
    ratstr = paste("pla ",repstr(" ",mlcol[1]-3),"| mean ",repstr(" ",mlcol[2]-4),
        "| - ",repstr(" ",mlcol[3]-1),"| + ",repstr(" ",mlcol[4]-1),
        "| low ",repstr(" ",mlcol[5]-3),"| high",repstr(" ",mlcol[6]-3),"\n",sep="")
    ratstr = paste(ratstr,repstr("-",nchar(ratstr)),"\n",sep="")
    for (i in 1:Npla) {
        ratstr = paste(ratstr,itostr[i],repstr(" ",mlcol[1]-nchar(itostr[i])+1),sep="")
        for (j in 1:5) {
            ratstr = paste(ratstr,"| ",ratings[i,j],repstr(" ",mlcol[j+1]-nchar(ratings[i,j])+1),sep="")
        }
        ratstr = paste(ratstr,"\n",sep="")
    }
    writeLines(ratstr,con=paste(args$root_dir,"/",args$ratings_write_file,sep=""))
} else {
    cat("Bad file extension\n")
}

#Graphing
plotResults = function() {
    openGraph(x=0,y=0,height=7.9,width=10)
    par(mfrow=c(Npla,Npla))
    for (i in 1:Npla) {
        for (j in 1:Npla) {
            if (i == j) plotPost(sampleMat[,i],cenTend="mean",main=itostr[i],xlab="rating")
            else {
                dif = sampleMat[,i]-sampleMat[,j]
                plotPost(dif,cenTend="mean",main=paste(itostr[i],"-",itostr[j]),xlab="rating")
            }
        }
    }
}
if (Npla < 10 && args$utilities_path != "") {
    source(args$utilities_path)
    plotResults()
}
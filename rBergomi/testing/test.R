data <- read.table(file = "./debug4.txt", header = FALSE)
data.MT <- data.frame(m = data[1:200,1], payoff = data[1:200,2])
data.ST <- data.frame(m = data[201:400,1], payoff = data[201:400,2])
data.MT <- data.MT[order(data.MT$m),]
## Compare
head(data.MT$payoff - data.ST$payoff)
head(data.MT$payoff == data.ST$payoff)
data.MT$payoff == data.ST$payoff
## Count number of wrong results
sum(abs(data.MT$payoff - data.ST$payoff) > 1e-5)

##Same with old code
data <- read.table(file = "./debug5.txt", header = FALSE)
data.MT <- data.frame(m = data[1:200,1], payoff = data[1:200,2])
data.ST <- data.frame(m = data[201:400,1], payoff = data[201:400,2])
data.MT <- data.MT[order(data.MT$m),]
## Compare
head(data.MT$payoff - data.ST$payoff)
head(data.MT$payoff == data.ST$payoff)
data.MT$payoff == data.ST$payoff
## Count number of wrong results
sum(abs(data.MT$payoff - data.ST$payoff) > 1e-5)

##Old code, but explicitly calling omp_get_thread_num() everywhere
data <- read.table(file = "./debug6.txt", header = FALSE)
data.MT <- data.frame(m = data[1:200,1], payoff = data[1:200,2])
data.ST <- data.frame(m = data[201:400,1], payoff = data[201:400,2])
data.MT <- data.MT[order(data.MT$m),]
## Compare
head(data.MT$payoff - data.ST$payoff)
head(data.MT$payoff == data.ST$payoff)
ind <- data.MT$payoff != data.ST$payoff
## Count number of wrong results
sum(abs(data.MT$payoff - data.ST$payoff) > 1e-5)

## New code baed on Convolve class
data.MT <- read.table(file = "./debug7.txt", header = FALSE)
data.MT <- data.frame(m = data.MT[,1], payoff = data.MT[,2])
data.MT <- data.MT[order(data.MT$m),]
data.ST <- read.table(file = "./debug8.txt", header = FALSE)
data.ST <- data.frame(m = data.ST[,1], payoff = data.ST[,2])

data.ST$payoff[data.MT$m+1] - data.MT$payoff
data.MT$payoff[13]
data.ST$payoff[14]

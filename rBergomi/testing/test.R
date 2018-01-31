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
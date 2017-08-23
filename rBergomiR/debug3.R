####################################################
## Check the OMP implementation in the C++ code
####################################################

## Check that normals generated in st and in mt environments
## have right distributions
normal.st <- scan(file = "./Release/normal_ST.txt")
normal.mt <- scan(file = "./Release/normal_MT.txt")

shapiro.test(normal.st)
shapiro.test(normal.mt)
## Neither are rejected
ks.test(normal.st, normal.mt)
## We cannot reject that both are samples from the same distribution.
summary(normal.st)
summary(normal.mt)
qqnorm(normal.st)
qqline(normal.st)
qqnorm(normal.mt)
qqline(normal.mt)

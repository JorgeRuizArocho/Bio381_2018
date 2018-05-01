#### PACKAGE PRESENTATION ####
# AGRICOLAE
# BY Jorge A Ruiz-Arocho
# May 2, 2018

# Preliminaries
install.packages("agricolae")
library(agricolae) 

# Open datasets from 'agricolae' 
data(sweetpotato) 
data(corn)
head(corn) # a numeric vector
str(corn)

# Some basic functions 

x<-seq(10,40,5)
y<-c(2,6,8,7,3,4)

# Poligon and frecuency

h<-graph.freq(x,counts=y,col=colors()[86],xlab=" ", ylab="Frequency",
              axes=FALSE)
axis(1,x,las=2)
axis(2,0:10)
polygon.freq(h,col="red",xlab=" ", ylab="")
title( main="Histogram and polygon", xlab="Variable X")

# Imagine that we have other variables to graph
axis(4,0:10) 

## Experimental designs
# 4 treatments and 5 blocks
trt<-c("A","B","C","D")

# Generating Randomized Complete Block Design (standard design for agricultural experiments. The field or orchard is divided into units to account for any variation in the field. Treatments are then assigned at random to the subjects in the blocks-once in each block.)

outdesign <-design.rcbd(trt,
                        r = 5,serie=2,
                        seed = 45,kinds = "Super-Duper") # Arguments -> r = Replications or blocks, series = number plot,
                                                         # kinds = method for to randomize
outdesign
rcbd <- outdesign$book
rcbd # field book

sketch <- outdesign$sketch
sketch

# Generating Latin Square Design (The Latin square design is used where the researcher desires to control the variation in an experiment that is related to rows and columns in the field.)

varieties<-c("green","blue","white","red")
outdesign <-design.lsd(varieties,serie=2,seed=23)
lsd <- outdesign$book 
print(outdesign$sketch)
print(lsd) # field book.
plots <-as.numeric(lsd[,1])
print(matrix(plots,byrow = TRUE, ncol = 4))


# Write on hard disk.
write.table(lsd,"lsd.txt", row.names=FALSE, sep="\t")
file.show("lsd.txt")

# Several design options are accompanied by complementary analysis functions

Genotype<-paste("geno",1:30,sep="") # the number of treatments must be multiple of k (size block) 
r<- 2 # replication
k<- 3 # block size
plan_alpha <- design.alpha(Genotype,k,r,seed=5)

# Agricolae allows you to conduct a bunch of multiple parametric and non-paprametric comparisons. 
# Parametric comparisons (e.g. Difference test (HSD.test()), the Duncan test (duncan.test()), the Scheffe test (scheffe.test()), among others) 

# Non-Parametric comparisons (e.g. kruskal(), waerden.test(), friedman() and durbin.test()). 
# An important noteto avoid possible confusion: R itself provides the functions kruskal.test() and friedman.test() - these report only a single value; whereas the functions provided by agricolae perform multiple comparisons

data(corn)
kruskal.test(corn$observation, corn$method)

agri_krusk = kruskal(corn$observation, corn$method, group = TRUE, main = "corn") # Main = title
agri_krusk

agri_krusk = kruskal(corn$observation, corn$method, group = FALSE, main = "corn") # Main = title
agri_krusk

# Same issue with Tukey
data("sweetpotato")
head(sweetpotato)

Potato = aov(yield~virus, data= sweetpotato)
summary(Potato)

TukeyHSD(Potato) # default TukeyHSD test on the anova

agricolaeTukey= HSD.test(Potato, "virus", group = TRUE, main = "Yield of sweet poptato with different virus")
agricolaeTukey

agricolaeTukey= HSD.test(Potato, "virus", group = FALSE, main = "Yield of sweet poptato with different virus")
agricolaeTukey

# Graphs for multiple comparisons
data(sweetpotato)
model<-aov(yield~virus,data=sweetpotato)
out <- waller.test(model,"virus", console=TRUE,
                   main="Yield of sweetpotato\ndealt with different virus")

bar.err(out$means,variation="SE",horiz=TRUE,xlim=c(0,45),
        bar=FALSE,col=0)
bar.err(out$means,variation="range",ylim=c(0,45),bar=FALSE,col="grey",
        main="range")

############### Stability analysis
##AMMI

# Significant interaction exists when the first two component terms of a principal component analysis sum to more than 50%. I provide the function plot.AMMI() to visualize the interaction. The parameter type specifies if a biplot or triplot is created.

options(digit=2)

v1 <- c(10.2,8.8,8.8,9.3,9.6,7.2,8.4,9.6,7.9,10,9.3,8.0,10.1,9.4,10.8,
        6.3,7.4)
v2 <- c(7,7.8,7.0,6.9,7,8.3,7.4,6.5,6.8,7.9,7.3,6.8,8.1,7.1,7.1,6.4,
        4.1)
v3 <- c(5.3,4.4,5.3,4.4,5.5,4.6,6.2,6.0,6.5,5.3,5.7,4.4,4.2,5.6,5.8,
        3.9,3.8)
v4 <- c(7.8,5.9,7.3,5.9,7.8,6.3,7.9,7.5,7.6,5.4,5.6,7.8,6.5,8.1,7.5,
        5.0,5.4)
v5 <-c(9,9.2,8.8,10.6,8.3,9.3,9.6,8.8,7.9,9.1,7.7,9.5,9.4,9.4,10.3,
       8.8,8.7)

study <- data.frame(v1, v2, v3, v4, v5)
rownames(study) <- LETTERS[1:17]
study # We assigned letters to each component of each vector

output <- stability.par(study, rep=4, MSerror=2)# This procedure calculates the stability variations as well as the statistics of selection for the yield and the stability. The averages of the genotype through the different environment repetitions are required for the calculations. The mean square error must be calculated from the joint variance analysis.

altitude<-c(1200, 1300, 800, 1600, 2400)
stability <- stability.par(study,rep=4,MSerror=2, cova=TRUE,
                           name.cov= "altitude",
                           file.cov=altitude)

rdto <- c(study[,1], study[,2], study[,3], study[,4], study[,5])
environment <- gl(5,17) 

genotype <- rep(rownames(study),5)
model<-AMMI(ENV=environment, GEN=genotype, REP=4, Y=rdto, MSE=2,
            console=TRUE) # Additive Main Effects and Multiplicative Interaction Models (AMMI) are widely used to analyze main effects and genotype by environment (GEN, ENV) interactions in multilocation variety trials. Furthermore, this function generates data to biplot, triplot graphs and analysis.

plot(model,type=1,las=1)
plot(model,type=2,las=1)

# Another related visualization is the addition of a contour line which adds a contour line proportional to the longest distance of a genotype and has a values between 0 and 1

data(sinRepAmmi)
Environment <- sinRepAmmi$ENV
Genotype <- sinRepAmmi$GEN
Yield <- sinRepAmmi$YLD
REP <- 3
MSerror <- 93.24224
model<-AMMI(Environment, Genotype, REP, Yield, MSerror)

plot(model)
AMMI.contour(model,distance=0.80,shape=8,col="red",lwd=2,lty=5)

#############################################################################################################

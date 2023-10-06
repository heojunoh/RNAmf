# RNAmf

This is the R package `RNAmf` implementing the recursive non-additive emulation for multi-fidelity data 
by incorporating the outputs from previous layer to the input,
which enables a comprehensive analysis of complex computer codes, 
capturing functional correlations between different fidelity levels and enhancing the modeling and prediction capabilities. 
It includes three proposed active learning strategies.
Squared exponential kernel and Matern kernel-1.5, 2.5 are available for covariance structure. 

Details are
> <arXiv:2309.11772>

Maintainer: Junoh Heo <heojunoh@msu.edu>

## Installation

``` r
install.packages("devtools")
library(devtools)
devtools::install_github("heojunoh/RNAmf")
```

## Example
### Nonlinear example function (Perdikaris et al., 2017)

``` r
install_github("cran/MuFiCokriging")
library(MuFiCokriging) # For nested space-filling design
library(Multifidelity) 
library(lhs)
library(ggplot2)
library(hrbrthemes)

### synthetic function ###
f1 <- function(x)
{
  sin(8*pi*x)
}

f2 <- function(x)
{ 
  (x-sqrt(2))*(sin(8*pi*x))^2
}

### training data ###
n1 <- 15; n2 <- 10

set.seed(1)

X1 <- maximinLHS(n1, 1)
X2 <- maximinLHS(n2, 1)

NestDesign <- NestedDesignBuild(design = list(X1,X2))

X1 <- NestDesign$PX
X2 <- ExtractNestDesign(NestDesign,2)

y1 <- f1(X1)
y2 <- f2(X2)

### test data ###
x <- seq(0,1,length.out=100)

### fitting and prediction ###
fit.RNAmf <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
predy <- predRNAmf(fit.RNAmf, x)$mu
predsig2 <- predRNAmf(fit.RNAmf, x)$sig2

### plotting ###
p1 <- ggplot(data.frame(x, predy), aes(x=x), color=group) +
  theme_ipsum() + 
  theme(axis.title.x = element_text(size=20,margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=20,margin = margin(t = 10), hjust=0.5),
        text=element_text(size=16,  family="serif")
  )+
  labs(x="x", y = "y")+ 
  geom_line(aes(y=f2(x)), color="green")+ 
  geom_line(aes(y=f1(x)), color="red")+ 
  geom_line(aes(y=predy), color="blue") + 
  geom_point(data=data.frame(X1,y1),aes(x=X1, y=y1),col="red", shape=2, size=1.5) +
  geom_point(data=data.frame(X2,y2),aes(x=X2, y=y2),col="green", shape=19, size=1.5) +
  geom_line(aes(y=predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2))), color = "blue", linetype = "dotted") + 
  geom_line(aes(y=predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2))), color = "blue", linetype = "dotted") +
  scale_y_continuous(limits = c(-1.5,1.2))
```

![Example plot](https://github.com/heojunoh/Multifidelity/assets/99292199/4c9858d0-28ab-417d-ba9d-05a2c87dafb4)

``` r
print(sqrt(mean((f2(x) - predy)^2)))
[1] 0.07253297
```

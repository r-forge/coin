# conditional inference 05 2005: FRIEDMAN ranks
library(coin)
library(multcomp)
##### complete block design
 RoundingTimes <- data.frame(
         times = c(5.40, 5.50, 5.55,
                   5.85, 5.70, 5.75,
                   5.20, 5.60, 5.50,
                   5.55, 5.50, 5.40,
                   5.90, 5.85, 5.70,
                   5.45, 5.55, 5.60,
                   5.40, 5.40, 5.35,
                   5.45, 5.50, 5.35,
                   5.25, 5.15, 5.00,
                   5.85, 5.80, 5.70,
                   5.25, 5.20, 5.10,
                   5.65, 5.55, 5.45,
                   5.60, 5.35, 5.45,
                   5.05, 5.00, 4.95,
                   5.50, 5.50, 5.40,
                   5.45, 5.55, 5.50,
                   5.55, 5.55, 5.35,
                   5.45, 5.50, 5.55,
                   5.50, 5.45, 5.25,
                   5.65, 5.60, 5.40,
                   5.70, 5.65, 5.55,
                   6.30, 6.30, 6.25),
         methods = factor(rep(c("Round Out", "Narrow Angle", "Wide Angle"), 22)),
         block = factor(rep(1:22, rep(3, 22))))
## H1_unrestricted: friedman blocks     
b111=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999,
              xtrafo = function(data) trafo(data, factor_trafo = function(x) model.matrix(~ x - 1)),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb111=pvalue(b111)
pvalue(b111, adjust=TRUE)
## H1_unrestricted friedman ; was ist ohne xtrafo ??? pb112 nicht=pb111!
b112= oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999,
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb112=pvalue(b112)
pvalue(b112, adjust=TRUE)
## H1_unrestricted friedman blocks quad; was ist Unterschied zu max??
## ist quad ein analogon vom Kruskal-Wallis, aber max ein analogon zu globaltest_Tukey all pairs?
b113=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "quad",distribution="approx",B=99999,
              xtrafo = function(data) trafo(data, factor_trafo = function(x) model.matrix(~ x - 1)),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb113=pvalue(b113)
pvalue(b113, adjust=TRUE)
## H1_unrestricted symmetry test mit friedman blocks; was ist Unterschied zu oneway_test??     
b114=symmetry_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999,
              xtrafo = function(data) trafo(data, factor_trafo = function(x) model.matrix(~ x - 1)),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb114=pvalue(b114)
pvalue(b114, adjust=TRUE)
## H1_unrestricted: als friedman_test; was ist Unterschied zu oneway_test und symmetry_test??     
b115=friedman_test(times ~ methods | block, data = RoundingTimes,distribution="approx",B=99999)
pb115=pvalue(b115)
pvalue(b115, adjust=TRUE)

######## MCP#############################################################
### MCP:all pairs Tukey contrast
b116=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Tukey"))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
              
pb116=pvalue(b116)
pvalue(b116, adjust=TRUE)
## MCP:all pairs Tukey contrast mit quad: ist da unsinn oder was ist Unterschied zu max???
b117=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "quad",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Tukey"))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
              
pb117=pvalue(b117)
pvalue(b117, adjust=TRUE)
## MCP:many-to-one two-sided Dunnett contrast
b118=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Dunnett"))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
              
pb118=pvalue(b118)
pvalue(b118, adjust=TRUE)
## MCP:many-to-one one-sided Dunnett contrast
b119=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Dunnett"))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb119=pvalue(b119)
pvalue(b119, adjust=TRUE) 
# one-sided p-value; so sicher falsch: wie multiple??
     1 - pnorm(sqrt(statistic(b119)))
## MCP:many-to-one one-sided Dunnett contrast mit alternative="less": geht nicht
b120=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, alternative="less",
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Dunnett"))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
 pb120=pvalue(b120)
pvalue(b120, adjust=TRUE)    
## MCP: comparison with mean; NEU
b121=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "AVE"))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb121=pvalue(b121)
pvalue(b121, adjust=TRUE)
## MCP: comparison with best; decomposed in k many-to-one comparisons; geht aber nur mit EINSEITIGEN p-Werten! NEU
b122a=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Dunnett", base=1))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb122a=pvalue(b122a)
pvalue(b122a, adjust=TRUE)
b122b=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Dunnett", base=2))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb122b=pvalue(b122b)
pvalue(b122b, adjust=TRUE)
b122c=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Dunnett", base=3))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb122c=pvalue(b122c)
pvalue(b122c, adjust=TRUE)

###### single contrast tests
## linear c.
clin=c(-1,0,1)
cl=rbind(clin)
b123=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(cl)),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb123=pvalue(b123)
## helmert
chelm=c(-2,1,1)
ch=rbind(chelm)
b124=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(ch)),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb124=pvalue(b124)
# usw

###### multiple contrasts
## williams; NEU
b125=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Williams"))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb125=pvalue(b125)
pvalue(b125, adjust=TRUE)
## marcus; NEU
b126=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Marcus"))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb126=pvalue(b126)
pvalue(b126, adjust=TRUE)
##mcdermott; NEU
b127=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "McDermott"))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb127=pvalue(b127)
pvalue(b127, adjust=TRUE)
## incremental; NEU
b128=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Sequen"))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb128=pvalue(b128)
pvalue(b128, adjust=TRUE)
## hirotsu; NEU
b129=oneway_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",distribution="approx",B=99999, 
              xtrafo = function(data)trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Changepoint"))),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb129=pvalue(b129)
pvalue(b129, adjust=TRUE)

###### user-defined multiple contrasts
## up and down ; NEU    
c1=c(-1,-1,2)
c2=c(-2,1,1)
cges=rbind(c1,c2)
b130=symmetry_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",
              xtrafo = function(data) trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(cges)),
              ytrafo = function(data) trafo(data, numeric_trafo = rank, block = RoundingTimes$block))
pb130=pvalue(b130)
pvalue(b130, adjust=TRUE)

################################################################
## median test: was ist das überhaupt
b131=symmetry_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",
              xtrafo = function(data) trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Tukey"))),
              ytrafo = function(data) trafo(data, numeric_trafo =median_trafo, block = RoundingTimes$block))
pb131=pvalue(b131)
pvalue(b131, adjust=TRUE)
## conover salsburg scores: das gibt es nur für random einweganlage; NEU
b132=symmetry_test(times ~ methods | block, data = RoundingTimes,
              teststat = "max",
              xtrafo = function(data) trafo(data, factor_trafo = function(x)
                      model.matrix(~x-1)%*%t(contrMat(table(x), "Tukey"))),
              ytrafo = function(data) trafo(data, numeric_trafo =consal_trafo, block = RoundingTimes$block))
pb132=pvalue(b132)
pvalue(b132, adjust=TRUE)

###################################################
### chunk number 1: coin-setup
###################################################
options(prompt = "R> ")
library("coin")
library("vcd")
set.seed(290875)


###################################################
### inspect C code and documentation
###################################################
browseURL(system.file("documentation", "html", "index.html", 
                       package = "coin"))


###################################################
### cosmetics on job satisfaction data
###################################################
js <- jobsatisfaction
dimnames(js)[[2]] <- c("VeryDiss", "ModDiss", "ModSat", "VerySat")
ftable(Job.Satisfaction ~ Gender + Income, data = js)


###################################################
### mosaic plot of job satisfaction data
###################################################
cotabplot(js, split_vertical = TRUE, gp = gpar(fill = rev(gray.colors(4))),
          spacing = spacing_highlighting, 
          labeling_args = list(rot_labels = 0, varnames = FALSE, 
                               just_labels = c("center", "right")), 
          panel_args = list(margins = c(3,1,2,3.5)))


###################################################
### conditional Cochran-Mantel-Haenszel test
###################################################
it <- independence_test(js, teststat = "quad", distribution = asymptotic())
it


###################################################
### extract linear statistic
###################################################
statistic(it, "linear")


###################################################
### the same
###################################################
margin.table(js, 1:2)


###################################################
### extract standardized linear statistic
###################################################
statistic(it, "standardized")


###################################################
### apply maximum-type test statistic
###################################################
independence_test(js, teststat = "max")


###################################################
### single-step adjusted p-values
###################################################
pvalue(independence_test(js, teststat = "max"), 
       method = "single-step")


###################################################
### linear-by-linear association test via Monte-Carlo
###################################################
it <- independence_test(js, 
    scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                  Income = c(3, 10, 20, 35)),
    distribution = approximate(B = 10000))
pvalue(it)

###################################################
### chunk number 1: coin-setup
###################################################
options(prompt = "R> ")
library("coin")
library("e1071")
set.seed(290875)
### get rid of the NAMESPACE
attach(asNamespace("coin"))

### use UML? http://argouml.tigris.org/

### extract slots of a class
c2t <- function(x) {        

    classdef <- getClassDef(x)

    extends <- names(classdef@contains)[1]
    if (!is.null(extends)) {
        eslots <- names(getClassDef(extends)@slots)
        slots <- classdef@slots[!names(classdef@slots) %in% eslots]
    } else {
        slots <- classdef@slots
    }

    RET <- cbind(names(slots), slots)
    attr(RET, "contains") <- extends 
    attr(RET, "name") <- x
    class(RET) <- "c2t"        
    RET
}
 
### pretty printing
toLatex.c2t <- function(object, center = TRUE, ...) {        

    RET <- c()

    if (center) RET <- c(RET, "\\begin{center}")

    ### class name
    RET <- c(RET, "\\begin{tabular}{ll}",
                  paste("\\multicolumn{2}{l}{Class \\Rclass{", 
                        attr(object, "name"), "}} \\\\", sep = "")) 

    ### extends?
    if (!is.null(attr(object, "contains")))
        RET <- c(RET, paste("\\multicolumn{2}{l}{Contains \\Rclass{", 
                            attr(object, "contains"), "}} \\\\", sep = ""))

    ### slots
    RET <- c(RET, " & \\\\", "Slot & Class \\\\ \\hline ",
             apply(object, 1, function(x) {
                 x <- cbind(paste("\\code{", x[1], "}", sep = ""),
                            paste("\\Rclass{", x[2], "}", sep = ""))
                 paste(paste(x, collapse = " & "), "\\\\")
             }),
             "\\hline")
    RET <- c(RET, "\\end{tabular}")

    if (center) RET <- c(RET, "\\end{center}")

    class(RET) <- "Latex"
    return(RET)
}


###################################################
### chunk number 2: Ex
###################################################
library("coin")
data("rotarod", package = "coin")
independence_test(time ~ group, data = rotarod, ytrafo = rank)


###################################################
### chunk number 3: IndependenceProblem
###################################################
toLatex(c2t("IndependenceProblem"))


###################################################
### chunk number 4: Ex-IndependenceProblem
###################################################
ip <- new("IndependenceProblem", y = rotarod["time"], x = rotarod["group"])


###################################################
### chunk number 5: IndependenceTestProblem
###################################################
toLatex(c2t("IndependenceTestProblem"))


###################################################
### chunk number 6: Ex-IndependenceTestProblem
###################################################
itp <- new("IndependenceTestProblem", ip)


###################################################
### chunk number 7: IndependenceLinearStatistic
###################################################
toLatex(c2t("IndependenceLinearStatistic"))


###################################################
### chunk number 8: IndependenceTestStatistic
###################################################
toLatex(c2t("IndependenceTestStatistic"))


###################################################
### chunk number 9: ScalarIndependenceTestStatistic
###################################################
toLatex(c2t("ScalarIndependenceTestStatistic"))


###################################################
### chunk number 10: Ex-IndependenceTestStatistic
###################################################
its <- new("IndependenceTestStatistic", itp)


###################################################
### chunk number 11: Ex-IndependenceTestStatistic-statistic
###################################################
statistic(its, "linear")


###################################################
### chunk number 12: Ex-IndependenceTestStatistic-statistic
###################################################
expectation(its)
variance(its)


###################################################
### chunk number 13: Ex-ScalarIndependenceTestStatistic
###################################################
sits <- new("ScalarIndependenceTestStatistic", its, 
             alternative = "two.sided")
statistic(sits, "standardized")


###################################################
### chunk number 14: MaxTypeIndependenceTestStatistic
###################################################
toLatex(c2t("MaxTypeIndependenceTestStatistic"))


###################################################
### chunk number 15: QuadTypeIndependenceTestStatistic
###################################################
toLatex(c2t("QuadTypeIndependenceTestStatistic"))


###################################################
### chunk number 16: PValue
###################################################
toLatex(c2t("PValue"))


###################################################
### chunk number 17: NullDistribution
###################################################
toLatex(c2t("NullDistribution"))


###################################################
### chunk number 18: Ex-NullDistribution-pvalue
###################################################
and <- AsymptNullDistribution(sits)
pvalue(and, statistic(sits))
qperm(and, 0.95)


###################################################
### chunk number 19: IndependenceTest
###################################################
toLatex(c2t("IndependenceTest"))


###################################################
### chunk number 20: IndependenceTest
###################################################
new("IndependenceTest", statistic = sits, distribution = and)


###################################################
### chunk number 21: Ex-distribution
###################################################
set.seed(2908)
tmp <- data.frame(x = rnorm(7), y = rnorm(7))
sexact <- function() {
      x <- object@xtrans  
      y <- object@ytrans  
      perms <- permutations(nrow(x))
      pstats <- apply(perms, 1, function(p) sum(x[p,] * y))
      pstats <- (pstats - expectation(object)) / sqrt(variance(object)) 
      p <- function(q) 1 - mean(pstats > q)
      new("PValue", p = p, pvalue = p)
 }


###################################################
### chunk number 22: Ex-distribution
###################################################
independence_test(y ~ x, data = tmp, alternative = "less",
                   distribution = sexact)


###################################################
### chunk number 23: Ex-wilcox
###################################################
ip <- new("IndependenceProblem", y = rotarod["time"], x = rotarod["group"])
itp <- new("IndependenceTestProblem", ip, 
            ytrafo = function(y) (rank(y) - (nrow(y) + 1) / 2)^2)
its <- new("IndependenceTestStatistic", itp)
sits <- new("ScalarIndependenceTestStatistic", its, 
             alternative = "two.sided")
new("ScalarIndependenceTest", statistic = sits,
     distribution = ExactNullDistribution(sits, algorithm = "split"))


###################################################
### chunk number 24: Ex-wilcox
###################################################
independence_test(time ~ group, data = rotarod, 
                   ytrafo = (rank(y) - (nrow(y) + 1) / 2)^2,
                   distribution = exact())


###################################################
### chunk number 25: js
###################################################
data("jobsatisfaction", package = "coin")
js <- jobsatisfaction
dimnames(js)[[2]] <- c("VeryDiss", "ModDiss", "ModSat", "VerySat")
ftable(Job.Satisfaction ~ Gender + Income, data = js)


###################################################
### chunk number 26: js-plot
###################################################
library("vcd")
cotabplot(js, split_vertical = TRUE, gp = gpar(fill = rev(gray.colors(4))),
          spacing = spacing_highlighting, 
          labeling_args = list(rot_labels = 0, varnames = FALSE, 
                               just_labels = c("center", "right")), 
          panel_args = list(margins = c(3,1,2,3.5)))


###################################################
### chunk number 27: jobsatisfaction-it
###################################################
it <- independence_test(js, teststat = "quad", distribution = asymptotic())
it


###################################################
### chunk number 28: jobsatisfaction-T
###################################################
statistic(it, "linear")


###################################################
### chunk number 29: jobsatisfaction-margin
###################################################
margin.table(js, 1:2)


###################################################
### chunk number 30: jobsatisfaction-stat
###################################################
statistic(it, "standardized")


###################################################
### chunk number 31: jobsatisfaction-max
###################################################
independence_test(js, teststat = "max")


###################################################
### chunk number 32: jobsatisfaction-minp
###################################################
pvalue(independence_test(js, teststat = "max"), 
        method = "single-step")


###################################################
### chunk number 33: jobsatisfaction-ordinal
###################################################
it <- independence_test(js, 
     scores = list(Job.Satisfaction = c(1, 3, 4, 5),
                   Income = c(3, 10, 20, 35)),
     distribution = approximate(B = 10000))
pvalue(it)


###################################################
### chunk number 34: coin-doxygen eval=FALSE
###################################################
## browseURL(system.file("documentation", "html", "index.html", 
##                        package = "coin"))


###################################################
### chunk number 35: motivation-perm-exact
###################################################
pvalue(independence_test(time ~ group, 
                          data = rotarod, distribution = exact()))


###################################################
### chunk number 36: motivation-wmw
###################################################
independence_test(time ~ group, data = rotarod, ytrafo = rank)
rt <- independence_test(time ~ group, data = rotarod, ytrafo = rank, 
                         distribution = exact())
rt


###################################################
### chunk number 37: motivation-wmw-sum
###################################################
statistic(rt, "linear")



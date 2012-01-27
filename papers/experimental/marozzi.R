
service <- data.frame(
time = c(
7.1, 2.2, 3.4, 1.7, 2.2, 1.6, 1.9, 1.8, 1.5, 2.5, 4.5, 5.8, 1.9, 2.2, 1.7,
1.2,
4.6, 1.7, 2.7, 1.9, 2.6, 1.9, 5.2, 1.1, 3.9, 2.8, 2.1, 1.5, 1.6, 3.4, 2.6,
5.2,
3.3, 2.8, 2.1, 1.6, 2.2, 1.5, 1.2, 1.6, 1.6, 3.5, 1.9, 3.4, 1.9, 2.2, 3.1,
1.2,
1.9, 1.5, 2.2, 1.9, 1.5, 2.8, 2.5, 1.5, 10.8, 1.5, 3.1, 2.2, 3.1, 2.3, 3.3,
1.1,
2.6, 4.4, 3.0, 1.7, 2.6, 1.3, 1.8, 3.0, 2.0, 1.9, 1.0,
5.4, 1.3, 2.5, 3.0, 1.4, 2.4, 1.2, 2.6, 1.5, 1.1, 1.2, 1.2, 1.2, 2.2, 1.0,
2.2,
1.8, 2.1, 6.2, 2.1, 3.5, 1.6, 1.5, 1.2, 1.6, 1.4, 1.0, 4.3, 1.4, 4.1, 4.2,
1.5,
1.2, 1.0, 2.8, 2.4, 3.9, 2.2, 1.3, 7.6, 2.3, 1.8, 1.1, 2.1, 1.5, 3.4, 1.4,
1.9,
2.6, 3.0, 1.5, 3.3, 0.7, 1.1, 0.9, 2.0, 1.2, 2.0, 1.2, 1.2, 1.6, 2.9, 5.5,
4.0,
1.4, 3.1, 1.2, 9.6, 1.4, 1.2, 10.7, 1.6, 1.8, 1.2, 1.5),
method = factor(c(rep("first", 16 * 4 + 11), rep("second", 16 * 4 + 11))))



airway <- data.frame(FVC = c(
3.45, 4.00, 4.00, 2.74, 3.95, 4.03, 3.80, 4.05, 4.66, 3.45, 3.49, 4.75,
3.55,
4.14, 3.15, 3.86, 3.85, 4.94, 3.10, 3.65, 4.44, 3.99, 4.13, 4.54, 4.60,
3.73,
3.94, 3.90, 3.20, 3.74, 3.87, 3.44, 4.44, 3.70, 3.10, 4.81, 3.41, 3.38,
3.39,
3.50, 3.62, 4.27, 3.55, 3.82, 4.20, 3.86, 4.34, 4.45, 4.05, 3.60, 4.21,
3.72,
4.73, 3.45, 4.78, 4.54, 3.86, 4.04, 4.46, 3.90, 3.66, 4.08, 3.84, 2.82,
3.24,
3.68, 3.94, 4.10, 4.22, 3.63, 3.42, 4.31, 4.24, 2.92, 4.05, 3.94, 4.10,
4.03, 3.69, 3.83, 3.99, 3.12, 3.43, 3.58, 3.95, 3.78, 3.63, 3.74, 3.84, 3.2,
3.65, 4.29, 4.38, 2.93, 4.77, 4.03, 4.48, 4.26, 3.45, 3.99, 3.78, 2.9, 3.94,
3.84, 3.33, 4.18, 2.70, 3.74, 3.65, 3.72, 4.69, 2.84, 3.34, 3.47, 4.14,
4.78,
4.36, 4.37, 3.20, 3.29, 3.40, 4.40, 3.36, 2.72, 4.21, 3.53, 5.48, 3.62,
3.51,
3.73, 3.40, 3.63, 3.68, 4.07, 3.95, 4.25,
3.04, 4.34, 3.50, 2.68, 3.10, 3.60, 4.93, 3.02, 3.12, 4.05, 4.33, 3.39,
4.24,
4.37, 4.21, 4.87, 4.02, 3.31, 4.25, 4.37, 2.97, 3.89, 3.80, 2.87, 3.89,
4.07,
3.64, 4.62, 4.64, 2.74, 4.34, 4.10, 3.75, 4.06, 3.67, 3.07, 4.59, 3.60),
smoking = factor(c(rep("non", 5 * 13 + 12), rep("light", 4 * 13 + 7),
rep("heavy", 2 * 13 + 12))))


steelball <- data.frame(diameter = c(
1.72, 0.77, 1.62, 1.44, 1.69, 1.29, 0.79, 1.96, 1.79, 0.99,
1.18, 1.09, 1.42, 1.53, 0.69, 1.02, 0.88, 1.19, 1.62, 1.32),
line = factor(c(rep("first", 10), rep("second", 10))))


library("coin")

layout(matrix(1:2, ncol = 2))
hist(subset(service, method = "first")$time)
hist(subset(service, method = "second")$time)

### Example 1: Service Calls
## t-test: 0.2583 OK
t.test(time ~ method, data = service, alt = "greater", var.equal = TRUE)
## T_A: 0.2553 OK
independence_test(time ~ method, data = service, alt = "greater", distribution = "exact")
## T_B: 0.0055 - seems too small
independence_test(time ~ method, data = service, alt = "greater", distribution = "exact", ytrafo =
function(data) trafo(data, numeric_trafo = function(x) as.numeric(x > median(x))))

## T_AB: 0.0101 - seems too small
marozzi <- function(x) cbind(x, x > median(x))
independence_test(time ~ method, data = service, alt = "greater",
                  ytrafo = function(data) trafo(data, numeric_trafo = marozzi))

### Example 2: Airway Function
## F-test: 0.5983 - seems too large
summary(aov(FVC ~ smoking, data = airway))
## T_A: 0.5963 - seems too large
oneway_test(FVC ~ smoking, data = airway, teststat = "quadtype")
## T_B: 0.2070 - seems too large
oneway_test(FVC ~ smoking, data = airway, teststat = "quadtype", 
            ytrafo = function(data) 
                trafo(data, numeric_trafo = function(x) as.numeric(x > median(x))))

## T_AB: 0.2961
independence_test(FVC ~ smoking, data = airway, 
    ytrafo = function(data) trafo(data, numeric_trafo = marozzi))


### Example 3: Steel Balls
## t-test: 0.1056: OK
t.test(diameter ~ line, data = steelball, alt = "greater", var.equal = TRUE)
## T_A: 0.1015: OK
independence_test(diameter ~ line, data = steelball, alt = "greater", distribution = "exact")
## T_B: 0.3243 OK
independence_test(diameter ~ line, data = steelball, alt = "greater", distribution = "exact", 
    ytrafo = function(data) trafo(data, numeric_trafo = function(x) as.numeric(x > median(x))))

## T_AB: 0.1265 maybe OK
independence_test(diameter ~ line, data = steelball, alt = "greater",
                        ytrafo = function(data) trafo(data, numeric_trafo = marozzi))

it <- independence_test(diameter ~ line, data = steelball, alt = "greater",
                        ytrafo = function(data) trafo(data, numeric_trafo =
marozzi), distribution = approximate(B = 200000))
pvalue(it)


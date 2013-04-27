source("hsicForSeminar.R")
set.seed(5002)
n.samples = 400
x1 <- rnorm(n.samples)
y1 <- rnorm(n.samples)
z1 <- x1 + y1
plot(x1,y1)
#dev.copy2pdf(file = "x1_vs_y1.pdf")
show(hsicTest(x1,y1))
show(cor.test(x1,y1))
plot(x1,z1)
#dev.copy2pdf(file = "x1_vs_z1.pdf")
show(hsicTest(x1,z1))
show(cor.test(x1,z1))
a <- runif(n.samples, -pi, pi)
b <- rnorm(n.samples, 1, 0.1)
x2 = cos(a)*b
y2 = sin(a)*b
plot(x2,y2)
#dev.copy2pdf(file = "x2_vs_y2.pdf")
show(hsicTest(x2,y2))
show(cor.test(x2,y2))
# import library
library(tidyverse)

# import data
df = read.delim("TableS5_kovats/ramp5_10.tsv",header = T,stringsAsFactors = F)

# plot ki = f(rt)
p1 <- ggplot(data = df,mapping = aes(x = rt,y=ki)) + 
  geom_point() +
  scale_x_continuous(limits=c(0,2000)) +
  scale_y_continuous(limits=c(0,2000))
p1

# fit a non-linear model using the nls function (non-linear weighted least-squares estimates of the parameters of a non-linear model)
# from the graph: a_start=600
#a_start = 600 # a when b = 0
#b_start = -2*log(a_start)    # b when a = 0
#model = nlsLM(ki ~ exp(a + b*rt),data = df,start=list(a=a_start,b=10))

###################################
# simulation of an exponential curve
# y = a * exp(bx)
# beta <- 0.05
# n <- 100
# temp <- data.frame(y = exp(beta * seq(n)) + rnorm(n), x = seq(n))
# mod = nls(y ~ exp(a + b*x),data = temp,start=list(a=0,b=0))
# plot(temp$x,temp$y)
# lines(temp$x,predict(mod,temp$x),type="l",col="red",lty=2,lwd=1)

######################
# simulation of a second degree polynomial curve
# n <- 100
# x <- seq(n)
# y <- rnorm(n, 50 + 30* x + 30*x^2, 1)
# temp <- data.frame(x, y)
# plot(y~x,temp)
# 
# mod = nls(y ~ c + a*x + b*x^2,data = temp,start = list(a=0,b=0,c=600))
# lines(temp$x,predict(mod,temp$x),type="l",col="red",lty=2,lwd=1)

# on our RT/KI
y = df$ki
x = df$rt
temp <- data.frame(x, y)
plot(x,y,xlim = c(0,2000),ylim = c(0,2000))
mod = nls(y ~ c + a*x + b*x^2,data = temp,start = list(a=0,b=0,c=600))

# R regression coefficient 
coef_cor = cor(y,predict(mod,x))
coef_cor = coef_cor^2
to_print1 = paste("The squared R coefficient is:",round(coef_cor,3))

# extract coefficients
a = round(summary(mod)$coefficients[1],digits = 2)
b = round(summary(mod)$coefficients[2],digits=5)
c = round(summary(mod)$coefficients[3],digits=2)
to_print2 = paste("y = ",a,".x + ",b,".x^2 +",c,sep = '')

# plot
p2 = p1 + geom_line(aes(x = rt,y=predict(mod,x)),colour="red")
p2 = p2 + annotate("text",x = 500,y=1800,label=to_print1)
p2 = p2 + annotate("text",x = 500,y=1500,label=to_print2) +
  labs(y="Theoretical Kovats Index",x="Retention Time (seconds)")
print(p2)

ggsave(plot = p2,filename = "TableS5_kovats/RT_vs_KI.pdf",width = 7,height = 5)



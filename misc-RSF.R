
rsf.2.df <- data.frame(y = values(mtnlion.count.rast), 
                       elev = scale(values(elevation.rast)), 
                       slope = scale(values(slope.rast)), 
                       exposure = scale(values(exposure.rast)))

poisson.reg <- glm(y ~ elev + slope + exposure,
                   family = poisson(link = "log"),
                   data = rsf.2.df)

```{r}
library(mvtnorm)

x = seq(-4, 4, by = 0.1)
y = seq(-4, 4, by = 0.1)
xy = expand.grid(x,y)
z = dmvnorm(xy, mean = c(0,0), sigma = diag(2))
persp(x,y,matrix(z * 20, length(x), length(y)), 
      scale = FALSE, 
      box = FALSE, 
      theta = 30, phi = 40)

```

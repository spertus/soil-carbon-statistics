# simulating correlated random fields (plots) using gstat package

library(gstat)

x <- 1:200
y <- 1:200
xy <- expand.grid(x, y)
colnames(xy) <- c("x", "y")

vario <- vgm(nugget = .01, psill = .05, range = 20, model = "Exp")
gstat_mod <- gstat(formula = z ~ 1 + x + y, locations = ~ x + y, dummy = TRUE, beta = c(-2, 0, .01), model = vario, nmax = 20) 

simulation <- predict(gstat_mod, newdata = xy, nsim = 1)

#transform to bounded (as proportions must be)
simulation <- simulation %>% 
  mutate(bounded_sim = pnorm(sim1))

sim_plot <- ggplot(data = simulation, aes(x = x, y = y, fill = bounded_sim)) + 
  geom_raster()
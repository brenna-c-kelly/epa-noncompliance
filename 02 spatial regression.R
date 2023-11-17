
library(INLA)
library(spdep)
library(spatialreg)


caa_contiguous <- caa_eval[which(!caa_eval$fac_state %in% c("AK", "HI",
                                                             "PR", "VI")), ]

xy <- data.matrix(caa_contiguous[, c("fac_long", "fac_lat")])

system.time(
  nb.gab <- spdep::graph2nb(spdep::relativeneigh(xy), sym = TRUE))
par(mar = c(0,0,0,0))
plot(nb.gab, xy)

#listw.gab <- spdep::nb2listw(nb.gab, zero.policy = TRUE)

nb2INLA("graph.adj", nb.gab)
g <- inla.read.graph(filename = "graph.adj")

caa_contiguous$id_area <- 1:nrow(caa_contiguous)
caa_contiguous$id_area_2 <- 1:nrow(caa_contiguous)

table(caa_contiguous$caa_violations == 0)

x_his_z_hi <- c(x_his_hi, nothing1)


caa_contiguous$caa_violations_bin <- ifelse(caa_contiguous$caa_violations == 0,
                                            0, 1)

n = nrow(caa_contiguous)

# outcome matrix
nothing1 <- rep(NA, n)
nothing2 <- rep(NA, n)
summary(nb)

b = as.vector(caa_contiguous$caa_violations_bin == 1)
nb = ifelse(b == 1, caa_contiguous$caa_violations, NA)

bNA = as.vector(c(b, nothing1))
nbNA = as.vector(c(nothing2, nb))

outcome.matrix <- matrix(c(bNA, nbNA), ncol = 2)
summary(outcome.matrix)

# intercept vectors
mu_b <- c(rep(1, n), nothing1) # Binomial 
mu_nb <- c(nothing2, rep(1,  n)) # Gamma 

# spatial error terms (for structured and unstructured error)
id_space <- caa_contiguous$id_area

id_space_b <- c(id_space, nothing1) # Binomial
id_space_nb <- c(nothing2, id_space) # Gamma
id_space_b2 <- c(id_space, nothing1) # Binomial
id_space_nb2 <- c(nothing2, id_space) # Gamma

# covariates
# use sf, assign census tract to facility; population, race, poverty, MHHI

table(caa_contiguous$industry_sector)
mean(caa_contiguous$caa_violations)
aggregate(caa_contiguous$caa_violations, by = list(caa_contiguous$industry_sector), FUN = mean)

# industry
ind_constr <- ifelse(caa_contiguous$industry_sector == "construction", 1, 0)
ind_edh <- ifelse(caa_contiguous$industry_sector == "education and health services", 1, 0)
ind_fin <- ifelse(caa_contiguous$industry_sector == "financial activities", 1, 0)
ind_info <- ifelse(caa_contiguous$industry_sector == "information", 1, 0)
ind_hosp <- ifelse(caa_contiguous$industry_sector == "leisure and hospitality", 1, 0)
ind_manu <- ifelse(caa_contiguous$industry_sector == "manufacturing", 1, 0)
ind_mine <- ifelse(caa_contiguous$industry_sector == "natural resources and mining", 1, 0)
ind_na <- ifelse(caa_contiguous$industry_sector == "nonclassifiable", 1, 0)
ind_other <- ifelse(caa_contiguous$industry_sector == "other services (except public admin)", 1, 0)
ind_prof <- ifelse(caa_contiguous$industry_sector == "professional and business services", 1, 0)
ind_pub <- ifelse(caa_contiguous$industry_sector == "public administration", 1, 0)
ind_trans <- ifelse(caa_contiguous$industry_sector == "trade, transportation and utilities", 1, 0)

ind_constr_b <- c(ind_constr, nothing1)
ind_constr_nb <- c(nothing2, ind_constr)
ind_edh_b <- c(ind_edh, nothing1)
ind_edh_nb <- c(nothing2, ind_edh)
ind_fin_b <- c(ind_fin, nothing1)
ind_fin_nb <- c(nothing2, ind_fin)
ind_info_b <- c(ind_info, nothing1)
ind_info_nb <- c(nothing2, ind_info)
ind_hosp_b <- c(ind_hosp, nothing1)
ind_hosp_nb <- c(nothing2, ind_hosp)
ind_manu_b <- c(ind_manu, nothing1)
ind_manu_nb <- c(nothing2, ind_manu)
ind_mine_b <- c(ind_mine, nothing1)
ind_mine_nb <- c(nothing2, ind_mine)
ind_na_b <- c(ind_na, nothing1)
ind_na_nb <- c(nothing2, ind_na)
ind_other_b <- c(ind_other, nothing1)
ind_other_nb <- c(nothing2, ind_other)
ind_prof_b <- c(ind_prof, nothing1) # referent
ind_prof_nb <- c(nothing2, ind_prof) # referent
ind_pub_b <- c(ind_pub, nothing1)
ind_pub_nb <- c(nothing2, ind_pub)
ind_trans_b <- c(ind_trans, nothing1)
ind_trans_nb <- c(nothing2, ind_trans)

ind_b <- c(caa_contiguous$industry_sector, nothing1)
ind_p <- c(nothing2, caa_contiguous$industry_sector)

# state
table(caa_contiguous$fac_state)

state_b <- c(caa_contiguous$fac_state, nothing1)
state_p <- c(nothing2, caa_contiguous$fac_state)


# put it all together
hurdle_dat <- list(outcome.matrix = outcome.matrix, 
                   id_space_b = id_space_b, id_space_nb = id_space_nb,
                   id_space_b2 = id_space_b2, id_space_nb2 = id_space_nb2,
                   state_b = state_b, state_p = state_p,
                   mu_b = mu_b, mu_nb = mu_nb,
                   ind_constr_b = ind_constr_b,
                   ind_constr_nb = ind_constr_nb,
                   ind_edh_b = ind_edh_b,
                   ind_edh_nb = ind_edh_nb,
                   ind_fin_b = ind_fin_b,
                   ind_fin_nb = ind_fin_nb,
                   ind_info_b = ind_info_b,
                   ind_info_nb = ind_info_nb,
                   ind_hosp_b = ind_hosp_b,
                   ind_hosp_nb = ind_hosp_nb,
                   ind_manu_b = ind_manu_b,
                   ind_manu_nb = ind_manu_nb,
                   ind_mine_b = ind_mine_b,
                   ind_mine_nb = ind_mine_nb,
                   ind_na_b = ind_na_b,
                   ind_na_nb = ind_na_nb,
                   ind_other_b = ind_other_b,
                   ind_other_nb = ind_other_nb,
                   ind_prof_b = ind_prof_b,
                   ind_prof_nb = ind_prof_nb,
                   ind_pub_b = ind_pub_b,
                   ind_pub_nb = ind_pub_nb,
                   ind_trans_b = ind_trans_b,
                   ind_trans_nb = ind_trans_nb,
                   ind_b = ind_b,
                   ind_p = ind_p)


formula <- outcome.matrix ~ mu_b + mu_nb + 
  #f(state_b, model = "iid") + f(state_p, model = "iid") + #industry_sector + #f(fac_state) +
  # ind_constr_b + ind_constr_nb + ind_edh_b + ind_edh_nb +
  # ind_fin_b + ind_fin_nb + ind_info_b + ind_info_nb +
  # ind_hosp_b + ind_hosp_nb + ind_manu_b + ind_manu_nb +
  # ind_mine_b + ind_mine_nb + ind_na_b + ind_na_nb +
  # ind_other_b + ind_other_nb + ind_prof_b + ind_prof_nb +
  # ind_trans_b + ind_trans_nb + #ind_pub_b+ ind_pub_nb + 
  f(ind_b, model = "iid") + f(ind_p, model = "iid") +
  f(id_space_b, model = "bym2", graph = g) +
  f(id_space_b2, model = "iid") +
  f(id_space_nb, model = "bym2", graph = g) +
  f(id_space_nb2, model = "iid") - 1

system.time(res <- inla(formula,
            data = hurdle_dat, #caa_contiguous,
            #E = E, 
            family = c("binomial", "poisson"),
            control.inla = list(int.strategy = "eb"),
            control.predictor = list(compute = TRUE),
            control.compute = list(dic = TRUE, waic = TRUE))
)
summary(res)

# re for ind, no spatial aic: 49100.46, dic: 54276.68

#local.plot.result(res)

res$summary.random

results <- round(res$summary.fixed, 2)
results$cred <- case_when(results$`0.025quant` > 0 & results$`0.975quant` > 0 ~ "+",
                          results$`0.025quant` < 0 & results$`0.975quant` < 0 ~ "-",
                          results$`0.025quant` >= 0 & results$`0.975quant` <= 0 ~ ".",
                          results$`0.025quant` <= 0 & results$`0.975quant` >= 0 ~ ".")
results

plot(res, plot.fixed.effects = FALSE,
     plot.random.effects = FALSE,
     plot.hyperparameters = TRUE,
     plot.predictor = FALSE, cex = 1.25)

table(caa_contiguous$caa_days_last_evaluation,
      caa_contiguous$fac_state)

plot(res, plot.fixed.effects = FALSE,
     plot.random.effects = TRUE,
     plot.hyperparameters = FALSE,
     plot.predictor = FALSE, cex = 1.25)

# just space, intercept   aic: 39648.01 nb   dic: 40088.13
# space, int, industry    aic: 39365.43 nb   dic: 39784.82
# space, int, industry    aic: 39359.55 pois dic: 39784.92

length(res$residuals$deviance.residuals)
res$marginals.random$id_space_b


summary(res$summary.random$id_space_b)

# hurdle_dat$re_b <- res$summary.random$id_space_b

caa_contiguous$re_b <- res$summary.random$id_space_b[1:37775, 2]
caa_contiguous$re_nb <- res$summary.random$id_space_nb[1:37775, 2]

z_st$re_st_b <- res$summary.random$state_b[["mean"]]
z_st$re_st_p <- res$summary.random$state_p[["mean"]]

tm_shape(z_st, proj = aea) +
  tm_polygons(col = "re_st_b", style = "cont", palette = "-inferno", lwd = 0)
tm_shape(z_st, proj = aea) +
  tm_polygons(col = "re_st_p", style = "cont", palette = "-inferno", lwd = 0)


caa_contiguous$y_hat_b <- res$summary.linear.predictor$mean[1:37775]
ufo_inla_rr$summary.fitted.values[, "mean"]
caa_contiguous$fitted <- res$summary.fitted.values$mean[1:37775]
caa_contiguous$fitted_p <- res$summary.fitted.values$mean[37776:75550]

res$summary.random


caa_sf <- st_as_sf(caa_contiguous, coords = c("fac_long", "fac_lat"),
                   crs = 4326, agr = "constant")

library(tmap)
aea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +ellps=GRS80 +datum=NAD83"

z_st <- get_acs(geography = "state",
                variables = c("B01001_001"), # total population
                #state = i,
                geometry = TRUE,
                year = 2021)
#z_st <- st_as_sf(z_st, crs = 4326, agr = "constant")
z_st <- st_transform(z_st, crs = aea)
z_st <- z_st[which(!z_st$NAME %in% c("Alaska", "Hawaii", "Puerto Rico")), ]

tm_shape(z_st, proj = aea) +
  tm_polygons(col = "white", alpha = 0.5, lwd = 0.5) +
  tm_shape(caa_sf, proj = aea) +
  tm_dots(col = "y_hat_b_t", style = "cont", palette = "-inferno", alpha = 0.9,
          breaks = seq(0, 1, by = 0.2)) +
  tm_layout(#aes.palette = list(seq = "-RdYlGn"), 
    title = "y hat, binomial", 
    title.position = c(0.6, 0.95))
caa_sf$y_hat_b_t <- exp(caa_sf$y_hat_b) / (1 + exp(caa_sf$y_hat_b))


tm_shape(z_st, proj = aea) +
  tm_polygons(col = "white", alpha = 0.5, lwd = 0.5) +
  tm_shape(caa_sf, proj = aea) +
  tm_dots(col = "re_b", style = "cont", palette = "-inferno", alpha = 0.9) +
  tm_layout(#aes.palette = list(seq = "-RdYlGn"), 
            title = "Spatial RE, binomial", 
            title.position = c(0.6, 0.95))
tm_shape(z_st, proj = aea) +
  tm_polygons(col = "white", alpha = 0.5, lwd = 0.5) +
  tm_shape(caa_sf, proj = aea) +
  tm_dots(col = "re_nb", style = "cont", palette = "seq", alpha = 0.9) +
  tm_layout(aes.palette = list(seq = "-RdYlGn"), title = "Spatial RE, poisson", 
            title.position = c(0.6, 0.95))

tm_shape(z_st, proj = aea) +
  tm_polygons(col = "white", alpha = 0.5, lwd = 0.5) +
  tm_shape(caa_sf, proj = aea) +
  tm_dots(col = "fitted", style = "cont", palette = "-inferno", alpha = 0.9,
          breaks = seq(0, 1, by = 0.2)) +
  tm_layout(#aes.palette = list(seq = "-RdYlGn"), 
            title = "Fitted binomial", 
            title.position = c(0.65, 0.95))
tm_shape(z_st, proj = aea) +
  tm_polygons(col = "white", alpha = 0.5, lwd = 0.5) +
  tm_shape(caa_sf, proj = aea) +
  tm_dots(col = "fitted_p_log", style = "cont", palette = "-inferno", alpha = 0.9) + #, 
          #breaks = seq(0, 1, by = 0.1)) +
  tm_layout(#aes.palette = list(seq = "-RdYlGn"), 
            title = "Fitted poisson", 
            title.position = c(0.65, 0.95))

caa_sf$fitted_p_log <- log(caa_sf$fitted_p + 1)
summary(caa_sf$fitted_p_exp)

plot(hurdle_dat$re_b)

summary(res$summary.random$id_space_nb)


dlist <- spdep::nbdists(nb.gab, xy)
listw.d1 <- spdep::nb2listw(nb.gab, style = "W", zero.policy = TRUE, glist=dlist)
dlist <- lapply(dlist, function(x) 1/x)
dlist <- lapply(dlist, function(x) 1/x^2)
listw.d2 <- spdep::nb2listw(nb.gab, style = "W", zero.policy = TRUE, glist=dlist)

#caa_contiguous$inspected <- ifelse(is.na(caa_contiguous$inspected), 0, 
#                                    caa_contiguous$inspected)

spdep::moran.test(caa_contiguous$caa_violations, listw.gab, zero.policy = TRUE) # neighbors
spdep::moran.test(caa_contiguous$caa_violations, listw.d1, zero.policy = TRUE) # inverse distance weights
spdep::moran.test(caa_contiguous$caa_violations, listw.d2, zero.policy = TRUE) # inverse squared distance weights

#model.lm <- nlme::gls(caa_violations ~ 1, data = caa_contiguous, method = "REML")
#semivario <- nlme::Variogram(model.lm)#, form = ~x  + y, resType = "normalized")

ggplot(data = semivario, aes(x = dist, y = variog)) + 
  geom_point() + 
  geom_smooth(se=FALSE) +
  geom_hline(yintercept=1) + 
  ylim(c(0, 2.5)) + 
  xlab("Distance") + 
  ylab("Semivariance")

inla()


glmbase <- glm(c(caa_eval$inspected) ~ 1, family = "binomial")
names(caa_contiguous)
system.time(errorsarlm(caa_violations ~ 1, data = caa_contiguous, listw = listw.gab,
                       etype = "error", method="eigen", Durbin = FALSE)
)
summary(listw.gab$neighbours)
require("spdep", quietly=TRUE)
data(hopkins, package="spData")
hopkins_part <- hopkins[21:36,36:21]
hopkins_part[which(hopkins_part > 0, arr.ind=TRUE)] <- 1
hopkins.rook.nb <- spdep::cell2nb(16, 16, type="rook")
glmbase <- glm(c(hopkins_part) ~ 1, family="binomial")
lw <- spdep::nb2listw(hopkins.rook.nb, style="B")
set.seed(123)
system.time(MEbinom1 <- ME(c(hopkins_part) ~ 1, family="binomial",
                           listw=lw, alpha=0.05, verbose=TRUE, nsim=49))
summary(MEbinom1)

glmME <- glm(c(hopkins_part) ~ 1 + fitted(MEbinom1), family="binomial")
#anova(glmME, test="Chisq")
coef(summary(glmME))

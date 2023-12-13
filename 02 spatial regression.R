
library(INLA)
library(spdep)
library(spatialreg)
library(RColorBrewer)

caa_eval <- read.csv("caa_eval.csv")

caa_contiguous <- caa_eval[which(!caa_eval$fac_state %in% c("AK", "HI",
                                                            "PR", "VI")), ]
# check <- check[which(!check$fac_state %in% c("AK", "HI",
#                                                       "PR", "VI")), ]
# 
# caa_contiguous <- merge(caa_eval, check, by = c("registry_id", "fac_name")) |>
#   filter(!fac_state.x %in% c("AK", "HI", "PR", "VI"))

xy <- data.matrix(caa_contiguous[, c("fac_long", "fac_lat")])

system.time(
  nb.gab <- spdep::graph2nb(spdep::relativeneigh(xy), sym = TRUE))
par(mar = c(0,0,0,0))
plot(nb.gab, xy)

#listw.gab <- spdep::nb2listw(nb.gab, zero.policy = TRUE)

nb2INLA("graph.adj", nb.gab)
g <- inla.read.graph(filename = "graph.adj")

caa_contiguous$id_area <- 1:nrow(caa_contiguous)

caa_contiguous$caa_violations_bin <- ifelse(caa_contiguous$caa_violations == 0,
                                            0, 1)

n = nrow(caa_contiguous)

# outcome matrix
nothing1 <- rep(NA, n)
nothing2 <- rep(NA, n)

b = as.vector(caa_contiguous$caa_violations_bin == 1) #ifelse(caa_contiguous$caa_violations == 0, 0, 1) #
p = ifelse(b == 1, caa_contiguous$caa_violations, NA)

bNA = as.vector(c(b, nothing1))
pNA = as.vector(c(nothing2, p))

outcome.matrix <- matrix(c(bNA, pNA), ncol = 2)

# intercept vectors
mu_b <- c(rep(1, n), nothing1) # Binomial 
mu_p <- c(nothing2, rep(1,  n)) # pois 

# spatial error terms (for structured and unstructured error)
id_space <- 1:nrow(caa_contiguous)

id_space_b <- c(id_space, nothing1) # Binomial
id_space_p <- c(nothing2, id_space) # Gamma
id_space_b2 <- c(id_space, nothing1) # Binomial
id_space_p2 <- c(nothing2, id_space) # Gamma

# covariates
# use sf, assign census tract to facility; population, race, poverty, MHHI
typeof(caa_contiguous$sector_naics)

# sector
caa_contiguous$sector_naics <- ifelse(caa_contiguous$sector_naics %in% 
                                        c("31", "32", "33"), "31", 
                                      caa_contiguous$sector_naics) # manufacturing
caa_contiguous$sector_naics <- ifelse(caa_contiguous$sector_naics %in% 
                                        c("44", "45"), "44", 
                                      caa_contiguous$sector_naics) # retail trade
caa_contiguous$sector_naics <- ifelse(caa_contiguous$sector_naics %in% 
                                        c("48", "49"), "48", 
                                      caa_contiguous$sector_naics) # warehousing
sec_b <- c(caa_contiguous$sector_naics, nothing1)
sec_p <- c(nothing2, caa_contiguous$sector_naics)

# industry
ind_b <- c(caa_contiguous$industry_naics, nothing1)
ind_p <- c(nothing2, caa_contiguous$industry_naics)

# state
state_b <- c(caa_contiguous$fac_state, nothing1)
state_p <- c(nothing2, caa_contiguous$fac_state)

# put it all together
hurdle_dat <- list(outcome.matrix = outcome.matrix, 
                   mu_b = mu_b, mu_p = mu_p, 
                   id_space_b = id_space_b, id_space_p = id_space_p,
                   id_space_b2 = id_space_b2, id_space_p2 = id_space_p2,
                   sec_b = sec_b,  sec_p = sec_p, 
                   ind_b = ind_b,  ind_p = ind_p,
                   state_b = state_b, state_p = state_p,
                   sec_b_2 = sec_b,  sec_p_2 = sec_p)


formula <- outcome.matrix ~ mu_b + mu_p + 
  f(state_b) + f(state_p) + #industry_sector + #f(fac_state) +, model = "iid")
  f(sec_b) + f(sec_p) +
  f(id_space_b, model = "bym2", graph = g) +
  f(id_space_b2, model = "iid") + - 1
  f(id_space_p, model = "bym2", graph = g) +
  f(id_space_p2, model = "iid") - 1

system.time(res <- inla(formula,
            data = hurdle_dat, #caa_contiguous,
            #E = E, 
            family = c("binomial", "poisson"),
            control.inla = list(int.strategy = "eb"),
            control.predictor = list(link = 1, compute = TRUE),
            control.compute = list(dic = TRUE, waic = TRUE)))

summary(res)
# re: state, space, sector:   39168
# re: space, sector:          39260.32
# re: space, sector:          45105.52
# re: sp, sec fe: sec         45249.93
# re: state, space, sector    42083.64

caa_contiguous$residual_b <- res$residuals$deviance.residuals[1:37775]
caa_contiguous$residual_p <- res$residuals$deviance.residuals[37776:75550]

caa_contiguous_sf <- st_as_sf(caa_contiguous, coords = c("fac_long", "fac_lat"),
                              crs = 4326, agr = "constant")
listw.gab <- spdep::nb2listw(nb.gab, zero.policy = TRUE)

spdep::moran.test(caa_contiguous_sf$residual_b, listw.gab, zero.policy = TRUE)

caa_contiguous_sf_pois <- caa_contiguous_sf[which(caa_contiguous_sf$caa_violations > 0), ]
caa_contiguous_sf_pois$caa_violations_log <- log(caa_contiguous_sf_pois$caa_violations)
summary(caa_contiguous_sf_pois$caa_violations)
var(caa_contiguous_sf_pois$caa_violations)

prop.table(table(caa_contiguous_sf_pois$caa_violations == 1))
prop.table(table(caa_contiguous_sf_pois$caa_violations == 2))
prop.table(table(caa_contiguous_sf_pois$caa_violations == 3))
prop.table(table(caa_contiguous_sf_pois$caa_violations == 4))
prop.table(table(caa_contiguous_sf_pois$caa_violations >= 5))

hist(caa_contiguous_sf_pois$caa_violations_log)
caa_contiguous_sf_pois[which(caa_contiguous_sf_pois$caa_violations > 100), ]

tm_shape(states, proj = aea) +
  tm_polygons(col = "white", border.col = "gray10") +
  tm_shape(caa_contiguous_sf_pois) +
  tm_dots(col = "caa_violations_log", style = "cont", palette = "Reds", legend.show = FALSE)

# re for ind, no spatial waic: 49100.46, dic: 54276.68
# re for ind, space      waic: 39348.86, dic: 39775.68

#local.plot.result(res)

predict()

aggregate(caa_contiguous$caa_violations > 0, by = list(caa_contiguous$industry_sector), FUN = mean)
aggregate(caa_contiguous$caa_violations, by = list(caa_contiguous$industry_sector), FUN = mean)

# results <- round(res$summary.fixed, 2)

results$cred <- case_when(results$`0.025quant` > 0 & results$`0.975quant` > 0 ~ "+",
                          results$`0.025quant` < 0 & results$`0.975quant` < 0 ~ "-",
                          results$`0.025quant` >= 0 & results$`0.975quant` <= 0 ~ ".",
                          results$`0.025quant` <= 0 & results$`0.975quant` >= 0 ~ ".")


sector_bin <- cbind(res$summary.random$sec_b$ID, 
                    round(res$summary.random$sec_b[, c("mean", "sd",
                                                       "0.025quant", "0.5quant",
                                                       "0.975quant", "mode")], 2))
sector_bin$cred <- case_when(sector_bin$`0.025quant` > 0 & sector_bin$`0.975quant` > 0 ~ "+",
                             sector_bin$`0.025quant` < 0 & sector_bin$`0.975quant` < 0 ~ "-",
                             sector_bin$`0.025quant` >= 0 & sector_bin$`0.975quant` <= 0 ~ ".",
                             sector_bin$`0.025quant` <= 0 & sector_bin$`0.975quant` >= 0 ~ ".")

sector_poi <- cbind(res$summary.random$sec_p$ID, 
                    round(res$summary.random$sec_p[, c("mean", "sd",
                                                       "0.025quant", "0.5quant",
                                                       "0.975quant", "mode")], 2))
sector_poi$cred <- case_when(sector_poi$`0.025quant` > 0 & sector_poi$`0.975quant` > 0 ~ "+",
                             sector_poi$`0.025quant` < 0 & sector_poi$`0.975quant` < 0 ~ "-",
                             sector_poi$`0.025quant` >= 0 & sector_poi$`0.975quant` <= 0 ~ ".",
                             sector_poi$`0.025quant` <= 0 & sector_poi$`0.975quant` >= 0 ~ ".")

creds_bin <- data.frame()
creds_poi <- data.frame()

for(i in 1:21) {
  # binomial
  if(sector_bin[i, "mean"] < 0) {
    h_cred <- inla.pmarginal(0.99, res$marginals.random$sec_b[[i]])
    m_cred <- inla.pmarginal(0.89, res$marginals.random$sec_b[[i]])
    w_cred <- inla.pmarginal(0.67, res$marginals.random$sec_b[[i]])
  }
  else{
    h_cred <- 1 - inla.pmarginal(0.99, res$marginals.random$sec_b[[i]])
    m_cred <- 1 - inla.pmarginal(0.89, res$marginals.random$sec_b[[i]])
    w_cred <- 1 - inla.pmarginal(0.67, res$marginals.random$sec_b[[i]])
  }
  creds_bin <- rbind(creds_bin, c(res$summary.random$sec_b$ID[i],
                                  round(h_cred, 3), 
                                  round(m_cred, 3), 
                                  round(w_cred, 3)))
  # poisson
  if(sector_poi[i, "mean"] < 0) {
    h_cred <- inla.pmarginal(0.99, res$marginals.random$sec_p[[i]])
    m_cred <- inla.pmarginal(0.89, res$marginals.random$sec_p[[i]])
    w_cred <- inla.pmarginal(0.67, res$marginals.random$sec_p[[i]])
  }
  else{
    h_cred <- 1 - inla.pmarginal(0.99, res$marginals.random$sec_p[[i]])
    m_cred <- 1 - inla.pmarginal(0.89, res$marginals.random$sec_p[[i]])
    w_cred <- 1 - inla.pmarginal(0.67, res$marginals.random$sec_p[[i]])
  }
  creds_poi <- rbind(creds_poi, c(res$summary.random$sec_p$ID[i],
                                  round(h_cred, 3), 
                                  round(m_cred, 3), 
                                  round(w_cred, 3)))
}

names(creds_bin) <- c("sector", "high", "moderate", "weak")
names(creds_poi) <- c("sector", "high", "moderate", "weak")

creds_bin
creds_poi

sector_bin <- sector_bin |>
  mutate_at(c("mean", "sd",
              "0.025quant", "0.5quant",
              "0.975quant", "mode"), ~ round(plogis(.), 2))
sector_poi <- sector_poi |>
  mutate_at(c("mean", "sd",
              "0.025quant", "0.5quant",
              "0.975quant", "mode"), ~ round(plogis(.), 2))


# checking for over/underdispersion
sqrt(var(res$summary.fitted.values$mean[1:(37775 + 1)]))
sqrt(var(caa_contiguous$caa_violations_bin))

sqrt(var(res$summary.fitted.values$mean[(37775 + 1):75550]))
sqrt(var(caa_contiguous[which(caa_contiguous$caa_violations_bin == 1), "caa_violations"]))

str(res$summary.fitted.values$mean[1:(75550/2)])

plot(res$residuals, res$summary.fitted.values$mean)

# m <- res$internal.marginals.hyperpar[[1]] 
# m.var <- inla.tmarginal(function(x) 1 / exp(x), m)
# inla.zmarginal(m.var)


library(tmap)
library(RColorBrewer)
caa_contiguous_sf$caa_violations_bin_tx <- ifelse(caa_contiguous_sf$caa_violations_bin == 1,
                                                  "yes", "no")

prop.table(table(caa_contiguous$caa_violations == 0))

states <- st_read("/Users/brenna/Documents/School/2023 Spring, Summer/GEOG 6960/us states/states.shp")
states <- states[which(!states$NAME %in% c("Alaska", "Hawaii", "Puerto Rico")), ]

tm_shape(states, proj = aea) +
  tm_polygons(col = "white", border.col = "gray10") +
  tm_shape(caa_contiguous_sf) +
  tm_dots(col = "caa_violations_bin_tx", palette = "-RdYlBu", alpha = 0.8)
  
count_state <- aggregate(caa_contiguous$caa_violations, by = list(caa_contiguous$fac_fips_code), FUN = sum)
count_state$Group.1 <- str_pad(count_state$Group.1, width = 5, side = "left", pad = "0")

county_shp <- st_read("/Users/brenna/Documents/School/Thesis/stream-pollution-analysis/data/counties.shp")
county_shp <- merge(county_shp, count_state, by.x = "GEOID", by.y = "Group.1", all.y = TRUE)
head(county_shp)
county_shp$Violations <- county_shp$x
county_shp$Count <- log(county_shp$Violations)
county_shp_pois <- county_shp[which(county_shp$Count != "-Inf"), ]
summary(county_shp_pois)
tm_shape(county_shp_pois, proj = aea) +
  tm_polygons(col = "Count", style = "cont", lwd = 0, palette = "BuPu")

ggplot(caa_contiguous, aes(x = caa_violations)) +
  geom_histogram(bins = 100)

sector_bin


rbind(creds, c(h_cred, m_cred, w_cred))

inla.pmarginal(0.89, res$marginals.random$sec_b$index.1)
1 - inla.pmarginal(0.99, res$marginals.random$sec_b$index.2)

inla.pmarginal(0.67, res$marginals.fixed$x_bla_z)
inla.pmarginal(0.89, res$marginals.fixed$x_bla_z)
inla.pmarginal(0.99, res$marginals.fixed$x_bla_z)



data.frame(sector_bin$ID, round(exp(sector_bin$mean), 2), sector_bin$cred)

# credible: 21, 22, 23, 31, 44, 48, 51, 56, 81

plot(res, plot.fixed.effects = FALSE,
     plot.random.effects = FALSE,
     plot.hyperparameters = TRUE,
     plot.predictor = FALSE, cex = 1.25)

table(caa_contiguous$caa_days_last_evaluation,
      caa_contiguous$fac_state)

plot(res, plot.fixed.effects = TRUE,
     plot.random.effects = FALSE,
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
caa_contiguous$fitted <- res$summary.fitted.values$mean[1:37775]
caa_contiguous$fitted_p <- res$summary.fitted.values$mean[37776:75550]

res$summary.random


caa_sf <- st_as_sf(caa_contiguous, coords = c("fac_long", "fac_lat"),
                   crs = 4326, agr = "constant")
caa_sf

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

z_st$state <- str_sub(z_st$GEOID, start = 0, end = 2)

aggregate(z_st$re_st_b, by = list(z_st$state), FUN = mean)


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

hi <- caa_sf[which(caa_sf$fitted > 0.6), ]
c(hi$fac_name, hi$fac_state)
table(hi$industry_sector)

huntsman <- caa_contiguous |>
  filter(grepl("HUNTSMAN", fac_name))

plot(hurdle_dat$re_b)

summary(res$summary.random$id_space_nb)

tx <- caa_contiguous |>
  filter(fac_state == "TX")

table(tx$caa_compliance_status)
tx$y_hat_b_t <- exp(tx$y_hat_b) / (1 + exp(tx$y_hat_b))
aggregate(tx$y_hat_b_t, by = list(tx$caa_compliance_status), FUN = mean)


ca <- caa_contiguous |>
  filter(fac_state == "CA")
ca$y_hat_b_t <- exp(ca$y_hat_b) / (1 + exp(ca$y_hat_b))
aggregate(ca$y_hat_b_t, by = list(ca$caa_compliance_status), FUN = mean)

table(ca$caa_compliance_status)
head(tx)

table(caa_sf$industry_sector_simple)
caa_sf$industry_sector_simple <- case_when(caa_sf$industry_sector == "manufacturing" ~ "manufacturing",
                                           caa_sf$industry_sector == "natural resources and mining" ~ "natural resources and mining",
                                           caa_sf$industry_sector == "trade, transportation and utilities" ~ "trade, transportation and utilities",
                                           caa_sf$industry_sector %in% c("construction",
                                                                         "education and health services",
                                                                         "financial activities",
                                                                         "information",
                                                                         "leisure and hospitality",
                                                                         "nonclassifiable",
                                                                         "other services (except public admin)",
                                                                         "professional and business services",
                                                                         "public administration") ~ "other")
# map of industry
tm_shape(z_st, proj = aea) +
  tm_polygons(col = "white", alpha = 0.5, lwd = 0.5) +
  tm_shape(caa_sf, proj = aea) +
  tm_dots(col = "industry_sector_simple", palette = "Dark2", alpha = 0.6, size = 0.05) + #, 
  #breaks = seq(0, 1, by = 0.1)) +
  tm_layout(#aes.palette = list(seq = "-RdYlGn"), 
    title = "Industry", 
    title.position = c(0.65, 0.95))

table(caa_contiguous$caa_compliance_status, caa_contiguous$fac_state)
aggregate(caa_contiguous$caa_compliance_status,
          by = list(caa_contiguous$fac_state), FUN = mean)

max(caa_contiguous$caa_violations)

caa_contiguous$caa_compliance_status_simple <- case_when(caa_contiguous$caa_compliance_status %in% 
                                                           c("Violation Addressed; EPA Has Lead Enforcement",
                                                             "Violation Unaddressed; EPA Has Lead Enforcement") ~
                                                           "epa lead enforcement",
                                                         caa_contiguous$caa_compliance_status %in% 
                                                           c("Violation Addressed; Local Has Lead Enforcement",
                                                             "Violation Unaddressed; Local Has Lead Enforcement") ~
                                                           "local lead enforcement",
                                                         caa_contiguous$caa_compliance_status %in% 
                                                           c("Violation Addressed; State Has Lead Enforcement",
                                                             "Violation Unaddressed; State Has Lead Enforcement") ~
                                                           "state lead enforcement",
                                                         caa_contiguous$caa_compliance_status %in% 
                                                           c("Violation Identified",
                                                             "Violation w/in 1 Year",
                                                             "Violation-Unresolved") ~
                                                           "no lead enforcement",
                                                         caa_contiguous$caa_compliance_status %in% 
                                                           c("No Violation Identified") ~
                                                           "no violation")

table(caa_contiguous$caa_compliance_status_simple)

caa_contiguous$caa_compliance_status_simple <- relevel(as.factor(caa_contiguous$caa_compliance_status_simple),
                                                       ref = "no violation")

summary(inla(caa_violations ~ caa_compliance_status_simple, 
             data = caa_contiguous, family = "poisson"))

aggregate(caa_contiguous$caa_violations, by = list(caa_contiguous$caa_compliance_status_simple), FUN = mean)
table(caa_contiguous$caa_compliance_status_simple, caa_contiguous$fac_state,
      caa_contiguous$caa_violations_bin)

head(caa_contiguous)

tail(caa_contiguous$fec_case_ids)

dlist <- spdep::nbdists(nb.gab, xy)
listw.d1 <- spdep::nb2listw(nb.gab, style = "W", zero.policy = TRUE, glist=dlist)
dlist <- lapply(dlist, function(x) 1/x)
dlist <- lapply(dlist, function(x) 1/x^2)
listw.d2 <- spdep::nb2listw(nb.gab, style = "W", zero.policy = TRUE, glist=dlist)

table(caa_contiguous$fec_case_ids != "", caa_contiguous$fac_state)

head(unique(caa_contiguous$caa_3yr_compl_qtrs_history))

table(caa_contiguous$caa_hpv_flag, caa_contiguous$industry_sector)

aggregate(caa_contiguous$caa_days_last_evaluation, by = list(caa_contiguous$caa_violations_bin), FUN = mean)


caa_contiguous$fac_derived_cb2010

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

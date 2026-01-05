library(tidyverse)
library(lubridate)
library(evd)
library(eva)
library(readxl)
library(ismev)

df <- read_excel("Base X Y.xlsx")

colnames(df) <- c("date", "X", "Y")

df <- df %>%
  mutate(
    date = as.Date(date),
    X = as.numeric(X),
    Y = as.numeric(Y)
  )

#Plot data
df_long <- df %>%
  pivot_longer(cols = c(X, Y), names_to = "Portefeuille", values_to = "Perte")

ggplot(df_long, aes(x = date, y = Perte, color = Portefeuille)) +
  geom_line() +
  facet_wrap(~Portefeuille, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(title = "Historique des pertes (X et Y)", x = "Date", y = "Montant")

print(head(df)) #print premières lignes
summary(df)  #print statistiques

#Loi normale 

mu_X <- mean(df$X)
sigma_X <- sd(df$X)

print(paste("Moyenne :", round(mu_X, 2)))
print(paste("Ecart-type :", round(sigma_X, 2)))

hist(df$X, freq = FALSE, breaks = 100, col = "lightblue", 
     main = "Distribution de X vs Loi Normale",
     xlab = "Montant des pertes X")

curve(dnorm(x, mean = mu_X, sd = sigma_X), 
      add = TRUE, col = "red", lwd = 3)


# Modèle GEV 1 (annuel)

df_gev <- df %>%
  mutate(annee = year(date)) %>%
  group_by(annee) %>%
  summarise(X_max = max(X))

print("Aperçu des maxima annuels :")
print(head(df_gev))

modele_gev <- fgev(df_gev$X_max, method = "Nelder-Mead") #fgev trouve les 3 paramètres (fitting GEV)

print(modele_gev)  #On obtient mu=32 040, sigma= 9631 et Xi=0.45, donc Domaine de Fréchet)

par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
plot(modele_gev, ask = FALSE) 
mtext("Analyse du Modèle GEV (Annuel)", 
      outer = TRUE,
      cex = 1.5,
      font = 2,
      line = 1)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))

# Modèle GEV 2 (semestriel)
df_gev_semestriel <- df %>%
  mutate(
    annee = year(date),
    semestre = semester(date) 
  ) %>%
  group_by(annee, semestre) %>%
  summarise(X_max = max(X), .groups = 'drop')

print("Aperçu des maxima semestriels (2 points par an) :")
print(head(df_gev_semestriel))

modele_gev_sem <- fgev(df_gev_semestriel$X_max, method = "Nelder-Mead")

print("Paramètres GEV (Semestriel)")
print(modele_gev_sem)

par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
plot(modele_gev_sem, ask = FALSE) 
mtext("Analyse du Modèle GEV (Semestriel)", 
      outer = TRUE,
      cex = 1.5,
      font = 2,
      line = 1)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))

#GPD

#Mean Excess Plot

par(mfrow = c(1, 1))
mrl.plot(df$X)

u_choisi <- 35000

print(paste("Seuil choisi pour le modèle GPD :", u_choisi))

modele_gpd <- gpd.fit(df$X, threshold = u_choisi, show = FALSE)

# Affichage des paramètres (Scale = Sigma, Shape = Xi)
print("Paramètres estimés du modèle GPD :")
print(modele_gpd$mle) 

gpd.diag(modele_gpd)

#Choix du modèle, on compare tous les graphs :

#Loi Normale : inadaptée: ne prend pas en compte la queue (fat tails), la fat tails est prouvée par le Xi positifi des GPD et GEV

#GEV annuel : incertitude énorme, densité qui ne prend pas assez en compte les fat tails, probability plot ok
#GEV semestrielle:Meilleure incertitude sur le QQ plot mais différence de 100k-180k sur le dernier quantile, meilleure densité mais probability plot incohérent

#On choisit un seuil à 35k pour la GPD (sur le rendu en mettre plusieurs pour comparer)
#GPD : Probability plot ok, return level plot ok, density plot ok (première barre d'histogramme qui ,ne match pas car seuil), QQ plot légèrement meilleure mais toujours gap avec le dernier point

#On choisit la GPD

#Quantiles :

p_cible <- 0.995 #quantile 99.5%
N_mois_an <- 12  #pour GEV annuel
N_mois_sem <- 6  #pour GEV semestrielle

#Loi normale
val_normale <- qnorm(p_cible, mean = mu_X, sd = sigma_X)
n_total <- length(df$X)
densite_norm <- dnorm(val_normale, mean = mu_X, sd = sigma_X)
se_normale <- (sigma_X / densite_norm) * sqrt((p_cible * (1 - p_cible)) / n_total)

ic_inf_norm <- val_normale - 1.96 * se_normale
ic_sup_norm <- val_normale + 1.96 * se_normale


# GEV annuel
prob_gev_annuel <- p_cible^N_mois_an  #car les données sont des max annuels

val_gev_annuel <- qgev(prob_gev_annuel, 
                       loc = modele_gev$estimate[1], 
                       scale = modele_gev$estimate[2], 
                       shape = modele_gev$estimate[3])

#calcul "à la main", formule du cours :

mu_an   <- modele_gev$estimate[1] # mu
sigma_an <- modele_gev$estimate[2] # sigma
xi_an    <- modele_gev$estimate[3] # xi
cov_mat  <- modele_gev$var.cov

p_adj <- p_cible^N_mois_an  

term_A <- -log(p_adj) 

val_gev_annuel_main <- mu_an - (sigma_an / xi_an) * (1 - term_A^(-xi_an)) #car Xi>0

dg_dmu <- 1 #dérivée partielle mu

dg_dsigma <- -(1 / xi_an) * (1 - term_A^(-xi_an)) #dérivée partielle sigma

term_xi <- term_A^(-xi_an) # (-ln(1-alpha))^(-xi)

dg_dxi <- (sigma_an / xi_an^2) * (1 - term_xi) - (sigma_an / xi_an) * log(term_A) * term_xi #dérivée partielle Xi

Dg <- c(dg_dmu, dg_dsigma, dg_dxi) #gradient

# Variance de l'estimateur : Var = Dg' * Sigma * Dg
variance_var <- t(Dg) %*% cov_mat %*% Dg
se_var <- sqrt(as.numeric(variance_var)) # Ecart-type

# Intervalle de confiance 95%
ic_inf_gev1 <- val_gev_annuel_main - 1.96 * se_var  #1.96 quantile à 95% loi normale
ic_sup_gev1 <- val_gev_annuel_main + 1.96 * se_var

#GEV semestriel
prob_gev_sem <- p_cible^N_mois_sem  #car les données sont des max semestriels

val_gev_sem <- qgev(prob_gev_sem, 
                    loc = modele_gev_sem$estimate[1], 
                    scale = modele_gev_sem$estimate[2], 
                    shape = modele_gev_sem$estimate[3])

mu_an2    <- modele_gev_sem$estimate[1]
sigma_an2 <- modele_gev_sem$estimate[2]
xi_an2    <- modele_gev_sem$estimate[3]
cov_mat2   <- modele_gev_sem$var.cov

p_adj2 <- p_cible^N_mois_sem 

term_A2 <- -log(p_adj2)
val_gev_sem_main <- mu_an2 - (sigma_an2 / xi_an2) * (1 - term_A2^(-xi_an2))

term_xi2 <- term_A2^(-xi_an2)

dg_dmu2 <- 1 #d mu
dg_dsigma2 <- -(1 / xi_an2) * (1 - term_xi2) #d sigma
dg_dxi2    <- (sigma_an2 / xi_an2^2) * (1 - term_xi2) - (sigma_an2 / xi_an2) * log(term_A2) * term_xi2 #d xi

Dg2 <- c(dg_dmu2, dg_dsigma2, dg_dxi2)
variance_var2 <- t(Dg2) %*% cov_mat2 %*% Dg2
se_var2 <- sqrt(as.numeric(variance_var2))

ic_inf_gev2 <- val_gev_sem_main - 1.96 * se_var2
ic_sup_gev2 <- val_gev_sem_main + 1.96 * se_var2

#GPD
scale_gpd <- modele_gpd$mle[1]
shape_gpd <- modele_gpd$mle[2]
zeta      <- modele_gpd$rate 
cov_mat3   <- modele_gpd$cov

val_gpd <- u_choisi + (scale_gpd/shape_gpd) * ( ((1 - p_cible)/zeta)^(-shape_gpd) - 1 )

prob_ratio <- (1 - p_cible) / zeta

term_A <- prob_ratio^(-shape_gpd)

dg_dsigma3 <- (1 / shape_gpd) * (term_A - 1) 

dg_dxi3 <- -(scale_gpd / shape_gpd^2) * (term_A - 1) - (scale_gpd / shape_gpd) * term_A * log(prob_ratio)

Dg3 <- c(dg_dsigma3, dg_dxi3)

variance_gpd <- t(Dg3) %*% cov_mat3 %*% Dg3
se_gpd <- sqrt(as.numeric(variance_gpd))

ic_inf_gpd <- val_gpd - 1.96 * se_gpd
ic_sup_gpd <- val_gpd + 1.96 * se_gpd

df_recap_final <- data.frame(
  Modele = c("Loi Normale", "GEV (Annuel)", "GEV (Semestriel)", "GPD (POT - Retenu)"),

  VaR_99.5 = round(c(val_normale, val_gev_annuel_main, val_gev_sem_main, val_gpd), 2),

  IC_Bas_95 = round(c(ic_inf_norm, ic_inf_gev1, ic_inf_gev2, ic_inf_gpd), 2),

  IC_Haut_95 = round(c(ic_sup_norm, ic_sup_gev1, ic_sup_gev2, ic_sup_gpd), 2),
  
  Largeur_IC = round(c(ic_sup_norm - ic_inf_norm, 
                       ic_sup_gev1 - ic_inf_gev1, 
                       ic_sup_gev2 - ic_inf_gev2, 
                       ic_sup_gpd - ic_inf_gpd), 2)
)

print(df_recap_final)

#Perte annuelle periode de retour 200 ans

T_annees <- 200
N_obs_an <- 12
N_periode <- T_annees * N_obs_an

prob_tail_200 <- 1 / N_periode 

prob_ratio_200 <- prob_tail_200 / zeta

val_200_ans <- u_choisi + (scale_gpd / shape_gpd) * ( prob_ratio_200^(-shape_gpd) - 1 )

#intervalle de confiance, delta methode
term_A_200 <- prob_ratio_200^(-shape_gpd)

dg_dsigma_200 <- (1 / shape_gpd) * (term_A_200 - 1) 

dg_dxi_200 <- -(scale_gpd / shape_gpd^2) * (term_A_200 - 1) - (scale_gpd / shape_gpd) * term_A_200 * log(prob_ratio_200)

Dg_200 <- c(dg_dsigma_200, dg_dxi_200)

variance_200 <- t(Dg_200) %*% cov_mat3 %*% Dg_200
se_200 <- sqrt(as.numeric(variance_200))

# Interval de confiance à 95%
ic_inf_200 <- val_200_ans - 1.96 * se_200
ic_sup_200 <- val_200_ans + 1.96 * se_200

cat(" Niveau de retour 200 ans\n")
cat("Modèle utilisé : GPD (Seuil =", u_choisi, ")\n")
cat("Probabilité mensuelle équivalente : 1 /", N_periode, "\n\n")

df_retour_200 <- data.frame(
  Estimation_Central = round(val_200_ans, 2),
  IC_Bas_95 = round(ic_inf_200, 2),
  IC_Haut_95 = round(ic_sup_200, 2),
  Largeur_IC = round(ic_sup_200 - ic_inf_200, 2)
)


print(df_retour_200)

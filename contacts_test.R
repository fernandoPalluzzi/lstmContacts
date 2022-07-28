# R

source("/home/fernando/contactsCore/contacts.R")

#names(contacts.map)
#lapply(contacts.map, function(x) names(x))

ab.library <- unlist(lapply(contact.map, function(x) names(x)))
ab.info <- aa.info(ab.library)
ab.groups <- aa.group(ab.info$aa)

# Mock
x <- c("h.D105", "Q449", "W455", "M492", "L493", "K494")

# Real
#x <- c("h.R105", "Y449", "L455", "L492", "Q493", "S494")
#x <- c("h.Y33", "K417", "Y421")
#x <- c("h.Y33", "K417", "Y421", "L455", "F456")
#x <- c("Y33", "V483", "T478", "F486")
#x <- c("h.R50", "V483", "E484")
#x <- c("l.Y92", "Y380", "G381", "P412", "D427", "F429", "D482")

x.ag <- x[2:length(x)]

X <- ab.search(x[1], ab.library)
X

#ab <- ab.select(X)
ab <- ab.select(X, keepAll = TRUE)

res <- res.fetch(ab)

R <- res.search(res, x.ag)
R <- residue.class(R, verbose = TRUE)
R

Ri <- res.extract(R, X, x)
Ri

K <- unique(names(Ri$ab.residues))
J <- which(names(contact.map) %in% paste0("ab.", Ri$antibody))
current.map <- contact.map[J]

variants <- vector()
for (i in 1:length(current.map)) {
	ki <- Ri$ab.residues[names(Ri$ab.residues) == K[i]]
	ci <- current.map[[i]]
	ci <- ci[names(ci) %in% ki]
	hot <- contact.mutations[contact.mutations$wt %in% unlist(ci),]
	mut0 <- contact.mutations[contact.mutations$mutant %in% x.ag,]
	mut1 <- contact.mutations[contact.mutations$mutant %in% Ri$class.residues,]
	if (nrow(mut0) > 0) {
		v <- unique(mut0$variant)
	} else if (nrow(mut1) > 0) {
		v <- unique(mut1$variant)
	} else if (nrow(hot) > 0) {
		v <- unique(hot$variant)
	} else {
		v <- NULL
	}
	variants <- c(variants, v)
}

Ri$variant <- variants





# R

source("/home/fernando/contactsCore/contacts.R")
source("/home/fernando/euNet/euNet_core/rutils.R")


#data <- prepareLibrary(a0 = 0.88)
#write.delim(data, "/home/fernando/contactsCore/contactLibrary_t10.txt")

#data <- prepareLibrary(chunk = 5, a0 = 0.88)
#write.delim(data, "/home/fernando/contactsCore/contactLibrary_t5.txt")

P <- extractProfiles()
group <- c(0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#group <- ifelse(apply(P, 2, median) > 0.88, 0, 1)
#group <- ifelse(apply(P, 2, median) > 0.8, 0, 1)
#group <- rep(c(0, 1, 1, 1, 1), 7)
#save(A, file = "~/GemelliBioinfoUnit/Stratibody/Stratibody_affinity_groups.RData")

A0 <- apply(P[, group == 0], 1, mean)
A1 <- apply(P[, group == 1], 1, mean)


pdf("~/GemelliBioinfoUnit/Stratibody/Stratibody_affinity_threshold_auc.pdf", width = 20, height = 10)
plot(density(A0, bw = 0.025), col = "blue", lwd = 3, lty = 1,
     main = "",
     xlim = c(0.4, 1),
     xlab = "Affinity score",
     cex.axis = 1.8,
     cex.lab = 1.4)
lines(density(A1, bw = 0.06), col = "red2", lwd = 3)
abline(v = 0.8815, lty = 2, lwd = 3, col = "black")
text(x = 0.8, y = 10,
     "Affinity cutoff = 0.88\nAUC = 0.941\nCI95% = (0.904, 0.979)",
     cex = 1.8)
legend("topleft", fill = c("blue", "red2", "black"), bg = "white",
legend = c("Stable clusters distribution",
           "Unstable clusters distribution",
           "Affinity score thrshold (0.88)"),
lty = c(1, 1, 2),
cex = 1.6)
dev.off()


library(OptimalCutpoints)

load("~/GemelliBioinfoUnit/Stratibody/Stratibody_affinity_groups.RData")

cutpoints <- function(data, xmin = 0.75, xmax = 0.95, xstep = 0.01,
                      healthy = 0) {
	cutoffs <- vector()
	auc <- vector()
	lower <- vector()
	upper <- vector()
	for (x in seq(xmin, xmax, xstep)) {
		group <- ifelse(apply(data, 2, median) > x, 0, 1)
		A0 <- apply(data[, group == 0], 1, mean)
		A1 <- apply(data[, group == 1], 1, mean)
		A <- data.frame(affinity = c(A1, A0), y = c(rep(1, length(A1)),
		                                            rep(0, length(A0))))
		w <- optimal.cutpoints(X = "affinity", status = "y",
		                       tag.healthy = healthy,
		                       methods = "SpEqualSe",
		                       data = A)
		cutoffs <- c(cutoffs, w$SpEqualSe$Global$optimal.cutoff$cutoff[1])
		auc <- c(auc, w$SpEqualSe$Global$measures.acc$AUC[1])
		lower <- c(lower, w$SpEqualSe$Global$measures.acc$AUC[2])
		up <- w$SpEqualSe$Global$measures.acc$AUC[3]
		if (up > 1) up <- 1
		upper <- c(upper, up)
	}
	return(data.frame(cutoffs = cutoffs, auc = auc,
	                  lower = lower, upper = upper))
}


cutoffs <- cutpoints(P, healthy = 1)

#cutoffs$cutoffs <- cutoffs$cutoffs - 0.0045697


opt <- optimal.cutpoints(X = "affinity", status = "y",
                         tag.healthy = 1,
                         methods = "SpEqualSe",
                         data = A)
opt

# Optimal cutoff: 0.8824
# AUC: 0.941 (0.904, 0.979)





# R

source("/home/fernando/contactsCore/contacts.R")
source("/home/fernando/euNet/euNet_core/rutils.R")

P <- extractProfiles()

if (FALSE) {
#           Best: 7l7d-WT
#          Worst: 7kmg-BETA
# False positive: 7l7d-DELTA

# Mock
x <- c("h.D105", "Q449", "W455", "M492", "L493", "K494")

# Real
#x <- c("h.R105", "Y449", "L455", "L492", "Q493", "S494")
#x <- c("h.Y33", "K417", "Y421")
#x <- c("h.Y33", "K417", "Y421", "L455", "F456")
#x <- c("h.Y33", "V483", "T478", "F486")
#x <- c("h.R50", "V483", "E484")
#x <- c("l.Y92", "Y380", "G381", "P412", "D427", "F429", "D482")

R <- preprocess(x)
R
}


# Complex

# 7kmg
x <- list(x1 = c("h.R50", "V483", "E484"),
          x2 = c("h.L55", "L452", "T470", "F490"),
          x3 = c("h.Y101", "E484", "F490"),
          x4 = c("h.R104", "Q493", "S494"),
          x5 = c("l.Y32", "F486", "Y489"),
          x6 = c("l.Y92", "F486", "Y489"),
          x7 = c("l.R96", "V483", "E484"))

R <- contacts(x)

profile <- extractProfiles(data = contact.data,
                           antibody = R$antibody,
                           variants = R$variant,
                           residues = R$ab.residues,
                           stochastic = TRUE)

worst <- profile$ab.7kmg.beta


# 7l7d
x <- list(x1 = c("h.W50", "F486", "N487"),
          x2 = c("h.N107", "F486", "N487", "Y489"),
          x3 = c("h.F110", "F486"),
          x4 = c("h.Y33", "V483", "T478", "F486"),
          x5 = c("h.Y92", "F486"),
          x6 = c("h.W98", "F486"))

R <- contacts(x)

profile <- extractProfiles(data = contact.data,
                           antibody = R$antibody,
                           variants = R$variant,
                           residues = R$ab.residues,
                           stochastic = TRUE)

best <- profile$ab.7l7d.wt
fp <- profile$ab.7l7d.delta

#if (class(profile) == "data.frame") profile <- apply(profile, 1, "median")
#plot(profile, type = "l", lwd = 2, col = "blue", ylim = c(0, 1))

plot(best, type = "l", lwd = 4, col = "green4", ylim = c(0, 1))
lines(worst, type = "l", lwd = 4, col = "red2")
lines(fp, type = "l", lwd = 4, col = "darkorange")

median(best)
median(worst)
median(fp)


# LSTM-RNN results

P <- read.csv("~/GemelliBioinfoUnit/Stratibody/RISULTATI_STRATIBODY/STRATIBODY.csv")
names(P)[4] <- "EC"

worst <- c(100, 91, 84, 81, 88, 84, 83, 75, 79, 73, 74, 73, 71, 72, 71,
           72, 69, 60, 67, 75, 73, 60, 66, 62, 65, 63, 64, 68, 64, 64,
           66, 66, 64, 69, 65, 65, 65, 68, 61, 65, 66, 65, 51, 55, 53,
           55, 54, 53, 54, 51, 50, 52, 52, 51, 49, 52, 48, 53, 53, 51,
           52, 55, 55, 55, 54, 55, 52, 59, 56, 46, 55, 51, 54, 55, 55,
           48, 56, 48, 53, 47, 52, 56, 53, 53, 56, 55, 59, 53, 48, 56,
           49, 43, 45, 52, 56, 55, 52, 47, 44, 62)

best <- c(100, 100, 95, 98, 96, 100, 100, 95, 98, 96, 100, 100, 95, 98,
          96, 100, 100, 95, 98, 96, 100, 99.5, 95, 98, 96, 100, 100, 95,
          98, 96, 100, 100, 95, 98, 96, 100, 100, 95, 98, 96, 100, 100,
          95, 98, 96, 100, 100, 95, 98, 96, 100, 100, 95, 98, 96, 100,
          100, 95, 98, 96, 100, 100, 95, 98, 96, 100, 100, 95, 98, 96,
          100, 100, 98, 98, 96, 100, 100, 95, 98, 96, 100, 97, 95, 97.5,
          96, 99.5, 97, 96.5, 98, 96, 100, 98.5, 95, 98, 96, 98.5, 98.5,
          95, 98, 96)

fp <- c(100, 98.5, 98, 98, 96, 100, 97, 95, 98, 96, 100, 98, 98, 98, 96,
        94.5, 97, 95, 98, 96, 98, 95.5, 95, 98, 93.5, 91.5, 91, 93, 98,
        92, 91, 97, 88, 95, 91, 95, 94, 91.5, 94, 89.5, 97, 98, 95, 97,
        96, 96, 95, 97, 95, 95.5, 94, 96, 90.5, 93, 85, 94, 91, 92, 90,
        94, 93, 92, 94, 93, 92, 94, 93, 93, 94, 91, 87, 88, 88, 92, 89,
        92, 92, 88, 86, 93, 90, 92, 93, 86, 85, 95, 88, 84, 88, 84, 88,
        93, 86, 91, 89, 86, 89, 91, 89, 82)

pdf("~/GemelliBioinfoUnit/Stratibody/Stratibody_molDynamics_vs_LSTM_prediction.pdf", width = 20, height = 10)
plot(P$score[P$EC == "BETA-LyCov555"], type = "l",
     lwd = 4, col = "darkred", lty = 2,
     ylim = c(0.2, 1),
     xlab = "nanoseconds",
     ylab = "Affinity score",
     cex.axis = 1.8,
     cex.lab = 1.4)
lines(P$score[P$EC == "WT-AZD8895"], type = "l", lwd = 4, col = "darkgreen", lty = 2)
lines(P$score[P$EC == "DELTA-AZD8895"], type = "l", lwd = 4, col = "darkorange3", lty = 2)
lines(best/100, type = "l", lwd = 4, col = "green3")
lines(worst/100, type = "l", lwd = 4, col = "red2")
lines(fp/100, type = "l", lwd = 4, col = "darkorange")
abline(h = 0.88, lwd = 5, lty = 3)
text(x = -0.2, y = 0.86, "0.88", cex = 1.8)
legend("bottomleft", fill = c("green3", "darkgreen", "darkorange",
                              "darkorange3", "red3", "darkred", "black"),
                     bg = "white",
legend = c("WT-AZD8895 LSTM prediction",
           "WT-AZD8895 molecular dynamics",
           "DELTA-AZD8895 LSTM prediction",
           "DELTA-AZD8895 molecular dynamics",
           "BETA-LyCov555 LSTM prediction",
           "BETA-LyCov555 molecular dynamics",
           "Affinity score threshold (0.88)"),
lty = c(1, 2, 1, 2, 1, 2, 3),
cex = 1.6)
dev.off()





# MERGE FUNCTION





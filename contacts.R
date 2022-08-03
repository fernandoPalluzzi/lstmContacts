# R 

load("~/lstmContacts/RData/Contacts_AminoAcids.RData")
load("~/lstmContacts/RData/Contacts_AminoAcid_groups.RData")
load("~/lstmContacts/RData/Contacts_map.RData")
load("~/lstmContacts/RData/Contacts.RData")
load("~/lstmContacts/RData/Contacts_abClass.RData")
load("~/lstmContacts/RData/Contacts_mutations.RData")
load("~/lstmContacts/RData/affinity_groups.RData")
load("~/lstmContacts/RData/affinity_cutoffs.RData")
load("~/lstmContacts/RData/molDynamics/md_7kmg_beta.RData")
load("~/lstmContacts/RData/molDynamics/md_7l7d_wt.RData")
load("~/lstmContacts/RData/molDynamics/md_7l7d_beta.RData")
load("~/lstmContacts/RData/molDynamics/md_7l7d_delta.RData")

aa.info <- function(x, ab = TRUE) {
	aa <- vector()
	aj <- vector()
	if (ab) {
		ch = vector()
		for (i in 1:length(x)) {
			a <- strsplit(x[i], "")[[1]]
			ch <- c(ch, a[1])
			aa <- c(aa, a[3])
			aj <- c(aj, paste0(a[3:length(a)], collapse = ""))
		}
		return(list(chain = ch, aa = aa, residue = aj))
	} else {
		for (i in 1:length(x)) {
			a <- strsplit(x[i], "")[[1]]
			aa <- c(aa, a[1])
			aj <- c(aj, x[i])
		}
		return(list(chain = NULL, aa = aa, residue = aj))
	}
}

aa.group <- function(x) {
	v <- vector()
	for (i in 1:length(x)) {
		v <- c(v, contact.groups$group[contact.groups$code == x[i]])
	}
	return(v)
}

ab.search <- function(res, lib = NULL) {
	
	ab <- unlist(lapply(names(lib),
				 function(x) paste0(strsplit(x, "")[[1]][4:7],
				 collapse = "")))
	names(lib) <- ab
	
	aa <- lib[lib == res]
	if (length(aa) == 0) aa <- NULL
	
	rma <- lib[!(names(lib) %in% names(aa))]
	chain.res <- unlist(lapply(rma,
						function(x) paste0(strsplit(x, "")[[1]][1:3],
						collapse = "")))
	chain.ab <- paste0(strsplit(res, "")[[1]][1:3], collapse = "")
	
	lv3 <- rma[chain.res == chain.ab]
	if (length(lv3) == 0) lv3 <- NULL
	
	rm3 <- rma[!(names(rma) %in% names(lv3))]
	rm3.group <- unlist(lapply(res.neat(rm3),
				   function(x) contact.groups$group[contact.groups$code == x]))
	
	x.neat <- res.neat(res)
	x.group <- unlist(lapply(x.neat,
				 function(x) contact.groups$group[contact.groups$code == x]))
	
	x.chain <- strsplit(res, "")[[1]][1]
	rescued <- (rm3.group %in% x.group) & startsWith(rm3, x.chain)
	lv2 <- rm3[rescued]
	if (length(lv2) == 0) lv2 <- NULL
	
	rm2 <- rm3[!(names(rm3) %in% names(lv2))]
	rm2.neat <- res.neat(rm2)
	
	lv1 <- rm2[rm2.neat %in% x.neat]
	if (length(lv1) == 0) lv1 <- NULL
	
	w <- 0
	if (is.null(aa)) w <- w + 8
	if (is.null(lv3)) w <- w + 4
	if (is.null(lv2)) w <- w + 2
	if (is.null(lv1)) w <- w + 1
	if (w == 8) {
		warning(paste0("No exact antibody match found for: ",
	                   res), collapse = "")
	} else if (w > 8 & w < 15) {
		warning("Partial antibody match.")
	} else if (w == 15) {
		stop("No suitable antibody was found in the current library: exiting.")
	}
	
	return(list(exact.match = aa, chain.match = lv3, group.match = lv2,
	            aa.match = lv1, warning.level = w))
}

ab.select <- function(ab, exact = 2, level = 1, keepAll = FALSE) {
	if (!keepAll & length(ab$exact.match) >= exact) {
		v <- ab$exact.match
	} else if (level == 1 | keepAll) {
		v <- c(ab$exact.match, ab$chain.match, ab$aa.match, ab$group.match)
	} else if (level == 2) {
		v <- c(ab$exact.match, ab$chain.match, ab$aa.match)
	} else if (level == 3) {
		v <- c(ab$exact.match, ab$chain.match)
	}
	return(v)
}

res.fetch <- function(v) {
	u <- names(v)
	res <- vector()
	var.library <- contact.map[names(contact.map) %in% paste0("ab.", u)]
	for (i in 1:length(v)) {
		command <- paste0("x <- var.library$ab.", u[i], "$", v[i],
		                  collapse = "")
		eval(parse(text = command))
		a <- rep(u[i], length(x))
		names(x) <- a
		res <- c(res, x)
	}
	return(res)
}

res.neat <- function(v) {
	R <- vector()
	for (x in v) {
		x <- strsplit(x, "\\.")[[1]]
		if (length(x) > 1) {
			x <- x[2]
		} else {
			x <- x[1]
		}
		R <- c(R, strsplit(x, "")[[1]][1])
	}
	names(R) <- names(v)
	return(R)
}

res.search <- function(R, ag) {
	
	lv3 <- R[R %in% ag]
	rm3 <- R[!(names(R) %in% names(lv3))]
	if (length(lv3) == 0) lv3 <- NULL
	
	rm3.neat <- res.neat(rm3)
	x.neat <- res.neat(ag)
	lv2 <- rm3[rm3.neat %in% x.neat]
	if (length(lv2) == 0) lv2 <- NULL
	
	rm2 <- rm3[!(names(rm3) %in% names(lv2))]
	rm2.group <- unlist(lapply(res.neat(rm2),
                   function(x) contact.groups$group[contact.groups$code == x]))
    
    x.group <- unlist(lapply(x.neat,
                 function(x) contact.groups$group[contact.groups$code == x]))
    
    rescued <- rm2.group %in% x.group
    lv1 <- rm2[rescued]
    if (length(lv1) == 0) lv1 <- NULL
    
    lv1.group <- rm2.group[rescued]
    #h <- table(c(ag[!(ag %in% R)], ag[!(x.group %in% rm2.group)]))
    if (is.null(lv3)) {
		f <- table(c(lv1, lv2))
	} else {
		f <- table(lv3)
	}
    h <- table(c(ag[!(ag %in% R)], names(f[f < max(f)])))
    if (dim(table(h))[1] == 0) h <- NULL
    names(x.group) <- ag
    
    w <- 0
	if (is.null(lv3)) w <- w + 8
	if (is.null(lv2)) w <- w + 4
	if (is.null(lv1)) w <- w + 2
	if (length(lv3) < length(ag)) w <- w + 1
	if (w == 8) {
		warning("No perfect residue match found.", collapse = "")
	} else if (w > 8 & w < 15) {
		warning("Modeling possible suboptimal contacts.")
	} else if (w == 15) {
		stop("No suitable contact was found in the current library: exiting.")
	}
	
    return(list(group = x.group, exact.match = lv3, aa.match = lv2,
                group.match = lv1, hotspots = h, warning.level = w))
}

residue.class <- function(R, lib = contact.class, freq = "OR", verbose = FALSE) {
	C <- unlist(lib)
	names(C) <- unlist(lapply(names(C),
	                   function(x) strsplit(x, "\\.")[[1]][2]))
	
	#rclass <- list(class = names(C)[C %in% R$exact.match], level = 3)
	rclass <- list(class = C[C %in% R$exact.match], level = 3)
	
	for (i in 3:4) {
		if (length(rclass$class) < 1) {
			#rclass$class <- names(C)[C %in% R[[i]]]
			rclass$class <- C[C %in% R[[i]]]
			rclass$level <- rclass$level - 1
		}
	}
	if (length(rclass$class) < 1) {
		rclass$class <- NULL
		rclass$level <- 0
	}
	#classes <- unique(rclass$class)
	classes <- rclass$class
	
	n <- c(length(lib$I), length(lib$II), length(lib$III), length(lib$IV))
	
	k <- c(length(rclass$class[names(rclass$class) == "class1"]),
		   length(rclass$class[names(rclass$class) == "class2"]),
		   length(rclass$class[names(rclass$class) == "class3"]),
		   length(rclass$class[names(rclass$class) == "class4"]))
	
	if (verbose) cat("### p\n")
	w <- vector()
	for (i in 1:length(n)) {
		p <- k[i]
		q <- sum(k[-i])
		if (freq == "relative") {
			p <- p/n[i]
			q <- q/sum(n[-i])
		} else if (freq %in% c("or", "OR")) {
			p <- p/(n[i] - p)
			q <- q/(sum(n[-i]) - q)
			p <- (p + 0.001)/(q + 0.001)
			q <- 0
		}
		#if (verbose) print(c(p, q))
		if (verbose) cat(paste0(p, "\n"))
		w <- c(w, ifelse(p > q, p, 0))
	}
	w.max <- which(w == max(w))
	if (verbose) cat(paste0("### Majority class: ", w.max, "\n"))
	K <- paste0("class", w.max)
	
	return(list(group = R$group, exact.match = R$exact.match,
	            aa.match = R$aa.match, group.match = R$group.match,
	            hotspot = R$hotspots, warning.level = R$warning.level,
	            contact.level = rclass$level, contact.classes = classes,
	            majority.class = K))
}

res.extract <- function(R, A, x = NULL) {
	if (length(R$majority.class) ==  2) {
		f <- table(names(R$contact.classes))
		f.max <- names(f)[f == max(f)][1]
		f.min <- names(f)[names(f) != f.max]
		a <- R$contact.classes[names(R$contact.classes) == f.max]
		b <- R$contact.classes[names(R$contact.classes) == f.min]
		if (sum(b %in% a) == length(b)) {
			R$contact.classes <- a
			R$majority.class <- f.max
		} else {
			warning("residue class multiplicity.")
		}
	} else if (length(R$majority.class) > 2) {
		warning("residue class multiplicity.")
	}
	lv <- c(4, 3, 2)
	f <- table(names(R[[lv[R$contact.level]]]))
	ab <- names(f)[f == length(x[2:length(x)])]
	if (length(ab) != 1) ab <- unique(names(R[[lv[R$contact.level]]]))
	a <- c(A$exact.match, A$chain.match, A$group.match, A$aa.match)
	a <- a[names(a) %in% ab]
	res <- R$contact.classes[names(R$contact.classes) %in% R$majority.class]
	names(res) <- NULL
	if (length(R$aa.match) > 0) {
		sim <- R$aa.match
	} else if (length(R$group.match) > 0) {
		sim <- R$group.match
	} else {
		sim <- "None"
	}
	return(list(input = x, antibody = ab, ab.residues = a,
				exact.match = R$exact.match, similarity.match = sim,
	            class.residues = res, class = R$majority.class,
	            warning.level = R$warning.level,
	            contact.level = R$contact.level))
}

preprocess <- function(x, map = contact.map, keepAll = TRUE,
                          show.probs = TRUE, verbose = TRUE) {
	
	ab.library <- unlist(lapply(map, function(x) names(x)))
	
	x.ag <- x[2:length(x)]
	X <- ab.search(x[1], ab.library)
	
	ab <- ab.select(X, keepAll = keepAll)
	res <- res.fetch(ab)
	R <- res.search(res, x.ag)
	R <- residue.class(R, verbose = show.probs)
	R <- res.extract(R, X, x)
	
	K <- unique(names(R$ab.residues))
	J <- which(names(map) %in% paste0("ab.", R$antibody))
	current.map <- map[J]
	
	variants <- vector()
	for (i in 1:length(current.map)) {
		ki <- R$ab.residues[names(R$ab.residues) == K[i]]
		ci <- current.map[[i]]
		ci <- ci[names(ci) %in% ki]
		hot <- contact.mutations[contact.mutations$wt %in% unlist(ci),]
		mut0 <- contact.mutations[contact.mutations$mutant %in% x.ag,]
		mut1 <- contact.mutations[contact.mutations$mutant %in% R$class.residues,]
		if (nrow(mut0) > 0) {
			v <- unique(mut0$variant)
		} else if (nrow(mut1) > 0) {
			v <- unique(mut1$variant)
		} else if (nrow(hot) > 0) {
			v <- unique(hot$variant)
		} else {
			v <- "wt"
		}
		variants <- c(variants, v)
	}
	R$variant <- unique(variants)
	ires <- unlist(lapply(x.ag, function(x) strsplit(x, "")[[1]][1]))
	cres <- unlist(lapply(R$class.residues,
	               function(x) strsplit(x, "")[[1]][1]))
	R$hotspot <- x.ag[!(ires %in% cres)]
	if (length(R$hotspot) == 0) R$hotspot <- "None"
	R$exact.match <- R$exact.match[names(R$exact.match) %in% R$antibody]
	if (length(R$exact.match) == 0) R$exact.match <- "None"
	R$similarity.match <- R$similarity.match[names(R$similarity.match) %in% R$antibody]
	if (length(R$similarity.match) == 0) R$similarity.match <- "None"
	if (R$warning.level <= 8) {
		warn <- "ok"
	} else if (R$warning.level > 8 & R$warning.level < 15) {
		warn <- "suboptimal"
	} else if (R$warning.level >= 15) {
		warn <- "unreliable"
	}
	if (verbose) {
		if (R$contact.level == 3 & length(R$exact.match) == length(x.ag)) {
			match <- "exact"
		} else if (R$contact.level == 3) {
			match <- "partial"
		} else {
			match <- "similarity"
		}
		message(paste0("###  Warning level: ", R$warning.level,
		               "/15 [", warn, "]"))
		message(paste0("###  Contact match: ", match))
		#message(paste0("###       Hotspots: ", R$hotspot))
	}
	return(R)
}

contacts <- function(x) {
	antibody = vector()
	variants = vector()
	residues = vector()
	class = vector()
	warningLevel = vector()
	contactLevel = vector()
	for (contact in x) {
		suppressMessages(R <- preprocess(contact, show.probs = FALSE,
		                                 verbose = FALSE))
		antibody <- c(antibody, R$antibody)
		variants <- c(variants, R$variant)
		residues <- c(residues, R$ab.residues)
		class <- c(class, R$class)
		warningLevel <- c(warningLevel, R$warning.level)
		contactLevel <- c(contactLevel, R$contact.level)
	}
	antibody <- unique(antibody)
	variants <- unique(variants)
	residues <- unique(residues)
	class <- unique(class)
	warningLevel <- unique(warningLevel)
	contactLevel <- unique(contactLevel)
	return(list(antibody = antibody, variant = variants,
	            ab.residues = residues, class = class,
	            warning.level = warningLevel, contact.level = contactLevel))
}

randomScalar <- function(data) {
	s <- apply(data, 1, sd)
	s <- unlist(lapply(s, function(x) runif(n = 1, min = -x, max = x)))
	return(s)
}

extractProfiles <- function(data = contact.data,
                            antibody = contact.data$antibody,
                            variants = contact.data$variants,
                            residues = NULL, stochastic = FALSE,
                            f.rep = "mean", f.res = "median",
                            digits = 2) {
	n <- nrow(data$affinity)
	Q <- data.frame(init = rep(NA, n))
	for (ab in unique(antibody)) {
		for (ag in unique(variants)) {
			A <- data$affinity[, data$antibody == ab & data$variant == ag]
			if (!is.null(residues)) {
				res <- unlist(lapply(names(A),
				              function(x) strsplit(x, "\\.\\d+")[[1]][1]))
				A <- A[, res %in% residues]
			}
			k <- ncol(A)/3
			if (k > 0) {
				M <- matrix(nrow = n, ncol = k)
				for (i in 1:k) {
					R <- A[, c(i, i + k, i + 2*k)]
					c <- paste0("M[, i] <- apply(R, 1, ", f.rep, ")")
					eval(parse(text = c))
					if (stochastic) {
						s <- randomScalar(R)
						M[, i] <- M[, i] + s
					}
				}
				c <- paste0("q <- apply(M, 1, ", f.res,")")
				w <- paste0("ab.", ab, ".", ag)
				eval(parse(text = c))
				q <- round(q, digits = digits)
				eval(parse(text = paste0("Q$", w, " <- q")))
			}
		}
	}
	return(Q[, -1])
}

mergeLast <- function(x) {
	if (length(x[[length(x)-1]]) > length(x[[length(x)]])) {
		x[[length(x)-1]] <- c(x[[length(x)-1]], x[[length(x)]])
		x <- x[-length(x)]
	}
	return(x)
}

truncSplit <- function(x) {
	if (length(x[[length(x)-1]]) > length(x[[length(x)]])) {
		x <- x[-length(x)]
	}
	return(x)
}

prepareLibrary <- function(data = contact.data, chunk = 5,
                           residues = NULL, stochastic = FALSE,
                           f.rep = "mean", f.res = "median", a0 = 0) {
	
	affinity <- data$affinity
	antibody <- data$antibody
	variants <- data$variants
	P <- extractProfiles(affinity, antibody, variants, residues,
	                     stochastic, f.rep, f.res)
	if (a0 > 0) {
		group <- ifelse(apply(P, 2, median) > a0, 0, 1)
		Q <- data.frame(a = vector(), y = vector(), res = vector(),
	                    antibody = vector(), variant = vector())
	} else {
		group <- NULL
		Q <- data.frame(a = vector(), y = vector(), res = vector(),
	                    antibody = vector(), variant = vector(),
	                    group = vector())
	}
	w <- unlist(strsplit(names(P), "\\."))
	w <- split(w, ceiling(seq_along(w)/3))
	
	for (i in 1:length(P)) {
		y <- P[, i]
		y <- split(y, ceiling(seq_along(y)/chunk))
		y <- truncSplit(y)
		
		A <- affinity[antibody == w[[i]][2] & variants == w[[i]][3]]
		
		for (j in 1:ncol(A)) {
			a <- split(A[, j], ceiling(seq_along(A[, j])/chunk))
			a <- truncSplit(a)
			n <- length(y)
			q <- data.frame(a = rep(NA, n), y = rep(NA, n))
			for (k in 1:n) {
				q$a[k] <- paste(a[[k]], collapse = ",")
				q$y[k] <- paste(y[[k]], collapse = ",")
			}
			q$res <- strsplit(names(A)[j], "\\.\\d+")[[1]]
			q$antibody <- w[[i]][2]
			q$variant <- w[[i]][3]
			if (a0 > 0) {
				contact <- paste0(w[[i]], collapse = ".")
				q$group <- as.numeric(group[names(group) == contact])
			}
			Q <- rbind(Q, q)
		}
	}
	return(Q)
}

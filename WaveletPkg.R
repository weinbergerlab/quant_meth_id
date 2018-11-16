### continuous wavelet transform functions

normalize.series <- function(x) {
	x <- x - mean(x)
	x <- x / sd(x)
}

cwt <- function(x, dt, dj, pad=1, s0=(-1), J=(-1), mother="morlet", param=(-1)) {
	x <- normalize.series(x)
	n <- length(x)
	if (s0 == -1) s0 <- 2 * dt
	if (J == -1) J <- floor(1/dj * log2(n * dt / s0))
	if (pad == 1) x <- c(x, rep(0, 2^(1+ceiling(log2(n)))-n)) # pad with zeros power to make 2 if desired
	padded.n <- length(x)
	k <- 0:(padded.n - 1)
	omega <- 2 * pi * k / (padded.n * dt)
	omega[k > padded.n/2] <- (-omega[k > padded.n/2]) # [Eqn.5]
	scales <- s0 * 2^((0:J) * dj) # construct scale vector
	fftx <- fft(x) # [Eqn(3)] compute FFT of the (padded) time series
	wavet <- matrix(0+0i, nrow=J + 1, ncol=padded.n)   # define the complex wavelet array
	for (i in 1:(J + 1)){ # loop through all scales and compute transform
	    w.base <- wavelet.bases(mother, omega, scales[i], param, dt)
		wavet[i, ] <- fft((fftx * w.base$daughter), inverse=T) / length(fftx) # wavelet transform[Eqn(4)]
	}
	periods <- w.base$fourier.factor * scales
	coi <- w.base$coi * dt * c(1E-5, 1:((n + 1) / 2 - 1), (n / 2 - 1):1, 1E-5) # COI [Sec.3g]
	wavet <- wavet[ , 1:n]       # get rid of padding before returning
	return(list("transform"=wavet, "periods"=periods, "scales"=scales, "coi"=coi))
}

wavelet.bases <- function(mother, omega, scale, param, dt) {
	mother <- tolower(mother)
	if (mother == "morlet") {
		if (param == -1) param <- 6
	   daughter <- pi^(-1/4) * as.numeric(omega > 0) * 
	   		exp(-(scale * omega-param)^2 / 2) # Table 1
		fourier.factor <- (4 * pi) / (param + sqrt(2 + param^2))
				# Scale-->Fourier [Sec.3h] Table 1 without scale
		coi <- fourier.factor / sqrt(2)  # Cone-of-influence [Sec.3g]
	} else if (mother == "paul") {
		if (param == -1) param <- 4
	   m <- param
		daughter <- 2^m / sqrt(m * factorial(2 * m - 1)) * 
				as.numeric(omega > 0) * (scale * omega)^m * exp(-scale * omega)
		fourier.factor <- 4 * pi / (2 * m + 1)
		coi <- fourier.factor * sqrt(2)
	} else if (mother == "dog") {
		if (param == -1) param <- 2
	   m <- param
		daughter <- (-1i^param) / sqrt(gamma(m + 1 / 2)) * 
				(scale * omega)^m * exp(-(scale * omega)^2 / 2)
		fourier.factor <- 2 * pi / sqrt(m + 1 / 2)
		coi <- fourier.factor / sqrt(2)
	} else stop("wavelet must be morlet, paul, or dog")
	daughter <- sqrt(2 * pi * scale / dt) * daughter        # normalize [Eqn.6]
	w.base <- list("daughter"=daughter, "fourier.factor"=fourier.factor, "coi"=coi)
	return(w.base)
}

global.mean <- function(mat, periods, coi) {
	global <- numeric(length(periods))
	for (i in 1:length(periods)) {
	    global[i] <- mean(mat[i, ][coi >= periods[i]])
	}
	return(global)
}

power.chisquare <- function(series, dt, dj, red.noise=T, sig=0.95, mother="morlet", param=6, ...) {
	if (tolower(mother) != "morlet" | param != 6) stop("this significance test currently only works for the morlet(6) wavelet")
	n <- length(series)
	param.cwt <- cwt(1:n, dt, dj, ...)
	scales <- param.cwt$scales
	periods <- param.cwt$periods
	na <- numeric(length(scales))
	for (i in 1:length(scales)) na[i] <- sum(param.cwt$coi > periods[i])
	if (red.noise) {
		if(ar(series, order.max=1)$order == 1) {
			alpha <- ar(series, order.max=1)$ar
		} else alpha <- 0
	}
	if (!red.noise) alpha <- 0
	gamma.fac <- 2.32  # for morlet(6)
	freq <- dt / periods        # freq=k/N for the scales used in the transform 
	fft.theor <- (1 - alpha^2) / (1 + alpha^2 - 2 * alpha * cos(freq * 2 * pi))   # [Eqn(16)] - gives vector P(k)
	local.sig.tmp <- qchisq(1 - (1 - sig) / 2, 2) * fft.theor / 2	
	local.sig <- matrix(rep(local.sig.tmp, n), nrow=length(periods), ncol=n)
	v <- 2 * sqrt(1 + (na * dt / (gamma.fac * scales))^2)     # [Eqn(23)]
	fft.theor[na == 0] <- v[na == 0] <- NA
	global.sig <- matrix(nrow=length(periods), ncol=3, dimnames=list(1:length(periods), c("lower", "mean", "upper")))
	global.sig[ , "upper"] <- qchisq(1 - (1 - sig) / 2, v) * fft.theor / v
	global.sig[ , "mean"] <- fft.theor
	global.sig[ , "lower"] <- qchisq((1 - sig) / 2, v) * fft.theor / v
	local.sig <- qchisq(1 - (1 - sig) / 2, 2) * fft.theor / 2
	return(list("local"=matrix(local.sig, nrow=length(scales), ncol=n), "global"=as.data.frame(global.sig)))
}

red.noise.sim <- function(ar.object) { # generates a red noise series vector of length n and form x[t]=alpha*x[t-1]+N(0, var)
	x <- numeric(ar.object$n.used)
	x[1] <- rnorm(1, sd=sqrt(ar.object$var.pred))
	for (i in 2:ar.object$n.used) {
		x[i] <- ar.object$ar * x[i - 1] + 
				rnorm(1, sd=sqrt(ar.object$var.pred))
	}
	return(x)
}

power.mc <- function(series, dt, dj, sim=1000, sig=.95, red.noise=T, ...) { # use original time series (not normalized)
	n <- length(series)
	param.cwt <- cwt(1:n, dt, dj, ...)
	scales <- param.cwt$scales
	periods <- param.cwt$periods
	coi <- param.cwt$coi
	na <- numeric(length(scales))
	for (i in 1:length(scales)) na[i] <- sum(coi > periods[i])
	mc.array <- array(dim=c(length(scales), n, sim))
	for (i in 1:sim) {
		if (red.noise) {
			if(ar(series, order.max=1)$order==1) {
				wavet <- cwt(red.noise.sim(ar(series, order.max=1)), dt, dj, ...)
			} else wavet <- cwt(rnorm(n), dt, dj, ...)
		}
		if (!red.noise) wavet <- cwt(rnorm(n), dt, dj, ...)
		mc.array[ , , i] <- (abs(wavet$transform))^2
	}
	global.mc <- matrix(nrow=length(scales), ncol=3, dimnames=list(1:length(scales), c("lower", "mean", "upper")))
	local.mc <- rep(NA, length(scales))
	for (i in 1:length(scales)) {
		if (na[i] > 0) {
			global.mc[i, "mean"] <- mean(apply(mc.array[i, coi > periods[i], ], 2, mean))
			global.mc[i, "upper"] <- 
					quantile(apply(mc.array[i, coi > periods[i], ], 2, mean), probs=1 - ((1 - sig) / 2), type=8)
			global.mc[i, "lower"] <- quantile(apply(mc.array[i, coi > periods[i], ], 2, mean), probs=(1 - sig) / 2, type=8)
			local.mc[i] <- quantile(mc.array[i, coi > periods[i], ], probs=1 - ((1 - sig) / 2), type=8)
		}
	}
	return(list("local"=matrix(local.mc, nrow=length(scales), ncol=n), "global"=as.data.frame(global.mc)))
}

complete.cwt <- function(series, dj=(-1), dt=(-1), mc.sim=-1, sig=0.95, ...) { ### use mc.sim=-1 for chi-square (non-simulated) significance
	if (dj == -1) stop("dj must be specified by the user")
	if (dt == -1) dt <- series$t[2] - series$t[1]
	wavet <- cwt(series$x, dt, dj, ...)
	w.power <- (abs(wavet$transform))^2
	global <- global.mean(w.power, wavet$periods, wavet$coi)
	if (mc.sim > 0) power.sig <- power.mc(series$x, dt, dj, sim=mc.sim, sig=sig, ...)
	if (mc.sim == -1) power.sig <- power.chisquare(series$x, dt, dj, sig=sig, ...)
	return(list("transform"=wavet$transform, "local"=w.power, "local.sig"=power.sig$local, "global"=global, "global.sig"=power.sig$global, "time"=series$t, "periods"=wavet$periods, "scales"=wavet$scales, "coi"=wavet$coi, "x"=series$x, "dj"=dj, "dt"=dt))
}

### wavelet coherence functions

smooth.cwt <- function(wavet, dt, dj, scales) {
	n <- dim(wavet)[2]
	npad <- 2^ceiling(log2(n))	
	k <- c(0, 1:floor(npad / 2), (-ceiling(npad / 2 - 1)):-1)
	snorm <- scales / dt
	for (i in 1:dim(wavet)[1]) {
		time.filter <- exp(-k^2 / (2 * snorm[i]^2))
		time.filter <- time.filter / sum(time.filter) # normalize
		smooth <- convolve(time.filter, c(wavet[i, ], rep(0, npad - n)), 
				conj=F) # filter in time
		wavet[i, ] <- smooth[1:n]
	}
	dj0 <- 0.60
	sf.width <- dj0 / (dj * 2)
	scale.filter <- c(sf.width - floor(sf.width), 
			rep(1, round(2 * sf.width - 1)), sf.width - floor(sf.width))
	scale.filter <- scale.filter / sum(scale.filter) # normalize
	for (i in 1:dim(wavet)[2]) {
		conv.tmp <- convolve(wavet[ , i], scale.filter, 
				type='open') # filter in scale
		pad <- (length(conv.tmp) - length(scales)) / 2
	   wavet[ , i] <- conv.tmp[(pad + 1):(pad + length(scales))]
	}
	return(wavet)
}

normalize.E <- function(mat, scales) {
	out <- array(dim=dim(mat))
	for (i in 1:dim(out)[2]) {
		out[ , i] <- scales^(-1) * mat[ , i] # normalize
	}
	return(out)
}

coherence.calc <- function(cwt1, cwt2, dt, dj, start.t, 
		end.t, series1.t, series2.t) {
	scales <- cwt1$scales
	T1 <- start.t <= series1.t & series1.t <= end.t
	T2 <- start.t <= series2.t & series2.t <= end.t
	Wxy <- cwt1$transform[ , T1] * Conj(cwt2$transform[ , T2])
	nWxy <- array(dim=dim(Wxy))
	for (i in 1:dim(Wxy)[2]) {
		nWxy[ , i] <- cwt1$scales^(-1) * Wxy[ , i] # normalize
	}
	numer <- abs(smooth.cwt(nWxy, dt, dj, cwt1$scales))^2
	denom <- smooth.cwt(normalize.E(abs(cwt1$transform)^2, scales), 
			dt, dj, scales)[ , T1] * 
			smooth.cwt(normalize.E(abs(cwt2$transform)^2, scales), 
			dt, dj, scales)[ , T2]
	return(numer / denom)
}

coherence.mc <- function(n, x1, x2, dt, dj, J, sim=1000, sig=.95, red.noise=F, ...) {
	param.cwt <- cwt(1:n, dt, dj, J=J, ...)
	periods <- param.cwt$periods
	coi <- param.cwt$coi
	global.mat <- local.mat <- matrix(nrow=J+1, ncol=sim)
	for (i in 1:sim) { # generate wavelet transforms for red noise time series
		if (red.noise) {
			cwt1 <- cwt(red.noise.sim(ar(x1, order.max=1))[1:n], dt, dj, J=J, ...)
			cwt2 <- cwt(red.noise.sim(ar(x2, order.max=1))[1:n], dt, dj, J=J, ...)
      } else {
			cwt1 <- cwt(rnorm(n), dt, dj, J=J, ...)
			cwt2 <- cwt(rnorm(n), dt, dj, J=J, ...)
		}
		coh.tmp <- coherence.calc(cwt1, cwt2, dt, dj, 1, n, 1:n, 1:n)
		for (i2 in 1:(J + 1)) {
			global.mat[i2, i] <- mean(coh.tmp[i2, coi > periods[i2]])
		}
		local.mat[ , i] <- coh.tmp[ , round(n / 2)]
	}
	global.mc <- matrix(nrow=J + 1, ncol=3, dimnames=list(1:(J+1), c("lower", "mean", "upper")))
	local.mc <- rep(NA, J + 1)
	for (i in 1:(J + 1)) {
		if (sum(coi > periods[i]) > 0) {
			global.mc[i, "mean"] <- mean(global.mat[i, ])
			global.mc[i, c("upper", "lower")] <- quantile(global.mat[i, ], 
					probs=c(1 - ((1 - sig) / 2), (1 - sig) / 2), type=8)
			local.mc[i] <- quantile(local.mat[i, ], probs=1 - ((1 - sig) / 2), type=8)
		}
	}
	return(list("local"=matrix(local.mc, nrow=J+1, ncol=n), "global"=as.data.frame(global.mc)))
}

complete.coh <- function(series1, series2, dj=(-1), mc.sim=1000, dt=(-1), s0=(-1), J=(-1), sig=0.95, ...) {
	if (dj == -1) stop("dj must be specified by the user")
	if (dt == -1) dt <- series1$t[2] - series1$t[1]
	start.t <- max(c(min(series1$t), min(series2$t)))
	end.t <- min(c(max(series1$t), max(series2$t)))
	shared.t <- seq(start.t, end.t, by=dt)
	shared.cwt <- cwt(shared.t, dt, dj, J=J, ...)
	if (s0 == -1) s0 <- 2 * dt
	if (J == -1) J <- floor(1 / dj * log2(length(shared.t) * dt / s0)) # sets the scales of the transforms
	cwt1 <- cwt(series1$x, dt, dj, J=J, ...)
	cwt2 <- cwt(series2$x, dt, dj, J=J, ...)
	coherence <- coherence.calc(cwt1, cwt2, dt, dj, start.t, end.t, series1$t, series2$t)
	global <- global.mean(coherence, shared.cwt$periods, shared.cwt$coi)
	coherence.sig <- coherence.mc(length(shared.t), series1$x, series2$x, dt, dj, J, sim=mc.sim, sig=sig, ...)
	coherence.sig$local <- coherence.sig$local[ , 1:length(shared.t)]
	return(list("local"=coherence, "local.sig"=coherence.sig$local, 
	"global"=global, "global.sig"=coherence.sig$global, "time"=shared.t, 
	"periods"=shared.cwt$periods, "coi"=shared.cwt$coi, "dj"=dj, "dt"=dt))
}

### additional functions

reconstruct.series <- function(cwt.object, sp1=-1, sp2=-1, 
		mother="morlet", param=6, ...) { 
		# sp1 and sp2 are the range of periods for which reconstruction is desired
	if (mother != "morlet" | param != 6) {
		stop("reconstruction is only setup for morlet(6) wavelets")
	}
	if (sp1 == -1 & sp2 == -1) {
		sp1 <- min(cwt.object$periods)
		sp2 <- max(cwt.object$periods)
	}
	recon <- numeric(length(cwt.object$time))
	dt <- cwt.object$time[2] - cwt.object$time[1]
	Cdelta <- 0.776                    # for the MORLET wavelet
	selected <- sp1 <= cwt.object$periods & cwt.object$periods <= sp2
	for (i in 1:length(recon)) {
		recon[i] <- (cwt.object$dj * sqrt(dt) / (Cdelta * (pi^(-1 / 4)))) * 
		sum(Re(cwt.object$transform[selected, i]) / sqrt(cwt.object$scales[selected]))
	}
	return(recon)
}

phase <- function(cwt.object, s=NA) {
	if (length(s) == 1) {
		if (!is.na(s)) selected <- s # 1 scale
	} else if (length(s) == 2) {
		selected <- (1:length(cwt.object$periods))[
				s[1] <= cwt.object$periods & cwt.object$periods <= s[2]] # set of scales
	} else stop("ERROR: scale is not specified correctly")
	out <- numeric(length(cwt.object$x))
	for (i in 1:length(out)) {
			out[i] <- atan2(sum(Im(cwt.object$transform[selected, i])), 
					sum(Re(cwt.object$transform[selected, i])))
	}
	return(out)
}

cohPhaseDiff <- function(series1, series2, dj=(-1), 
		s=NA, dt=(-1), s0=(-1), J=(-1), ...) { # Torrence 1999 Eqn. A3
	if (dj == -1) stop("dj must be specified by the user")
	if (dt == -1) dt <- series1$t[2] - series1$t[1]
	startT <- max(c(min(series1$t), min(series2$t)))
	endT <- min(c(max(series1$t), max(series2$t)))
	T1 <- startT <= series1$t & series1$t <= endT
	T2 <- startT <= series2$t & series2$t <= endT
	if (s0 == -1) s0 <- 2 * dt
	if (J == -1) J <- floor(1 / dj * log2(sum(T1) * dt / s0))
	cwt1 <- cwt(series1$x, dt, dj, J=J, ...)
	cwt2 <- cwt(series2$x, dt, dj, J=J, ...)
	Wxy <- cwt1$transform[ , T1] * Conj(cwt2$transform[ , T2])
	sWxy <- smooth.cwt(normalize.E(Wxy, cwt1$scales), 
			dt, dj, cwt1$scales)
	phase.coh <- atan2(Im(sWxy), Re(sWxy))	
	if (length(s) == 1) {
		if (!is.na(s)) {
			out <- phase.coh[s, ]
		}
	} else if (length(s) == 2) {
		out <- apply(phase.coh[s[1] <= cwt1$periods & cwt1$periods <= s[2], ], 
				2, mean)
	} else stop("ERROR: scale is not specified correctly")
	return(out)
}

phaseToTime <- function(phase, p, dt) {
	phase / (2 * pi) * p * 1/dt
}

### local length calcs
local.length <- function(coh, p1, p2) {
	local.sig <- coh$local / coh$local.sig
	local.sig[local.sig >= 1] <- 1
	local.sig[local.sig < 1] <- 0
	for (i in 1:length(coh$periods)) local.sig[i, coh$coi < max(coh$periods[i])] <- 0
	selected <- (1:length(coh$periods))[p1 <= coh$periods & coh$periods <= p2]
	max.length <- rep(NA, length(coh$periods))
	for (i in selected) {
		max.length[i] <- max.tmp <- last.tmp <- 0
		for (i2 in 1:dim(local.sig)[2]) {
			if (last.tmp == 1) {
				if (local.sig[i, i2] == 1) {
					max.tmp <- max.tmp + 1
				}
				if (local.sig[i, i2] == 0) {
					if (max.tmp > max.length[i]) max.length[i] <- max.tmp
					last.tmp <- 0
				}
			}
			if (last.tmp == 0) {
				if (local.sig[i, i2] == 1) {
					last.tmp <- 1
				}
			}
		}
	}
	return(max.length)
}

local.length.sim <- function(coh, sim=1000, p1, p2, s0=(-1), J=(-1), ...) {
	dj <- coh$dj
	dt <- coh$time[2]-coh$time[1]
	if (s0==-1) s0 <- 2*dt
	if (J==-1) J <- floor(1/dj*log2(length(coh$time)*dt/s0))
	n <- length(coh$time)
	max.length <- matrix(NA, nrow=length(coh$periods), ncol=sim)
	coh.tmp <- coh # this preserves the significance calculation for the original coherence
	for (i in 1:sim) {
		series1 <- rnorm(n)
		series2 <- rnorm(n)
		cwt1 <- cwt(series1, dt, dj)#, ...)
		cwt2 <- cwt(series2, dt, dj)#, ...)
		coh.tmp$local <- coherence.calc(cwt1, cwt2, dt, dj, 1, n, 1:n, 1:n)
		max.length[ , i] <- local.length(coh.tmp, p1, p2)
	}
	return(max.length)
}

local.length.ci <- function(simulations, sig=.95, ...) {
    global.length <- matrix(NA, nrow=dim(simulations)[1], ncol=3, dimnames=list(1:dim(simulations)[1], c("lower", "mean", "upper")))
    for (i in 1:dim(simulations)[1]) {
        if (!is.na(mean(simulations[i, ]))) {
            global.length[i, ] <- quantile(simulations[i, ], probs=c((1-sig)/2, 0.5, 1-((1-sig)/2)), type=1)
        }
    }
    return(as.data.frame(global.length))
}

local.length.p <- function(max.length, simulations) {
    global.p <- rep(NA, dim(simulations)[1])
    for (i in 1:dim(simulations)[1]) {
        if (!is.na(mean(simulations[i, ]))) global.p[i] <- sum(simulations[i, ]<=max.length[i])/dim(simulations)[2]
    }
    return(round(global.p, 5))
}

local.length.composite <- function(coh, sim, p1, p2, sig=.95, ...) {
	local.length.tmp <- local.length(coh, p1, p2)
	local.length.sim.tmp <- local.length.sim(coh, sim, p1, p2, ...)
	local.length.ci.tmp <- local.length.ci(local.length.sim.tmp, ...)
	local.length.p.tmp <- local.length.p(local.length.tmp, local.length.sim.tmp)
	return(cbind("periods"=coh$periods, "actual"=local.length.tmp, "p"= local.length.p.tmp, local.length.ci.tmp))
}


####################################################################################
### plots
tim.64 <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", 
        "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", 
        "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", 
        "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", 
        "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", 
        "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
        "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", 
        "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", 
        "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", 
        "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", 
        "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", 
        "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", 
        "#AF0000", "#9F0000", "#8F0000", "#800000")

tim.12 <- c("#00008F", "#0000DF", "#0030FF", "#0080FF", "#00CFFF", "#20FFDF", "#70FF8F",
 			"#BFFF40", "#FFEF00", "#FF9F00", "#FF5000", "#FF0000", "#AF0000")
 
plot.cwt <- function(cwt.obj) {
	split.screen(rbind(c(0, 0.75, 0, 1), c(0.75, 1, 0, 1)))
	
	screen(1)
		par(mar=c(3, 4, 1, .5))
		image(cwt.obj$time, log2(cwt.obj$periods), t(cwt.obj$local), col=tim.64, breaks=quantile(cwt.obj$local, probs=seq(0, 1, by=1/64)), axes=F, xlab="", ylab="Period (years)", xlim=(range(cwt.obj$time)), ylim=range(log2(cwt.obj$periods)))
		axis(side=1, at=round(cwt.obj$time)+.5, labels=F)
		axis(side=1, at=seq(ceiling(min(cwt.obj$time)), floor(max(cwt.obj$time)), by=floor(length(cwt.obj$time)*cwt.obj$dt/5))+.5, labels=seq(ceiling(min(cwt.obj$time)), floor(max(cwt.obj$time)), by=floor(length(cwt.obj$time)*cwt.obj$dt/5)), tcl=-.8)
		axis(side=2, at=(-1:floor(log2(max(cwt.obj$periods)))), labels=2^(-1:floor(log2(max(cwt.obj$periods)))))
		contour(cwt.obj$time, log2(cwt.obj$periods), t(cwt.obj$local/cwt.obj$local.sig), add=T, levels=c(1), drawlabels=F, col="black", lwd=2, lty=1)	
		polygon(c(min(cwt.obj$time), min(cwt.obj$time), cwt.obj$time, max(cwt.obj$time), max(cwt.obj$time)), c(max(cwt.obj$periods), min(log2(cwt.obj$coi)),  log2(cwt.obj$coi), min(log2(cwt.obj$coi)), max(cwt.obj$periods)), border=NA, col=hsv(0, 0, .5, alpha=.5))
		box()
	
	screen(2)
		par(mar=c(3, .5, 1, 1))
		ci<-cbind(x=c(cwt.obj$global.sig$upper[!is.na(cwt.obj$global.sig$upper)], rev(cwt.obj$global.sig$lower[!is.na(cwt.obj$global.sig$lower)])), y=log2(c(cwt.obj$periods[!is.na(cwt.obj$global.sig$upper)], rev(cwt.obj$periods[!is.na(cwt.obj$global.sig$lower)]))))
		plot(cwt.obj$global, log2(cwt.obj$periods), type='l', lwd=2, xlab="", axes=F, xaxs='i', yaxs='i', xlim=c(0, 2*max(cwt.obj$global, na.rm=T)), ylab="", ylim=range(log2(cwt.obj$periods)), main="", lty=2)
		polygon(ci, border='dark grey', col='light grey')
		lines(cwt.obj$global.sig$mean, log2(cwt.obj$periods), col="dark grey", lwd=2)
		lines(cwt.obj$global, log2(cwt.obj$periods), col='black', lwd=2, lty=1)
		axis(side=1)
		box()
		
	close.screen(all = TRUE)
}


plot.coh <- function(coh.obj) {
	split.screen(rbind(c(0, 0.75, 0, 1), c(0.75, 1, 0, 1)))
	
	screen(1)
		par(mar=c(3, 4, 1, .5))
		image(coh.obj$time, log2(coh.obj$periods), t(coh.obj$local), 
				col=tim.64, breaks=quantile(coh.obj$local, 
				probs=seq(0, 1, by=1/64)), axes=F, 
				xlab="", ylab="Period (years)", 
				xlim=(range(coh.obj$time)), 
				ylim=range(log2(coh.obj$periods)))
		axis(side=1, at=round(coh.obj$time)+.5, labels=F)
		axis(side=1, at=seq(ceiling(min(coh.obj$time)), 
				floor(max(coh.obj$time)), 
				by=floor(length(coh.obj$time)*coh.obj$dt/5))+.5, 
				labels=seq(ceiling(min(coh.obj$time)), 
				floor(max(coh.obj$time)), 
				by=floor(length(coh.obj$time)*coh.obj$dt/5)), tcl=-.8)
		axis(side=2, at=(-1:floor(log2(max(coh.obj$periods)))), 
				labels=2^(-1:floor(log2(max(coh.obj$periods)))))
		contour(coh.obj$time, log2(coh.obj$periods), 
				t(coh.obj$local/coh.obj$local.sig), add=T, 
				levels=c(1), drawlabels=F, col="black", lwd=2, lty=1)	
		polygon(c(min(coh.obj$time), min(coh.obj$time), coh.obj$time, 
				max(coh.obj$time), max(coh.obj$time)), c(max(coh.obj$periods), 
				min(log2(coh.obj$coi)),  log2(coh.obj$coi), 
				min(log2(coh.obj$coi)), max(coh.obj$periods)), 
				border=NA, col=hsv(0, 0, .5, alpha=.5))
		box()
	
	screen(2)
		par(mar=c(3, .5, 1, 1))
		ci<-cbind(x=c(coh.obj$global.sig$upper[!is.na(coh.obj$global.sig$upper)], 
				rev(coh.obj$global.sig$lower[!is.na(coh.obj$global.sig$lower)])), 
				y=log2(c(coh.obj$periods[!is.na(coh.obj$global.sig$upper)], 
				rev(coh.obj$periods[!is.na(coh.obj$global.sig$lower)]))))
		plot(coh.obj$global, log2(coh.obj$periods), type='l', 
				lwd=2, xlab="", axes=F, xaxs='i', yaxs='i', 
				xlim=c(0, 2*max(coh.obj$global, na.rm=T)), 
				ylab="", ylim=range(log2(coh.obj$periods)), 
				main="", lty=2)
		polygon(ci, border='dark grey', col='light grey')
		lines(coh.obj$global.sig$mean, log2(coh.obj$periods), 
				col="dark grey", lwd=2)
		lines(coh.obj$global, log2(coh.obj$periods), 
				col='black', lwd=2, lty=1)
		axis(side=1)
		box()
		
	close.screen(all = TRUE)
}




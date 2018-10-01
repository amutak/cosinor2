#' Acrophase Correction
#'
#' @encoding UTF-8
#' @param x An object of the \code{cosinor.lm} class.
#' @description Corrects the value of the acrophase parameter of the cosinor model, placing it in the appropriate quadrant.
#' @details The acrophase parameter of a cosinor model is found by solving an equation with inverse tangent of an expression which contains linearized cosinor parameters. However, multiple numeric entities may result in a same value of tangent and just calculating the inverse tangent may result with the wrong value of the acrophase. This function corrects the acrophase from the \code{cosinor.lm} object according to the procedure from Bingham et al. (1982).
#' More specifically, the acrophase is calculated as: \deqn{K + g * arctan \vert\frac{\gamma}{\beta}\vert}
#' where values of \eqn{K} and \eqn{g} depend on the signs of \eqn{\beta} and \eqn{\gamma} and can be derived from the following table:
#' \tabular{cccc}{
#' sign \eqn{\beta} \tab sign \eqn{\gamma} \tab K \tab g \cr
#' + \tab + \tab 0 \tab -1 \cr
#' + \tab - \tab -2\eqn{\pi} \tab 1 \cr
#' - \tab + \tab -\eqn{\pi} \tab 1 \cr
#' - \tab - \tab -\eqn{\pi} \tab -1
#' }
#' @examples
#' fit.temperature<-cosinor.lm(Temperature~time(Time), period = 24, data = temperature_zg)
#' correct.acrophase(fit.temperature)
#' @references Bingham, C., Arbogast, B., Guillaume Cornélissen, G., Lee, J.K. & Halberg, F. (1982). Inferential Statistical Methods for Estimating and Comparing Cosinor Parameters. \emph{Chronobiologia}, \emph{9(4)}, 397-439.
#' @export

correct.acrophase<-function(x){
  coefs<-data.frame(t(x$fit$coefficients))
  if(coefs$rrr>0 & coefs$sss>0){
    acrophase<-0+(-1*atan(abs(coefs$sss/coefs$rrr)))
  }
  else if(coefs$rrr>0 & coefs$sss<0){
    acrophase<-2*-1*pi+(1*atan(abs(coefs$sss/coefs$rrr)))
  }
  else if(coefs$rrr<0 & coefs$sss>0){
    acrophase<-pi*-1+(1*atan(abs(coefs$sss/coefs$rrr)))
  }
  else {
    acrophase<-pi*-1+(-1*atan(abs(coefs$sss/coefs$rrr)))
  }
  return(acrophase)
}

#' Percent Rhythm
#'
#' @param x An object of the \code{cosinor.lm} or \code{population.cosinor.lm} class.
#' @description Calculates Percent Rhythm, the measure of the relative strength of a rhythm.
#' @details Percent Rhythm is the coefficient of determination obtained by squaring the correlation between observed and estimated values.
#' @examples
#' fit.temperature<-cosinor.lm(Temperature~time(Time), period = 24, data = temperature_zg)
#' cosinor.PR(fit.temperature)
#'
#' fit.november<-population.cosinor.lm(data = PANAS_november, time = PANAS_time,
#' period = 7)
#' cosinor.PR(fit.november)
#' @importFrom Hmisc rcorr
#' @importFrom graphics par plot
#' @importFrom stats cor cov lm pf sd qt as.formula na.omit
#' @export

cosinor.PR<-function(x){
  if (class(x) == "population.cosinor.lm") {
    corr<-rcorr(x$fitted.values,x$emp.mean)
    PR<-data.frame(cbind(corr$r[[2]],corr$r[[2]]^2,corr$P[[2]]))
    colnames(PR)<-c("r","Percent rhythm","p-value")
  }
  else {
    corr<-rcorr(x$fit$fitted.values,x$fit$model[,1])
    PR<-data.frame(cbind(corr$r[[2]],corr$r[[2]]^2,corr$P[[2]]))
    colnames(PR)<-c("r","Percent rhythm","p-value")
  }
  return(PR)
}

#' Rhythm Detection Test
#'
#' @param x An object of the \code{cosinor.lm} or \code{population.cosinor.lm} class.
#' @description Performs the rhythm detection test, a global test for the significance of the estimated model for single cosinor and population-mean cosinor.
#' @details The rhythm detection test, also called the zero-amplitude test, tests the overall significance of the cosinor model. The test is actually an F-ratio and is calculated as following (according to the procedure described in Cornélissen, 2014):
#' \deqn{F = \frac{\frac{\sum(\widehat{Y}_i - \bar{Y})^2}{2}}{\frac{\sum(Y_i - \bar{Y})^2}{N - 3}}} with \eqn{df_1 = 2} and \eqn{df_2 = N - 3}
#' , where \eqn{\widehat{Y}_i} is the \eqn{i}th estimated value, \eqn{Y_i} is the \eqn{i}th observed value, \eqn{\bar{Y}} is the arithmetic mean of observed values and \eqn{N} is the number of timepoints.
#' For the population-mean cosinor model, the test is calculated according to the procedure described in Bingham et al. (1982) as follows:
#' \deqn{F = \frac{k(k-2)}{2(k-1)}\frac{1}{1-(\frac{\widehat{\sigma}_{\beta\gamma}}{\widehat{\sigma}_\beta \widehat{\sigma}_\gamma})^2}[\frac{\beta^2}{\widehat{\sigma}^2_\beta}-2\frac{\widehat{\sigma}_{\beta\gamma}}{\widehat{\sigma}_\beta \widehat{\sigma}_\gamma}\frac{\beta \gamma}{\widehat{\sigma}_\beta \widehat{\sigma}_\gamma}+\frac{\gamma^2}{\widehat{\sigma}^2_\gamma}]} with \eqn{df_1 = 2} and \eqn{df_2 = k - 2}
#' , where \eqn{k} is the number of subjects in the population, \eqn{\widehat{\sigma}_\beta} and \eqn{\widehat{\sigma}_\gamma} are standard deviations of population \eqn{\beta} and \eqn{\gamma} coefficients and \eqn{\widehat{\sigma}_{\beta\gamma}} is the covariance of population \eqn{\beta} and \eqn{\gamma} coefficients.
#' @examples
#' fit.temperature<-cosinor.lm(Temperature~time(Time), period = 24, data = temperature_zg)
#' cosinor.detect(fit.temperature)
#'
#' fit.november<-population.cosinor.lm(data = PANAS_november, time = PANAS_time,
#' period = 7)
#' cosinor.detect(fit.november)
#' @references Cornélissen, G. (2014). Cosinor-Based Rhythmometry. \emph{Theoretical Biology and Medical Modeling}, \emph{11}, Article 16.
#'
#' Bingham, C., Arbogast, B., Guillaume Cornélissen, G., Lee, J.K. & Halberg, F. (1982). Inferential Statistical Methods for Estimating and Comparing Cosinor Parameters. \emph{Chronobiologia}, \emph{9(4)}, 397-439.
#' @importFrom stats cor cov lm pf sd qt as.formula na.omit
#' @export

cosinor.detect<-function(x){
  if(class(x) == "population.cosinor.lm"){
    betas<-t(x$pop.mat[2,])
    gammas<-t(x$pop.mat[3,])
    beta<-mean(betas)
    gamma<-mean(gammas)
    sdb<-(sd(betas))
    sdy<-(sd(gammas))
    r<-cor(betas,gammas)
    k=ncol(x$pop.mat)
    frac1<-(k*(k-2))/(2*(k-1))
    frac2<-1/(1-r^2)
    frac3<-beta^2/sdb^2
    frac4<-(beta*gamma)/(sdb*sdy)
    frac5<-gamma^2/sdy^2
    brack<-frac3-(2*r*frac4)+frac5
    Fvalue<-frac1*frac2*brack
    df2<-k-2
    pvalue<-pf(q=Fvalue,df1=2,df2=df2,lower.tail=F)
    detection_test<-cbind(Fvalue,2,k-2,pvalue)
    colnames(detection_test)<-c("F","df1","df2","p")
    rownames(detection_test)<-NULL
  }
  else{
    RSS<-sum((x$fit$residuals)^2)
    MSS<-sum((x$fit$fitted.values-mean(x$fit$fitted.values+x$fit$residuals))^2)
    Fvalue<-(MSS/2)/(RSS/x$fit$df.residual)
    pvalue<-pf(q=Fvalue,df1=2,df2=x$fit$df.residual,lower.tail=F)
    detection_test<-cbind(Fvalue,2,x$fit$df.residual,pvalue)
    colnames(detection_test)<-c("F","df1","df2","p")
  }
  return(detection_test)
}

#' Population-Mean Cosinor
#'
#' @param data A data frame containing responses of subjects collected over time, with subjects in the rows and timepoints in the columns.
#' @param time A vector containing the times at which the data was collected.
#' @param period Duration of one cycle of rhythm.
#' @param na.action Action to be performed on missing values. Defaults to \code{na.omit}.
#' @param alpha Significance level for calculating population cosinor parameters confidence intervals. Defaults to .05 (confidence intervals are 5\% risk intervals).
#' @param plot Logical, display plot after calculation? Defaults to \code{TRUE}.
#' @description Calculates the population-mean cosinor.
#' @note If the confidence interval of the population amplitude includes zero, confidence interval of the acrophase cannot be calculated reliably. If this case occurs while using this function, the user will be warned and acrophase confidence interval limits will be set to NA.
#' @details According to the procedure described in Cornélissen (2014), to calculate population-mean cosinor, single cosinors are first performed on each subject and linearized parameters are averaged, which allows for calculation of delinearized parameters. After such a procedure is completed, confidence intervals of population-mean cosinor parameters can be calculated as described in Bingham et al. (1982) using the following formulae:
#' \deqn{\widehat{M} \pm \frac{t_{1-\frac{\alpha}{2}}\widehat{\sigma}_M}{\sqrt{k}}}
#' \deqn{\widehat{\phi}+arctan(\frac{c_{23} t_{1-\frac{\alpha}{2}}^2 \pm t_{1-\frac{\alpha}{2}}\sqrt{c_{33}} \sqrt{\widehat{A}^2-\frac{(c_{22}c_{33}-c_{23}^2)t_{1-\frac{\alpha}{2}}^2}{c_{33}}}}{\widehat{A}^2 - c_{22} t_{1-\frac{\alpha}{2}}^2})}
#' \deqn{\widehat{A} \pm t_{1-\frac{\alpha}{2}} \sqrt{c_{22}}}
#' where \eqn{c_{22}}, \eqn{c_{23}} and \eqn{c_{33}} are elements of the sampling scheme matrix, calculated as follows:
#' \deqn{c_{22}=\frac{\widehat{\sigma}^2_{\beta}\widehat{\beta}^2+2\widehat{\sigma}_{\beta \gamma}\widehat{\beta}\widehat{\gamma}+\widehat{\sigma}^2_{\gamma}\widehat{\gamma}^2}{k\widehat{A}^2}}
#' \deqn{c_{23}=\frac{-(\widehat{\sigma}^2_{\beta}-\widehat{\sigma}^2_{\gamma})(\widehat{\beta}\widehat{\gamma})+\widehat{\sigma}_{\beta \gamma}(\widehat{\beta}^2-\widehat{\gamma}^2)}{k\widehat{A}^2}}
#' \deqn{c_{33}=\frac{\widehat{\sigma}^2_{\beta}\widehat{\gamma}^2-2\widehat{\sigma}_{\beta \gamma}\widehat{\beta}\widehat{\gamma}+\widehat{\sigma}^2_{\gamma}\widehat{\beta}^2}{k\widehat{A}^2}}
#' where \eqn{\widehat{M}}, \eqn{\widehat{A}}, \eqn{\widehat{\phi}}, \eqn{\widehat{\beta}} and \eqn{\widehat{\gamma}} are population-mean cosinor parameters, \eqn{\widehat{\sigma}_M}, \eqn{\widehat{\sigma}_{\beta}} and \eqn{\widehat{\sigma}_{\gamma}} are the standard deviations of the single cosinor parameters, \eqn{\widehat{\sigma}_{\beta \gamma}} is the covariance of the single cosinor \eqn{\beta} and \eqn{\gamma} coefficients, \eqn{k} is the number of subjects in a population and \eqn{t_{1-\frac{\alpha}{2}}} is the two-tailed inverse of the t-distribution with \eqn{\alpha} level of significance and \eqn{k - 1} degrees of freedom.
#' @examples
#' population.cosinor.lm(data = PANAS_november, time = PANAS_time,
#' period = 7, na.action = "na.exclude")
#' @references Cornélissen, G. (2014). Cosinor-Based Rhythmometry. \emph{Theoretical Biology and Medical Modeling}, \emph{11}, Article 16.
#'
#' Bingham, C., Arbogast, B., Guillaume Cornélissen, G., Lee, J.K. & Halberg, F. (1982). Inferential Statistical Methods for Estimating and Comparing Cosinor Parameters. \emph{Chronobiologia}, \emph{9(4)}, 397-439.
#' @return Object of the \code{population.cosinor.lm} class containing the following objects:
#'  \item{\code{single.cos}}{A list of objects containing all performed single cosinors.}
#'  \item{\code{pop.mat}}{A data frame containing the cosinor parameters of each subject in the population.}
#'  \item{\code{coefficients}}{Delinearized population-mean cosinor coefficients.}
#'  \item{\code{emp.mean}}{Empirical mean of the data across all timepoints.}
#'  \item{\code{fitted.values}}{Estimated values of the rhythm caclculated using the cosinor model.}
#'  \item{\code{residuals}}{The difference between empirical mean and the fitted values.}
#'  \item{\code{conf.int}}{Values of upper and lower limits of confidence intervals of delinearized cosinor parameters.}
#' @import cosinor ggplot2 matrixStats Hmisc
#' @importFrom magrittr %>%
#' @importFrom stringr str_c
#' @importFrom purrr map
#' @importFrom graphics par plot
#' @importFrom stats cor cov lm pf sd qt as.formula na.omit
#' @export

population.cosinor.lm<-function(data,time,period,na.action = na.omit,alpha=.05, plot = T){
  population.cosinor<-list()
  popmat<-data.frame(matrix(nrow=3,ncol=nrow(data)))
  coefficients<-data.frame(matrix(nrow=1,ncol=3))
  colnames(coefficients)<-c("MESOR","Amplitude","Acrophase")
  subject<-integer()
  amplitude<-integer()
  acrophase<-integer()
  fitted.values<-vector()
  residuals<-vector()
  data<-data.frame(cbind(t(data),time))
  colnames(data)<-c(paste0("Subj",1:(ncol(data)-1)), "time")
  cosinors<-str_c(paste0("Subj",1:(ncol(data)-1)), "time(time)", sep="~" ) %>%
    map(as.formula) %>%
    map(cosinor.lm, data = data, period = period, na.action = na.action)
  for (subject in 1:(ncol(data)-1)) {
    popmat[[subject]]<-cosinors[[subject]]$fit$coefficients
  }
  population.cosinor[[1]]<-cosinors
  population.cosinor[[2]]<-popmat
  names(population.cosinor)<-c("single.cos","pop.mat")
  coefs<-rowMeans(popmat)
  MESOR<-coefs[[1]]
  beta<-coefs[[2]]
  gamma<-coefs[[3]]
  amplitude<-sqrt(beta^2+gamma^2)
  if(beta>0 & gamma>0){
    acrophase<-0+(-1*atan(abs(gamma/beta)))
  }
  else if(beta>0 & gamma<0){
    acrophase<-2*-1*pi+(1*atan(abs(gamma/beta)))
  }
  else if(beta<0 & gamma>0){
    acrophase<-pi*-1+(1*atan(abs(gamma/beta)))
  }
  else {
    acrophase<-pi*-1+(-1*atan(abs(gamma/beta)))
  }
  coefficients[,1]<-MESOR
  coefficients[,2]<-amplitude
  coefficients[,3]<-acrophase
  population.cosinor[[3]]<-coefficients
  sdm<-(sd(popmat[1,]))
  sdb<-(sd(popmat[2,]))
  sdy<-(sd(popmat[3,]))
  covby<-(cov(t(popmat[2,]), t(popmat[3,])))
  k<-ncol(popmat)
  denom<-(amplitude^2)*k
  c22<-(((sdb^2)*(beta^2))+(2*covby*beta*gamma)+((sdy^2)*(gamma^2)))/denom
  c23<-(((-1*((sdb^2)-(sdy^2)))*(beta*gamma))+(covby*((beta^2)-(gamma^2))))/denom
  c33<-(((sdb^2)*(gamma^2))-(2*covby*beta*gamma)+((sdy^2)*(beta^2)))/denom
  t<-abs(qt(alpha/2, df = k-1))
  mesoru<-MESOR+((t*sdm)/sqrt(k))
  mesorl<-MESOR-((t*sdm)/sqrt(k))
  ampu<-amplitude+(t*sqrt(c22))
  ampl<-amplitude-(t*sqrt(c22))
  if (ampu > 0 & ampl < 0){
    fiu<-NA
    fil<-NA
    warning("Warning: Amplitude confidence interval contains zero. Acrophase confidence interval cannot be calculated and was set to NA.")
  } else {
    fiu<-acrophase+atan(((c23*(t^2))+((t*sqrt(c33))*sqrt((amplitude^2)-(((c22*c33)-(c23^2))*((t^2)/c33)))))/((amplitude^2)-(c22*(t^2))))
    fil<-acrophase+atan(((c23*(t^2))-((t*sqrt(c33))*sqrt((amplitude^2)-(((c22*c33)-(c23^2))*((t^2)/c33)))))/((amplitude^2)-(c22*(t^2))))
  }
  conf.m<-rbind(mesoru,mesorl)
  conf.a<-rbind(ampu,ampl)
  conf.fi<-rbind(fiu,fil)
  conf.mesor<-rbind(mesoru,mesorl)
  conf.ampl<-rbind(ampu,ampl)
  conf.acro<-rbind(fiu,fil)
  conf.ints<-cbind(conf.mesor,conf.ampl,conf.acro)
  rownames(conf.ints)<-c("Upper limit", "Lower limit")
  colnames(conf.ints)<-c("MESOR", "Amplitude", "Acrophase")
  names(population.cosinor)<-c("single.cos","pop.mat","coefficients")
  data2<-subset(data, select = -time)
  emp.mean<-rowMeans(data2, na.rm=T)
  for (n in 1:nrow(data)){
    fitted.values<-c(fitted.values,MESOR+(amplitude*cos(((2*pi*time[n])/period)+acrophase)))
  }
  residuals<-emp.mean-fitted.values
  population.cosinor[[4]]<-emp.mean
  population.cosinor[[5]]<-fitted.values
  population.cosinor[[6]]<-residuals
  population.cosinor[[7]]<-conf.ints
  names(population.cosinor)<-c("single.cos","pop.mat","coefficients","emp.mean","fitted.values","residuals","conf.ints")
  class(population.cosinor)<-"population.cosinor.lm"
  print(coefficients)
  if(plot == T){
    plot(ggplot(data.frame(cbind(time, emp.mean, fitted.values)), aes(x = time))+
           geom_line(aes(y = emp.mean, linetype = "Observed"))+
           geom_line(aes(y = fitted.values, linetype = "Estimated"))+
           labs(x = "Time", y = "Value", linetype = "")+
           scale_linetype_manual("",values=c("Estimated"=2,"Observed"=1)))
  }
  invisible(population.cosinor)
  }

#' Periodogram
#'
#' @param data A data frame containing responses of subjects collected over time, with subjects in the rows and timepoints in the columns.
#' @param time A vector containing the times at which the data was collected. If this vector includes midnight, it should be coded as 24 instead of 0.
#' @param periods A vector containing periods that are to be included in the periodogram. Defaults to the same periods as provided in the vector \code{time}.
#' @param na.action Action to be performed on missing values. Defaults to \code{na.omit}.
#' @param alpha Significance level for determining if a rhythm with a given period is significant or not. Defaults to .05.
#' @description Estimates the best-fitting period using iterative cosinor procedure.
#' @details Iterative cosinor procedure is performed as described in Klemfuss & Clopton (1993). Cosinor is performed iteratively with the period (\eqn{\tau}) increased by 1 in each iteration. Percent Rhythm is calculated in each iteration, which allows for an estimation of the best fitting period. A periodogram can be plotted, which shows Percent Rhythm (coefficient of determination) for each period. On the plot, periods with significant rhythm are shown as a point and periods with insignificant rhythm are shown as a cross.
#' @note The range of periods included in iterations starts from 3 (sinusoidality of the curve is not achieved for \eqn{\tau} < 3) and ends with the number of timepoints in the data.
#' @seealso \code{\link{cosinor.PR}}
#' @examples
#' periodogram(data = PANAS_november, time = PANAS_time)
#'
#' periodogram(data = t(data.frame(temperature_zg$Temperature)), time = temperature_zg$Time)
#' @references Klemfuss, H. & Clopton, P. L. (1993). Seeking Tau: A Comparison of Six Methods. \emph{Journal of Interdisciplinary Cycle Research}, \emph{24(1)}, 1-16.
#' @import cosinor ggplot2 matrixStats Hmisc
#' @importFrom graphics par plot
#' @importFrom stats cor cov lm pf sd qt as.formula na.omit
#' @export

periodogram<-function(data, time, periods = time, na.action = na.omit, alpha = .05){
  periodogram<-matrix()
  times<-periods
  signs<-rep(NA, 2)
  unsigns<-rep(NA, 2)
  if (nrow(data) == 1) {
    data<-data.frame(cbind(t(data),time))
    colnames(data)<-c("Subj","time")
    for (i in 3:length(times)) {
      cosinor<-cosinor.lm(Subj~time(time),data=data,na.action=na.action,period=times[i])
      periodogram[[i]]<-cosinor.PR(cosinor)[[2]]
      if (cosinor.detect(cosinor)[,4] < alpha){
        signs<-c(signs, cosinor.PR(cosinor)[[2]])
        unsigns<-c(unsigns, NA)
      } else {
        signs<-c(signs, NA)
        unsigns<-c(unsigns, cosinor.PR(cosinor)[[2]])
      }
    }
  }
  else {
      for (i in 3:length(times)){
      cosinor<-population.cosinor.lm(data = data, time = time, na.action = na.action, period = times[i], plot = F)
      periodogram[[i]]<-cosinor.PR(cosinor)[[2]]
      if (cosinor.detect(cosinor)[,4] < alpha){
        signs<-c(signs, cosinor.PR(cosinor)[[2]])
        unsigns<-c(unsigns, NA)
      } else {
        signs<-c(signs, NA)
        unsigns<-c(unsigns, cosinor.PR(cosinor)[[2]])
      }
    }
  }
  periodss<-as.numeric(periodogram)
  periodogram<-data.frame(periodogram)
  rows<-(nrow(periodogram))
  plot<-ggplot(periodogram,aes(x=times))+
    geom_point(aes(y=signs)) +
    geom_point(aes(y=unsigns), shape = 4) +
    geom_line(aes(y=periodogram)) +
    labs(x = "Period", y = "Coefficient of determination")+
    scale_x_continuous(labels = as.character(times), breaks = times)
  best<-(which(periodss == max(periodss,na.rm=T)))
  besttime<-times[best]
  print(paste0("The best fitting period is ",besttime,"."))
  return(plot)
}

#' Comparison of Cosinor Parameters of Two Populations
#'
#' @param pop1 An object of the \code{population.cosinor.lm} class calculated on the first population.
#' @param pop2 An object of the \code{population.cosinor.lm} class calculated on the second population.
#' @description Runs the tests that compare MESORs, amplitudes and acrophases of two different populations.
#' @details Bingham et al. (1982) describe tests for comparing population MESORs, amplitudes and acrophases. These tests are esentially F-ratios with \eqn{df_1 = m - 1} and \eqn{df_2 = K - m}, where \eqn{m} is the number of populations and \eqn{K} is the total number of subjects. The tests for MESOR, amplitude and acrophase differences respectively are calculated as follows:
#' \deqn{F_M = \frac{\sum_{j = 1}^{m}k_j(\widehat{M}_j - \widehat{M})^2}{(m-1)\widehat{\sigma}_M^2}}
#' \deqn{F_\phi = \frac{\frac{\sum_{j = 1}^{m}k_j A_j^2 * sin^2(\widehat{\phi}_j - \tilde{\phi})}{m - 1}} {\widehat{\sigma}_\beta^2 sin^2\tilde{\phi} + 2\widehat{\sigma}_{\beta \gamma} cos\tilde{\phi}sin\tilde{\phi} + \widehat{\sigma}_\gamma^2 cos^2\tilde{\phi}}}
#' \deqn{F_A = \frac{\frac{\sum_{j = 1}^{m}(\widehat{A}_j - \widehat{A})^2}{m - 1}}{\widehat{\sigma}^2_\beta cos^2\widehat{\phi} - 2\widehat{\sigma}_{\beta \gamma}cos\widehat{\phi}sin\widehat{\phi} + \widehat{\sigma}^2_\gamma sin^2 \widehat{\phi}}}
#' where \eqn{\widehat{M}}, \eqn{\widehat{A}} and \eqn{\widehat{\phi}} are weighted averages of parameters across populations calculated as:
#' \deqn{\widehat{M} = \frac{\sum_{j = 1}^{m}k_j\widehat{M}_j}{K}}
#' \deqn{\widehat{A} = \frac{\sum_{j = 1}^{m}k_j\widehat{A}_j}{K}}
#' \deqn{\widehat{\phi} = \frac{\sum_{j = 1}^{m}k_j\widehat{\phi}_j}{K}}
#' \eqn{\tilde{\phi}} is derived from the following expression:
#' \deqn{tan 2\tilde{\phi} = \frac{\sum_{j = 1}^{m}k_j\widehat{A}^2_j sin 2\widehat{\phi}_j}{\sum_{j = 1}^{m}k_j\widehat{A}^2_j cos 2\widehat{\phi}_j}}
#' where \eqn{2\tilde{\phi}} lies between \eqn{-\frac{\pi}{2}} and \eqn{\frac{\pi}{2}} if the denomanator is positive or between \eqn{\frac{\pi}{2}} and \eqn{\frac{3\pi}{2}} if the denominator is negative, \eqn{k_j} is the number of subjects in the \eqn{j}th population, \eqn{\widehat{M}_j}, \eqn{\widehat{A}_j} and \eqn{\widehat{\phi}_j} are the cosinor parameters of the \eqn{j}th population and \eqn{\widehat{\sigma}_\beta},\eqn{\widehat{\sigma}_\gamma} and \eqn{\widehat{\sigma}_{\beta \gamma}} are the estimates of pooled standard deviations (and covariance) calculated as following:
#' \deqn{\widehat{\sigma}_\beta = \frac{\sum_{j = 1}^{m} (k_j - 1)\widehat{\sigma}_{\beta_j}}{K - m}}
#' \deqn{\widehat{\sigma}_\gamma = \frac{\sum_{j = 1}^{m} (k_j - 1)\widehat{\sigma}_{\gamma_j}}{K - m}}
#' \deqn{\widehat{\sigma}_{\beta \gamma} = \frac{\sum_{j = 1}^{m} (k_j - 1)\widehat{\sigma}_{{\beta_j} {\gamma_j}}}{K - m}}
#' where \eqn{\widehat{\sigma}_{\beta_j}}, \eqn{\widehat{\sigma}_{\gamma_j}} and \eqn{\widehat{\sigma}_{{\beta_j} {\gamma_j}}} are the standard devations and covariance of \eqn{\beta} and \eqn{\gamma} parameters in the \eqn{j}th population.
#' @note These tests should only be performed on independent samples. If the acrophases of two populations are significantly different, the results of the amplitude difference test are not reliable and should not be interpreted. While it's possible to perform tests which compare more than two populations, this function can only compare two populations.
#' @examples
#' fit.extraverts<-population.cosinor.lm(data = PA_extraverts, time = PA_time,
#' period = 24)
#' fit.introverts<-population.cosinor.lm(data = PA_introverts, time = PA_time,
#' period = 24)
#' cosinor.poptests(pop1 = fit.extraverts, pop2 = fit.introverts)
#' @references Bingham, C., Arbogast, B., Guillaume Cornélissen, G., Lee, J.K. & Halberg, F. (1982). Inferential Statistical Methods for Estimating and Comparing Cosinor Parameters. \emph{Chronobiologia}, \emph{9(4)}, 397-439.
#' @import cosinor ggplot2  matrixStats Hmisc
#' @importFrom graphics par plot
#' @importFrom stats cor cov lm pf sd qt as.formula na.omit
#' @export

cosinor.poptests<-function(pop1,pop2){
  mesors1<-list()
  mesors2<-list()
  amplitudes1<-list()
  amplitudes2<-list()
  acrophases1<-list()
  acrophases2<-list()
  betas1<-list()
  betas2<-list()
  gammas1<-list()
  gammas2<-list()
  k1<-length(pop1$single.cos)
  k2<-length(pop2$single.cos)
  K<-(k1+k2)
  for (i in 1:k1){
    mesors1[[i]]<-pop1$single.cos[[i]]$coefficients[[1]]
    amplitudes2[[i]]<-pop1$single.cos[[i]]$coefficients[[2]]
    betas1[[i]]<-pop1$single.cos[[i]]$fit$coefficients[[2]]
    gammas1[[i]]<-pop1$single.cos[[i]]$fit$coefficients[[3]]
    acrophases1[[i]]<-correct.acrophase(pop1$single.cos[[i]])
  }
  for (i in 1:k2){
    mesors2[[i]]<-pop2$single.cos[[i]]$coefficients[[1]]
    amplitudes2[[i]]<-pop2$single.cos[[i]]$coefficients[[2]]
    betas2[[i]]<-pop2$single.cos[[i]]$fit$coefficients[[2]]
    gammas2[[i]]<-pop2$single.cos[[i]]$fit$coefficients[[3]]
    acrophases2[[i]]<-correct.acrophase(pop2$single.cos[[i]])
  }
  mesors1<-matrix(mesors1)
  mesors2<-matrix(mesors2)
  betas1<-matrix(betas1)
  betas2<-matrix(betas2)
  gammas1<-matrix(gammas1)
  gammas2<-matrix(gammas2)
  amplitudes1<-matrix(amplitudes1)
  amplitudes2<-matrix(amplitudes2)
  acrophases1<-matrix(acrophases1)
  acrophases2<-matrix(acrophases2)
  M = ((k1*pop1$coefficients[[1]])+(k2*pop2$coefficients[[1]]))/K
  A = ((k1*pop1$coefficients[[2]])+(k2*pop2$coefficients[[2]]))/K
  FI = ((k1*pop1$coefficients[[3]])+(k2*pop2$coefficients[[3]]))/K
  BETA = ((k1*(pop1$coefficients[[2]]*cos(pop1$coefficients[[3]])))+(k2*(pop2$coefficients[[2]]*cos(pop2$coefficients[[3]]))))/K
  GAMMA = ((k1*(-1*pop1$coefficients[[2]]*sin(pop1$coefficients[[3]])))+(k2*(-1*pop2$coefficients[[2]]*sin(pop2$coefficients[[3]]))))/K
  TM = (k1*((pop1$coefficients[[1]]-M)^2)) + (k2*((pop2$coefficients[[1]]-M)^2))
  tann<-((k1*(pop1$coefficients[[2]]^2))*sin(2*pop1$coefficients[[3]]))+((k2*(pop2$coefficients[[2]]^2))*sin(2*pop2$coefficients[[3]]))
  tand<-((k1*(pop1$coefficients[[2]]^2))*cos(2*pop1$coefficients[[3]]))+((k2*(pop2$coefficients[[2]]^2))*cos(2*pop2$coefficients[[3]]))
  if (tand > 0) {
    twofi = atan(tann/tand)
  }
  else {
    twofi = atan(tann/tand) + pi
  }
  FITILDE = twofi/2
  varm1<-rowVars(as.numeric(mesors1), dim. = c(1,k1))
  varm2<-rowVars(as.numeric(mesors2), dim. = c(1,k2))
  varb1<-rowVars(as.numeric(betas1), dim. = c(1,k1))
  varb2<-rowVars(as.numeric(betas2), dim. = c(1,k2))
  vary1<-rowVars(as.numeric(gammas1), dim. = c(1,k1))
  vary2<-rowVars(as.numeric(gammas2), dim. = c(1,k2))
  covby1<-cov(as.numeric(betas1),as.numeric(gammas1))
  covby2<-cov(as.numeric(betas2),as.numeric(gammas2))
  varm<-(((k1-1)*varm1)/(K-2))+(((k2-1)*varm2)/(K-2))
  varb<-(((k1-1)*varb1)/(K-2))+(((k2-1)*varb2)/(K-2))
  vary<-(((k1-1)*vary1)/(K-2))+(((k2-1)*vary2)/(K-2))
  covby<-(((k1-1)*covby1)/(K-2))+(((k2-1)*covby2)/(K-2))
  FM = TM/varm
  acrn<-(pop1$coefficients[[2]]^2+((sin(pop1$coefficients[[3]]-FITILDE))^2)) + (pop2$coefficients[[2]]^2+((sin(pop2$coefficients[[3]]-FITILDE))^2))
  acrd1<-varb*((sin(FITILDE))^2)
  acrd2<-2*covby*cos(FITILDE)*sin(FITILDE)
  acrd3<-vary*((cos(FITILDE))^2)
  acrd<-acrd1-acrd2+acrd3
  FFI<-acrn/acrd
  ampn<-((pop1$coefficients[[2]]-A)^2) + ((pop2$coefficients[[2]]-A)^2)
  ampd1<-varb*((cos(FI))^2)
  ampd2<-2*covby*cos(FI)*sin(FI)
  ampd3<-vary*((sin(FI))^2)
  ampd<-ampd1-ampd2+ampd3
  FA<-ampn/ampd
  df1=1
  df2=K-2
  PM<-pf(q=FM,df1=df1,df2=df2,lower.tail=F)
  PA<-pf(q=FA,df1=df1,df2=df2,lower.tail=F)
  PFI<-pf(q=FFI,df1=df1,df2=df2,lower.tail=F)
  resm<-cbind(FM,df1,df2,PM,pop1$coefficients[[1]],pop2$coefficients[[1]])
  resa<-cbind(FA,df1,df2,PA,pop1$coefficients[[2]],pop2$coefficients[[2]])
  resfi<-cbind(FFI,df1,df2,PFI,pop1$coefficients[[3]],pop2$coefficients[[3]])
  results<-rbind(resm,resa,resfi)
  colnames(results)<-c("F","df1","df2","p", "1st population", "2nd population")
  row.names(results)<-c("MESOR", "Amplitude", "Acrophase")
  if (PFI < 0.05) {
    warning("Results of population amplitude difference test are not reliable due to different acrophases.")
  }
  return(results)
}

#' Serial Sections
#'
#' @param data A data frame containing responses of subjects collected over time, with subjects in the rows and timepoints in the columns.
#' @param time A vector containing the times at which the data was collected.
#' @param period Duration of one cycle of rhythm.
#' @param interval Length of an interval (number of timepoints) on which cosinor analyses will be ran.
#' @param increment A number indicating for how much timepoints should the interval be displaced throughout the data. Note that the value of the increment cannot be higher than the value of the interval.
#' @param na.action Action to be performed on missing values. Defaults to \code{na.omit}.
#' @param alpha Significance level for calculating population cosinor parameters confidence intervals. Defaults to .05 (confidence intervals are 5\% risk intervals).
#' @description Performs serial section analysis of rhythmic data and fits non-stationary cosinor models.
#' @details Cornélissen (2014) describes procedures for rhythmometric analysis of non-stationary data. First, an interval of an user-specified length (\eqn{I}) is chosen and usual cosinor analysis (i.e. single cosinor or population-mean cosinor) is performed on the interval. The interval is then displaced throughout the data by an user-specified increment (\eqn{\Delta t}) and cosinor analysis is then performed on the new interval. Intervals can be overlapping (\eqn{\Delta t}<I) or non-overlapping (\eqn{\Delta t}=I). A rhythm detection test is also calculated in each interval. After values of cosinor parameters have been obtained for each interval, they and their confidence intervals can be plotted, along with the \emph{p}-values from the rhythm detection test.
#' @note The value of increment cannot be higher than the value of the interval.
#' @examples
#' ssections(data = PANAS_november, time = PANAS_time, period = 7,
#' interval = 7, increment = 1)
#' @references Cornélissen, G. (2014). Cosinor-Based Rhythmometry. \emph{Theoretical Biology and Medical Modeling}, \emph{11}, Article 16.
#' @return Object of the \code{Serial sections} class containing the following objects:
#'  \item{\code{coefficients}}{Cosinor coefficients in each of the intervals.}
#'  \item{\code{emp.mean}}{Empirical mean of the data across all timepoints.}
#'  \item{\code{p-values}}{\emph{p}-values from the rhythm detection test in each interval.}
#'  \item{\code{cosinors}}{A list containing all cosinor objects calculated in each interval.}
#'  \item{\code{plots}}{Stacked plots showing the empirical data, cosinor parameters and their confidence intervals, p-values of the rhythm detection test and number of measurements across time.}
#' @import cosinor ggplot2  matrixStats Hmisc
#' @importFrom cowplot plot_grid
#' @importFrom scales pretty_breaks
#' @importFrom graphics par plot
#' @importFrom stats cor cov lm pf sd qt as.formula na.omit
#' @export

ssections<-function(data,time,period,interval,increment,na.action = na.omit,alpha = .05){
  i<-1
  j<-1
  q<-0
  cosinors<-list()
  coefficients<-data.frame(matrix(nrow=3,ncol=1))
  lmcoeffs<-data.frame(matrix(nrow=2,ncol=1))
  pvalues<-list()
  cor.acrs<-list()
  acro4<-list()
  acro5<-list()
  acro6<-list()
  interm.coefs<-list()
  sections<-list()
  fitted.values<-list()
  emp.mean<-list()
  plots<-list()
  Ns<-vector()
  Mup<-vector()
  Mdown<-vector()
  Aup<-vector()
  Adown<-vector()
  Acrup<-vector()
  Acrdown<-vector()
  while (j <= ncol(data)){
    j=i+interval-1
    if(length(i:j) == interval & j <= ncol(data)){
      q=q+1
    }
    i=i+increment
  }
  lmfitted<-data.frame(matrix(nrow=q,ncol=1))
  if (increment>interval){
    stop("Value of increment cannot be higher than the value of interval.")
  }
  if (nrow(data) == 1){
    dataa<-data.frame(cbind(t(data),time))
    colnames(dataa)<-c("Subj","time")
    for (k in c(1:q)){
      cosinors[[k]]<-cosinor.lm(Subj~time(time),period=period,na.action=na.action,data=dataa[c(k:(k+interval-1)),])
      pvalues[[k]]<-cosinor.detect(cosinors[[k]])[[4]]
      coefficients[[k]]<-cosinors[[k]]$coefficients
      cor.acrs[[k]]<-correct.acrophase(cosinors[[k]])
      Ns<-c(Ns, sum(!is.na(dataa[c(k:(k+interval-1)),])))
      Mup<-c(Mup, summary(cosinors[[k]])$transformed.table$upper.CI[1])
      Mdown<-c(Mdown, summary(cosinors[[k]])$transformed.table$lower.CI[1])
      Aup<-c(Aup, summary(cosinors[[k]])$transformed.table$upper.CI[2])
      Adown<-c(Adown, summary(cosinors[[k]])$transformed.table$lower.CI[2])
      Acrup<-c(Acrup, cor.acrs[[k]] + (summary(cosinors[[k]])$transformed.table$upper.CI[3] - summary(cosinors[[k]])$transformed.table$estimate[3]))
      Acrdown<-c(Acrdown, cor.acrs[[k]] + (summary(cosinors[[k]])$transformed.table$estimate[3] - summary(cosinors[[k]])$transformed.table$lower.CI[3]))
      k=k+increment
    }
    coefficients<-data.frame(t(coefficients))
  }
  else {
    for (k in c(1:q)){
      cosinors[[k]]<-population.cosinor.lm(data=data[,c(k:(k+interval-1))],time=time[c(k:(k+interval-1))],period=period,na.action=na.action, plot = F)
      pvalues[[k]]<-cosinor.detect(cosinors[[k]])[[4]]
      coefficients[[k]]<-t(cosinors[[k]]$coefficients)
      cor.acrs[[k]]<-cosinors[[k]]$coefficients$Acrophase
      Ns<-c(Ns, sum(!is.na(data[,c(k:(k+interval-1))])))
      Mup<-c(Mup, cosinors[[k]]$conf.ints[1,1])
      Mdown<-c(Mdown, cosinors[[k]]$conf.ints[2,1])
      Aup<-c(Aup, cosinors[[k]]$conf.ints[1,2])
      Adown<-c(Adown, cosinors[[k]]$conf.ints[2,2])
      Acrup<-c(Acrup, cosinors[[k]]$conf.ints[1,3])
      Acrdown<-c(Acrdown, cosinors[[k]]$conf.ints[2,3])
      k=k+increment
      }
  }
  Ns<-data.frame(Ns)
  acro2<-as.numeric(cor.acrs)-2*pi
  acro3<-as.numeric(cor.acrs)+2*pi
  cor.acrs<-as.numeric(cor.acrs)
  if (nrow(data) != 1){
    coefficients<-t(coefficients)
  }
  coefficients<-cbind(coefficients,acro2,cor.acrs,acro3)
  coefficients<-data.frame(coefficients)
  colnames(coefficients)<-c("MESOR","Amplitude","Acrophase","acro2","cor.acrs","acro3")
  for (l in c(1:nrow(coefficients))){
    if (coefficients[l,4] >= -2*pi & coefficients[l,4] <= 0)
    {
      coefficients[l,3]<-coefficients[l,4]
    }
    else if (coefficients[l,5] >= -2*pi & coefficients[l,5] <= 0)
    {
      coefficients[l,3]<-coefficients[l,5]
    }
    else if (coefficients[l,6] >= -2*pi & coefficients[l,6] <= 0)
    {
      coefficients[l,3]<-coefficients[l,6]
    }
  }
  cor.acrs<-data.frame(t(cor.acrs))
  acro4[[1]]<-cor.acrs[[1]]
  for (m in c(2:nrow(coefficients))){
    if(abs(as.numeric(cor.acrs[[m]])-as.numeric(acro4[[(m-1)]]) < 1*pi)){
      acro4[[m]]<-cor.acrs[[m]]
    }
    else if(abs(as.numeric(acro2[[m]])-as.numeric(acro4[[(m-1)]]) < 1*pi)){
      acro4[[m]]<-acro2[[m]]
    }
    else {
      acro4[[m]]<-acro3[[m]]
    }
  }
  acro4<-data.frame(acro4)
  acro4<-t(acro4)
  acro5<-acro4+(2*pi)
  acro6<-acro4-(2*pi)
  coefficients<-cbind(coefficients,acro4,acro5,acro6)
  colnames(coefficients)<-c("MESOR","Amplitude","Acrophase","acro2","cor.acrs","acro3","acro4","acro5","acro6")
  for (n in c(1:nrow(coefficients))){
    if (coefficients[n,7] >= -2*pi & coefficients[n,7] <= 0)
    {
      coefficients[n,3]<-coefficients[n,7]
    }
    else if (coefficients[n,8] >= -2*pi & coefficients[n,8] <= 0)
    {
      coefficients[n,3]<-coefficients[n,8]
    }
    else if (coefficients[n,9] >= -2*pi & coefficients[n,9] <= 0)
    {
      coefficients[n,3]<-coefficients[n,9]
    }
  }
  acrophases<-data.frame(coefficients[,c(7,8,9)])
  coefficients<-data.frame(coefficients[,c(-4,-5,-6,-7,-8,-9)])
  row.names(coefficients)<-c(1:nrow(coefficients))
  lmcoeffs[[1]]<-lm(coefficients$MESOR~c(1:nrow(coefficients)))$coefficients
  lmcoeffs[[2]]<-lm(coefficients$Amplitude~c(1:nrow(coefficients)))$coefficients
  lmcoeffs[[3]]<-lm(acrophases[,1]~c(1:nrow(coefficients)))$coefficients
  lmcoeffs[[4]]<-lm(acrophases[,2]~c(1:nrow(coefficients)))$coefficients
  lmcoeffs[[5]]<-lm(acrophases[,3]~c(1:nrow(coefficients)))$coefficients
  lmfitted[[1]]<-lm(coefficients$MESOR~c(1:nrow(coefficients)))$fitted.values
  lmfitted[[2]]<-lm(coefficients$Amplitude~c(1:nrow(coefficients)))$fitted.values
  lmfitted[[3]]<-lm(acrophases[,1]~c(1:nrow(coefficients)))$fitted.values
  lmfitted[[4]]<-lm(acrophases[,2]~c(1:nrow(coefficients)))$fitted.values
  lmfitted[[5]]<-lm(acrophases[,3]~c(1:nrow(coefficients)))$fitted.values
  for (o in c(1:nrow(lmfitted))){
    if (lmfitted[o,3] >= -2*pi & lmfitted[o,3] <= 0)
    {
      lmfitted[l,3]<-lmfitted[l,3]
    }
    else if (lmfitted[o,4] >= -2*pi & lmfitted[o,4] <= 0)
    {
      lmfitted[o,3]<-lmfitted[o,4]
    }
    else if (lmfitted[o,5] >= -2*pi & lmfitted[o,5] <= 0)
    {
      lmfitted[o,3]<-lmfitted[o,5]
    }
  }
  acrolms<-data.frame(lmfitted[,c(3,4,5)])
  lmfitted<-data.frame(lmfitted[,c(-4,-5)])
  colnames(lmfitted)<-c("MESOR","Amplitude","Acrophase")
  if (lmcoeffs[[3]][[1]] >= 2*pi & lmcoeffs[[3]][[1]] <=0){
    lmcoeffs[[3]]<-lmcoeffs[[3]]
  }
  else if (lmcoeffs[[4]][[1]] >= 2*pi & lmcoeffs[[4]][[1]] <=0){
    lmcoeffs[[3]]<-lmcoeffs[[4]]
  }
  else if (lmcoeffs[[5]][[1]] >= 2*pi & lmcoeffs[[5]][[1]] <=0){
    lmcoeffs[[3]]<-lmcoeffs[[5]]
  }
  lmcoeffs<-data.frame(lmcoeffs[,c(-4,-5)])
  sections[[1]]<-coefficients
  sections[[2]]<-lmcoeffs
  sections[[3]]<-lmfitted
  names(sections)<-c("coefficients","lm.coeffs","lm.fitted")
  colnames(sections$lm.coeffs)<-c("MESOR","Amplitude","Acrophase")
  row.names(sections$lm.coeffs)<-c("b","a")
  for (p in 1:ncol(data)){
    fitted.values<-c(fitted.values,((sections$lm.coeffs$MESOR[[2]]*p)+sections$lm.coeffs$MESOR[[1]])+(((sections$lm.coeffs$Amplitude[[2]]*p)+sections$lm.coeffs$Amplitude[[1]])*cos(((2*pi*time[p])/period)+((sections$lm.coeffs$Acrophase[[2]]*p)+sections$lm.coeffs$Acrophase[[1]]))))
  }
  sections[[4]]<-data.frame(t(data.frame(fitted.values)))
  names(sections)<-c("coefficients","lm.coeffs","lm.fitted","fitted.values")
  data2<-data.frame(t(data))
  emp.mean<-rowMeans(data2, na.rm = T)
  sections[[5]]<-data.frame(emp.mean)
  pvalues<-data.frame(t(data.frame(pvalues)))
  rownames(pvalues)<-NULL
  colnames(pvalues)<-"pvalues"
  sections[[6]]<-pvalues
  sections[[7]]<-cosinors
  names(sections)<-c("coefficients","lm.coeffs","lm.fitted","fitted.values","emp.mean","p-values", "cosinors")
  plots[[1]]<-ggplot(data.frame(emp.mean),aes(y = emp.mean))+
    geom_line(aes(x = time))+
    labs(x = "Time", y = "Value")+
    scale_x_continuous(breaks = time, labels = time)+
    theme(axis.title = element_text(size = 12))
  plots[[2]]<-ggplot(data.frame(coefficients),aes(x=c(1:nrow(coefficients))))+
    geom_line(aes(y=coefficients[,1]))+
    geom_ribbon(aes(ymax = Mup, ymin = Mdown), alpha = .3)+
    labs(x = "Section", y = "MESOR")+
    theme(axis.line.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size = 12))
  plots[[3]]<-ggplot(data.frame(coefficients),aes(x=c(1:nrow(coefficients))))+
    geom_line(aes(y=coefficients[,2]))+
    geom_ribbon(aes(ymax = Aup, ymin = Adown), alpha = .3)+
    labs(x = "Section", y = "Amplitude")+
    theme(axis.line.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size = 12))
  acroupper<-coefficients[,3]
  acroupper[acroupper < -pi]<-NA
  acrolower<-coefficients[,3]
  acrolower[acrolower >= -pi]<-NA
  Acrupupper<-Acrup
  Acrupupper[Acrupupper < -pi]<-NA
  Acruplower<-Acrup
  Acruplower[Acruplower >= -pi]<-NA
  Acrdownupper<-Acrdown
  Acrdownupper[Acrdownupper < -pi]<-NA
  Acrdownlower<-Acrdown
  Acrdownlower[Acrdownlower >= -pi]<-NA
  acros<-data.frame(cbind(acroupper, acrolower))
  plots[[4]]<-ggplot(acros,aes(x=c(1:nrow(acros))))+
    geom_line(aes(y=acroupper))+
    geom_line(aes(y=acrolower))+
    geom_ribbon(aes(ymax = Acrupupper, ymin = Acrdownupper), alpha = .3)+
    geom_ribbon(aes(ymax = Acruplower, ymin = Acrdownlower), alpha = .3)+
    labs(x = "Section", y = "Acrophase")+
    theme(axis.line.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size = 12))
  plots[[5]]<-ggplot(pvalues, aes(y = pvalues, x = 1:q))+
    geom_line()+
    geom_point()+
    theme(axis.line.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size = 12))+
    labs(x = "Serial section", y = "p-value")+
    geom_line(aes(y=0.05, linetype = "Cutoff"), size = 0.7)+
    geom_line(aes(y=0.01, linetype = "Cutoff"), size = 0.7)+
    theme(legend.position="none", axis.title = element_text(size = 12))+
    scale_linetype_manual("",values=c("Cutoff"=2))
  plots[[6]]<-ggplot(Ns, aes(y = Ns, x = 1:q))+
    geom_line()+
    theme(axis.title = element_text(size = 12), axis.line.x = element_line(color = "black"))+
    labs(x = "Serial section", y = "Number of \n measurements")+
    scale_x_continuous(labels = paste0(1:q, "."), breaks = 1:q)+
    scale_y_continuous(breaks = pretty_breaks())
  sections[[8]]<-plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], ncol = 1, align = "v")
  names(sections)<-c("coefficients","lm.coeffs","lm.fitted","fitted.values","emp.mean","p-values", "cosinors", "plots")
  par(ask=F)
  sections<-sections[c(-2, -3, -4)]
  class(sections)<-"Serial sections"
  plot(sections[[5]])
  invisible(sections)
}

#' Daily air temperature in Zagreb
#'
#' A dataset containing the air temperature in Zagreb on July 1, 2015.
#' @details Measurements were taken every 30 minutes throughout the whole day.
#'
#' @format A data frame with two variables:
#'  \describe{
#'   \item{Temperature}{Air temperature in Zagreb on July 1, 2015.}
#'   \item{Time}{Time of the day when the temperature was measured.}
#'   }
#' @source Croatian Meteorological and Hydrological Service, \url{http://www.meteo.hr}
"temperature_zg"

#' Self-reported mood
#'
#' A dataset containing the responses of 19 subjects on the shortened version of the PANAS questionnaire (Watson, Clark & Tellegen, 1988) in November 2015.
#' @details Measurements were taken every day after 8 PM.
#'
#' @format A data frame with 19 rows and 30 variables:
#' \describe{
#'   \item{X01, X02, X03, X04, X05, X06, X07, X08, X09, X10, X11, X12, X13, X14, X15, X16, X17, X18, X19, X20, X21, X22, X23, X24, X25, X26, X27, X28, X29, X30}{Responses of subjects at 30 measurement points (days).}}
#' @source Mutak, A. i Vukasović Hlupić, T. (2017). Exogeneity of the Circaseptan Mood Rhythm and Its Relation to the Working Week. \emph{Review of Psychology}, \emph{24} (1-2), 15-28.
#' @note The data contained in this dataset has been reduced compared to the original data that included more subjects. This dataset contains only the subjects that have responded to the PANAS questionnaire on more than 85\% of the timepoints in both of the research cycles (July and November).
#' @references Watson, D., Clark, L. A. & Tellegen, A. (1988). Development and Validation of Brief Measures of Positive and Negative Affect: The PANAS Scales. \emph{Journal of Personality and Social Psychology}, \emph{54(6)}, 1063-1070.
"PANAS_november"

#' Self-reported positive affect of extraverts
#'
#' A dataset containing the responses of 24 subjects on the Positive Affect scale of the shortened version of the PANAS questionnaire (Watson, Clark & Tellegen, 1988) in January 2017.
#' @details Measurements were taken at 10 AM, 12 PM, 2 PM, 4 PM, 6 PM and 8 PM \eqn{\pm} 30 minutes in the period of January 16 - 22, 2017. The data contained in this dataset has been averaged for each hour across 7 days of measurement.
#'
#' @format A data frame with 24 rows and 6 variables:
#' \describe{
#'   \item{X1, X2, X3, X4, X5, X6}{Responses of subjects at 6 measurement points (hours).}
#'   }
#' @source Mutak, A., Pavlović, M. & Zibar, K. (2017, May). \emph{Postoje li razlike između introverata i ekstraverata u cirkadijurnim ritmovima raspoloženja?} [\emph{Are There Differences Between Introverts and Extraverts in Circadian Mood Rhythms?}]. Study presented at the 3rd \emph{Regionalni susret studenata psihologije - SPIRI} [\emph{Regional Meeting of Psychology Students - SPIRI}] conference, Rijeka, Croatia.
#' @references Watson, D., Clark, L. A. & Tellegen, A. (1988). Development and Validation of Brief Measures of Positive and Negative Affect: The PANAS Scales. \emph{Journal of Personality and Social Psychology}, \emph{54(6)}, 1063-1070.
"PA_extraverts"

#' Self-reported positive affect of introverts
#'
#' A dataset containing the responses of 29 subjects on the Positive Affect scale of the shortened version of the PANAS questionnaire (Watson, Clark & Tellegen, 1988) in January 2017.
#' @details Measurements were taken at 10 AM, 12 PM, 2 PM, 4 PM, 6 PM and 8 PM \eqn{\pm} 30 minutes in the period of January 16 - 22, 2017. The data contained in this dataset has been averaged for each hour across 7 days of measurement.
#'
#' @format A data frame with 29 rows and 6 variables:
#' \describe{
#'   \item{X1, X2, X3, X4, X5, X6}{Responses of subjects at 6 measurement points (hours).}
#'   }
#' @source Mutak, A., Pavlović, M. & Zibar, K. (2017, May). \emph{Postoje li razlike između introverata i ekstraverata u cirkadijurnim ritmovima raspoloženja?} [\emph{Are There Differences Between Introverts and Extraverts in Circadian Mood Rhythms?}]. Study presented at the 3rd \emph{Regionalni susret studenata psihologije - SPIRI} [\emph{Regional Meeting of Psychology Students - SPIRI}] conference, Rijeka, Croatia.
#' @references Watson, D., Clark, L. A. & Tellegen, A. (1988). Development and Validation of Brief Measures of Positive and Negative Affect: The PANAS Scales. \emph{Journal of Personality and Social Psychology}, \emph{54(6)}, 1063-1070.
"PA_introverts"

#' Measurement times of self-reported mood
#'
#' A dataset containing the measurement times (dates) of self reported mood contained in the data frame \code{PANAS_november}
#'
#' @format A numeric vector of length 30.
#' @source Mutak, A. i Vukasović Hlupić, T. (2017). Exogeneity of the Circaseptan Mood Rhythm and Its Relation to the Working Week. \emph{Review of Psychology}, \emph{24} (1-2), 15-28.
"PANAS_time"

#' Measurement times of self-reported positive affect
#'
#' A dataset containing the measurement times (hours) of self reported positive affect contained in the data frames \code{PA_extravers} and \code{PA_introverts}.
#'
#' @format A numeric vector of length 6.
#' @source Mutak, A., Pavlović, M. & Zibar, K. (2017, May). \emph{Postoje li razlike između introverata i ekstraverata u cirkadijurnim ritmovima raspoloženja?} [\emph{Are There Differences Between Introverts and Extraverts in Circadian Mood Rhythms?}]. Study presented at the 3rd \emph{Regionalni susret studenata psihologije - SPIRI} [\emph{Regional Meeting of Psychology Students - SPIRI}] conference, Rijeka, Croatia.
"PA_time"

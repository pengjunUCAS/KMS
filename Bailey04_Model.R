############
###################################################################################
### R program for simulating the comprehensive kinetic model by Bailey (2004).
###
### Reference: Bailey RM, 2004. Paper I-simulation of dose absorption in quartz
###            over geological timescales and its implications for the precision
###            and accuracy of optical dating. Radiation Measurements, 38: 299-310.
###
### Dependence: The R software; R package deSolve.
###
### Authors: Original Mathematica code by Vasilis Pagonis,
###          transform to R by Jun Peng.
####################################################################################
    Bailey04 <- function(t, inis, pars) {
        with(as.list(c(inis, pars)),{
            ###
            dn1 <-   A1*(N1-n1)*nc - n1*pValue*th1*exp(-E1th/k/(273+temp+hRate*t)) - n1*s1*exp(-E1/k/(273+temp+hRate*t))
            dn2 <-   A2*(N2-n2)*nc - n2*pValue*th2*exp(-E2th/k/(273+temp+hRate*t)) - n2*s2*exp(-E2/k/(273+temp+hRate*t))
            dn3 <-   A3*(N3-n3)*nc - n3*pValue*th3*exp(-E3th/k/(273+temp+hRate*t)) - n3*s3*exp(-E3/k/(273+temp+hRate*t))
            dn4 <-   A4*(N4-n4)*nc - n4*pValue*th4*exp(-E4th/k/(273+temp+hRate*t)) - n4*s4*exp(-E4/k/(273+temp+hRate*t))
            dn5 <-   A5*(N5-n5)*nc - n5*pValue*th5*exp(-E5th/k/(273+temp+hRate*t)) - n5*s5*exp(-E5/k/(273+temp+hRate*t))
            dn6 <-   A6*(N6-n6)*nc - n6*pValue*th6*exp(-E6th/k/(273+temp+hRate*t)) - n6*s6*exp(-E6/k/(273+temp+hRate*t))
            dn7 <-   A7*(N7-n7)*nc - n7*pValue*th7*exp(-E7th/k/(273+temp+hRate*t)) - n7*s7*exp(-E7/k/(273+temp+hRate*t))
            dn8 <-   A8*(N8-n8)*nc - n8*pValue*th8*exp(-E8th/k/(273+temp+hRate*t)) - n8*s8*exp(-E8/k/(273+temp+hRate*t))
            ###
            dn9 <- A9*(N9-n9)*nv - n9*s9*exp(-E9/k/(273+temp+hRate*t)) - nc*n9*B9
            dn10 <- A10*(N10-n10)*nv - n10*s10*exp(-E10/k/(273+temp+hRate*t)) - nc*n10*B10
            dn11 <- A11*(N11-n11)*nv - n11*s11*exp(-E11/k/(273+temp+hRate*t)) - nc*n11*B11
            dn12 <- A12*(N12-n12)*nv - n12*s12*exp(-E12/k/(273+temp+hRate*t)) - nc*n12*B12
            ###
            dnc <- rValue - (dn1+dn2+dn3+dn4+dn5+dn6+dn7+dn8) - nc*(n9*B9+n10*B10+n11*B11+n12*B12)
            dnv <- rValue - (dn9+dn10+dn11+dn12) - nc*(n9*B9+n10*B10+n11*B11+n12*B12)
            ###
            list(c(dn1, dn2, dn3, dn4, dn5, dn6, dn7, 
                   dn8, dn9, dn10, dn11, dn12, dnc, dnv))
        })
    } # end function Pagonis.
    ###
    ### 
    setPars <- function () {
            pars <<- c(N1=1.42e10,   E1=0.97,  s1=5e12,  A1=1e-8,   th1=1e-19,  E1th=0.1,
                       N2=1.5e9,     E2=1.55,  s2=5e14,  A2=1e-8,   th2=0,      E2th=0, 
                       N3=2.05e11,   E3=1.7,   s3=5e12,  A3=1e-9,   th3=1e-16,  E3th=0.1,
                       N4=7.04e10,   E4=1.72,  s4=5e13,  A4=8e-10,  th4=3e-17,  E4th=0.13,
                       N5=1.77e11,   E5=1.8,   s5=5e13,  A5=8e-10,  th5=4e-18,  E5th=0.2,
                       N6=2.53e11,   E6=1.65,  s6=5e13,  A6=5e-10,  th6=3e-19,  E6th=0.2,
                       N7=3.58e12,   E7=2.6,   s7=5e13,  A7=2e-10,  th7=2e-21,  E7th=0.2,
                       N8=1.28e13,   E8=2.0,   s8=1e10,  A8=1e-10,  th8=0,      E8th=0,
                       N9=8.83e13,   E9=5.0,   s9=1e13,  A9=1e-9,   B9=1e-10, 
                       N10=1.15e14,  E10=5.0,  s10=1e13, A10=1e-10, B10=1e-10,
                       N11=4.16e12,  E11=1.75, s11=5e14, A11=1e-9,  B11=5e-10,
                       N12=4.2e11,   E12=1.43, s12=5e13, A12=5e-8,  B12=5e-9,
                       k=8.617e-5, rValue=0, pValue=0, hRate=0, temp=0)
    } # end function setPars.
    ###
    ###
    setInis <- function() {
        inis <<- c(n1=0, n2=0, n3=0, n4=0, n5=0, 
                   n6=0, n7=0, n8=0, n9=0, n10=0, 
                   n11=0, n12=0, nc=0, nv=0)
    } # end function setInis.
    ###
    ###
    irradiate<- function(temp, tim, doseRate)  { 
        stopifnot(temp>0, tim>0, doseRate>0)
        pars["rValue"] <<- 2.5e10*doseRate
        pars["pValue"] <<- 0
        pars["hRate"] <<- 0
        pars["temp"] <<- temp
        times <- seq(from=0, to=tim, by=tim/100)
        res <- ode(y=inis, times=times,
                   func=Bailey04, parms=pars,
                   rtol=1e-3, atol=1e-3, hmax=1e13,
                   method="bdf", maxsteps=1e5)
        inis <<- res[length(times),-1L,drop=TRUE]
        invisible(res[-1L,,drop=FALSE])
    } # end function irradiate.
    ###
    ###
    heatAt <- function(temp, tim) {
        stopifnot(temp>0, tim>0)
        pars["rValue"] <<- 0
        pars["pValue"] <<- 0
        pars["hRate"] <<- 0
        pars["temp"] <<- temp
        times <- seq(from=0, to=tim, by=tim/100)
        res <- ode(y=inis, times=times,
                   func=Bailey04, parms=pars,
                   rtol=1e-3, atol=1e-3, hmax=1e13,
                   method="bdf", maxsteps=1e5)
        inis <<- res[length(times),-1L,drop=TRUE]
        invisible(res[-1L,,drop=FALSE])
    }  # end function heatAt.
    ###
    ###
    heatTo <- function(temp1, temp2, hRate) {
        stopifnot(temp1>0, temp2>0, 
                  sign(temp2-temp1)==sign(hRate))
        pars["rValue"] <<- 0
        pars["pValue"] <<- 0
        pars["hRate"] <<- hRate
        pars["temp"] <<- temp1
        times <- seq(from=0, 
                     to=(temp2-temp1)/hRate, 
                     by=(temp2-temp1)/hRate/100)
        res <- ode(y=inis, times=times,
                   func=Bailey04, parms=pars,
                   rtol=1e-3, atol=1e-3, hmax=1e13,
                   method="bdf", maxsteps=1e5)
        inis <<- res[length(times),-1L,drop=TRUE]
        invisible(res[-1L,,drop=FALSE])
    }  # end function heatTo.
    ###
    ###
    stimOSL<- function(temp, tim, pValue, nChannel) {
        stopifnot(temp>0, tim>0, pValue>0, 
                  nChannel>=10)
        pars["rValue"] <<- 0
        pars["hRate"] <<- 0
        pars["pValue"] <<- pValue
        pars["temp"] <<- temp
        times <- seq(from=0, to=tim, by=tim/nChannel)
        res <-  ode(y=inis, times=times,
                    func=Bailey04, parms=pars,
                    rtol=1e-3, atol=1e-3, hmax=1e13,
                    method="bdf", maxsteps=1e5)
        inis <<- res[length(times),-1L,drop=TRUE]
        oslx <- times
        osly <- (B9<-pars["B9"])*
                (n9<-res[,"n9",drop=TRUE])*
                (nc<-res[,"nc",drop=TRUE])/
                (1.0+2.8*1e7*exp(-0.64/pars["k"]/
                (273.0+pars["temp"])))
        outPut <- (cbind(res,oslx,osly))[-1L,,drop=FALSE]
        invisible(outPut)
    }  # end function stimOSL.
    ###
    ###
    stimTL <- function(lowTemp, upTemp, hRate, nChannel) {
        stopifnot(lowTemp>0, hRate>0,
                  upTemp>lowTemp, 
                  nChannel>=10) 
        pars["rValue"] <<- 0
        pars["pValue"] <<- 0
        pars["hRate"] <<- hRate
        pars["temp"] <<- lowTemp
        times <- seq(from=0, 
                     to=(upTemp-lowTemp)/hRate, 
                     by=(upTemp-lowTemp)/hRate/(nChannel-1))
        res <- ode(y=inis, times=times,
                   func=Bailey04, parms=pars,
                   rtol=1e-3, atol=1e-3, hmax=1e13,
                   method="bdf", maxsteps=1e5)
        inis <<- res[length(times),-1L,drop=TRUE]
        tlx <- pars["temp"] + pars["hRate"]*times
        tly <- (B9<-pars["B9"])*
               (n9<-res[,"n9",drop=TRUE])*
               (nc<-res[,"nc",drop=TRUE])/
               (1.0+2.8*1e7*exp(-0.64/pars["k"]/
               (273.0+pars["temp"] + pars["hRate"]*times)))
        outPut <- (cbind(res,tlx,tly))
        invisible(outPut)
    } # end function stimTL.
##############################
######################
#############

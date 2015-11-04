############
########################################################################################
### R program for simulating the comprehensive kinetic model by Pagonis et al. (2008).
###
### Reference: Pagonis V, Wintle AG, Chen R, Wang XL, 2008. A theoretical model 
###            for a new dating protocol for quartz based on thermally transferred 
###            OSL (TT-OSL). Radiation Measurements, 43: 704-708.
###
### Dependence: The R software; R package deSolve.
###
### Authors: Original Mathematica code by Vasilis Pagonis,
###          transform to R by Jun Peng.
########################################################################################
    Pagonis <- function(t, inis, pars) {
        with(as.list(c(inis, pars)),{
            ###
            dn1 <-   A1*(N1-n1)*nc - n1*pValue*th1*exp(-E1th/k/(273+temp+hRate*t)) -  n1*s1*exp(-E1/k/(273+temp+hRate*t))
            dn2 <-   A2*(N2-n2)*nc - n2*pValue*th2*exp(-E2th/k/(273+temp+hRate*t)) -  n2*s2*exp(-E2/k/(273+temp+hRate*t))
            dn3 <-   A3*(N3-n3)*nc - n3*pValue*th3*exp(-E3th/k/(273+temp+hRate*t)) -  n3*s3*exp(-E3/k/(273+temp+hRate*t))
            dn4 <-   A4*(N4-n4)*nc - n4*pValue*th4*exp(-E4th/k/(273+temp+hRate*t)) -  n4*s4*exp(-E4/k/(273+temp+hRate*t))
            dn5 <-   A5*(N5-n5)*nc - n5*pValue*th5*exp(-E5th/k/(273+temp+hRate*t)) -  n5*s5*exp(-E5/k/(273+temp+hRate*t))
            dn10 <-  A10*(N10-n10)*nc - n10*pValue*th10*exp(-E10th/k/(273+temp+hRate*t)) - n10*s10*exp(-E10/k/(273+temp+hRate*t))
            dn11 <-  A11*(N11-n11)*nc - n11*pValue*th11*exp(-E11th/k/(273+temp+hRate*t)) - n11*s11*exp(-E11/k/(273+temp+hRate*t)) 
            ###
            dn6 <- A6*(N6-n6)*nv - n6*s6*exp(-E6/k/(273+temp+hRate*t)) - nc*n6*B6
            dn7 <- A7*(N7-n7)*nv - n7*s7*exp(-E7/k/(273+temp+hRate*t)) - nc*n7*B7
            dn8 <- A8*(N8-n8)*nv - n8*s8*exp(-E8/k/(273+temp+hRate*t)) - nc*n8*B8
            dn9 <- A9*(N9-n9)*nv - n9*s9*exp(-E9/k/(273+temp+hRate*t)) - nc*n9*B9
            ###
            dnc <- rValue - (dn1+dn2+dn3+dn4+dn5+dn10+dn11) - nc*(n6*B6+n7*B7+n8*B8+n9*B9)
            dnv <- rValue - (dn6+dn7+dn8+dn9) - nc*(n6*B6+n7*B7+n8*B8+n9*B9)
            ###
            list(c(dn1, dn2, dn3, dn4, dn5, dn6, dn7, 
                   dn8, dn9, dn10, dn11, dnc, dnv))
        })
    } # end function Pagonis.
    ###
    ### 
    setPars <- function () {
            pars <<- c(N1=1.5e7,   E1=0.97,   s1=5e12,     A1=1e-8,    th1=0.75,   E1th=0.1,
                       N2=1e7,     E2=1.55,   s2=5e14,     A2=1e-8,    th2=0,      E2th=0, 
                       N3=4e7,     E3=1.73,   s3=6.5e13,   A3=5e-9,    th3=6.0,    E3th=0.1,
                       N4=2.5e8,   E4=1.8,    s4=1.5e13,   A4=5e-10,   th4=4.5,    E4th=0.13,
                       N5=5e10,    E5=2.0,    s5=1e10,     A5=1e-10,   th5=0,      E5th=0,
                       N10=5e9,    E10=1.65,  s10=6.5e13,  A10=1e-11,  th10=0.01,  E10th=0.2,
                       N11=4e9,    E11=1.6,   s11=5e12,    A11=6e-12,  th11=0,     E11th=0,
                       N6=3e8,     E6=1.43,   s6=5e13,     A6=5e-7,    B6=5e-9, 
                       N7=1e10,    E7=1.75,   s7=5e14,     A7=1e-9,    B7=5e-10,
                       N8=3e10,    E8=5.0,    s8=1e13,     A8=1e-10,   B8=1e-10,
                       N9=1.2e12,  E9=5.0,    s9=1e13,     A9=1e-14,   B9=3e-10,
                       k=8.617e-5, rValue=0,  pValue=0, hRate=0, temp=0)
    } # end function setPars.
    ###
    ###
    setInis <- function() {
        inis <<- c(n1=0, n2=0, n3=0, n4=0, n5=0, 
                   n6=0, n7=0, n8=0, n9=0, n10=0, 
                   n11=0, nc=0, nv=0)
    } # end function setInis.
    ###
    ###
    irradiate<- function(temp, tim, doseRate)  { 
        stopifnot(temp>0, tim>0, doseRate>0)
        pars["rValue"] <<- 5e7*doseRate
        pars["pValue"] <<- 0
        pars["hRate"] <<- 0
        pars["temp"] <<- temp
        times <- seq(from=0, to=tim, by=tim/100)
        res <- ode(y=inis, times=times,
                   func=Pagonis, parms=pars,
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
                   func=Pagonis, parms=pars,
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
                   func=Pagonis, parms=pars,
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
        res <- ode(y=inis, times=times,
                   func=Pagonis, parms=pars,
                   rtol=1e-3, atol=1e-3, hmax=1e13,
                   method="bdf", maxsteps=1e5)
        inis <<- res[length(times),-1L,drop=TRUE]
        oslx <- times
        osly <- (B8<-pars["B8"])*
                (n8<-res[,"n8",drop=TRUE])*
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
                   func=Pagonis, parms=pars,
                   rtol=1e-3, atol=1e-3, hmax=1e13,
                   method="bdf", maxsteps=1e5)
        inis <<- res[length(times),-1L,drop=TRUE]
        tlx <- pars["temp"] + pars["hRate"]*times
        tly <- (B8<-pars["B8"])*
               (n8<-res[,"n8",drop=TRUE])*
               (nc<-res[,"nc",drop=TRUE])/
               (1.0+2.8*1e7*exp(-0.64/pars["k"]/
               (273.0+pars["temp"] + pars["hRate"]*times)))
        outPut <- (cbind(res,tlx,tly))
        invisible(outPut)
    } # end function stimTL.
##############################
######################
#############

    
    rm(list=ls())
    library(deSolve)
    library(numOSL)
    source("Bailey01_Model.R")

    ### A function for simulating the SAR-OSL protocol 
    ### according to the kinetic model by Bailey (2001).
    sar <- function() {
        ### Initialize the model.
        setInis()

        ### Burial dose and test dose.
        nDose <- 200
        testDose <- 0.1*nDose

        ### A series of regenerative doses. 
        reDose <- c(1e-13, 0.4*nDose, 0.8*nDose, 1.2*nDose, 
                    1.6*nDose, 1e-13, 0.4*nDose)

        # Natural dose rate (1e-11 Gy/s).
        ndr <- 1e-11   

        # Laboratory dose rate (1 Gy/s).
        ldr <- 1     

        ### Normalized photon fulx.
        photonFlux <- 2.0

        ### 1000 Gy dose (the geological dose) using 1 Gy/s at 20 deg.C.
        irradiate(temp=20, tim=1000, doseRate=1)
        ### Relaxation period. 
        heatAt(temp=20, tim=60) 
       
        ### Heat to 350 deg.C (geological time).
        heatTo(temp1=20, temp2=350, hRate=5)

        ### Illumination at 200 deg.C for 100 s.
        heatTo(temp1=350, temp2=200, hRate=-5) 
        stimOSL(temp=200, tim=100, pValue=photonFlux, nChannel=1000)

        ### Given burial dose using natural dose rate 0.01 Gy/s at 20 deg.C. 
        heatTo(temp1=200, temp2=20, hRate=-5)
        irradiate(temp=20, tim=nDose/ndr, doseRate=ndr)
        ### Relaxation period.
        heatAt(temp=20, tim=60)
        ###
        len <- length(reDose)
        LxTx<-sLxTx<-vector(length=len)     
         
        ### Simulate the regenerative-dose and test-dose OSL responses.
        ### Calculate the sensitivity-corrected regenerative OSL signals (Lx/Tx).
        for (i in seq(len)) { 
            ### Administer regnerative dose.
            irradiate(temp=20, tim=reDose[i]/ldr, doseRate=ldr) 
            heatAt(temp=20, tim=60)
            heatTo(temp1=20, temp2=260, hRate=5) 
            ### Preheat temperature, 260 deg. C.
            heatAt(temp=260, tim=10)
            heatTo(temp1=260, temp2=125, hRate=-5) 
            Lxdat<-stimOSL(temp=125, tim=100, pValue=photonFlux, nChannel=1000)
            heatTo(temp1=125, temp2=20, hRate=-5)
            ###
            ###
            ### Administer test dose.
            irradiate(temp=20, tim=testDose/ldr, doseRate=ldr)
            heatAt(temp=20, tim=60) 
            heatTo(temp1=20, temp2=220, hRate=5) 
            ### Cut heat temperature, 220 deg. C.
            heatAt(temp=220, tim=10) 
            heatTo(temp1=220, temp2=125, hRate=-5) 
            Txdat<-stimOSL(temp=125, tim=100, pValue=photonFlux, nChannel=1000)
            heatTo(temp1=125, temp2=20, hRate=-5) 
            ###
            ### Integration interval for OSL signal: first five channels, 
            ### and last 100 channels for background subtraction.
            Lx <- (sum(Lxdat[1:5,"osly"]) - 5*mean(Lxdat[901:1000, "osly"]))
            Tx <- (sum(Txdat[1:5,"osly"]) - 5*mean(Txdat[901:1000, "osly"]))
            LxTx[i]<- Lx/Tx   
            ### Constant relative standard error for LxTx.               
            sLxTx<-0.03*LxTx 
        } # end for.
        ###
        ###
        ### Fit the growth curve and calculate 
        ### the equivalent dose value.
        Curvedata<-cbind(reDose[-1], LxTx[-1], sLxTx[-1])
        Ltx<-c(LxTx[1], sLxTx[1])
        res <- calED(Curvedata, Ltx, model="exp")
        return(res$ED)
    } # end function sar.
    #########################
    ##################


    #pdf("Fig1.pdf")
    ### Number of versions of random 
    ### kinetic parameters.
    nsim <- 100    
  
    ### Percent within the original
    ### kinetic parameters.
    perct <- 0.2

    EDmat <- matrix(nrow=nsim, ncol=2)
     
    ###
    iter <- 0
    setPars()
    system.time( 
         ### Random simulation.
         repeat  {  
             ### Change initial concentrations of electron 
             ### and hole randomly using a uniform distribution. 
             pars["N1"] <-  runif(1, 1.5e7*(1.0-perct), 1.5e7*(1.0+perct))
             pars["N2"] <-  runif(1, 1e7*(1.0-perct),   1e7*(1.0+perct))
             pars["N3"] <-  runif(1, 1e9*(1.0-perct),   1e9*(1.0+perct))
             pars["N4"] <-  runif(1, 2.5e8*(1.0-perct), 2.5e8*(1.0+perct))
             pars["N5"] <-  runif(1, 5e10*(1.0-perct),  5e10*(1.0+perct))
             pars["N6"] <-  runif(1, 3e8*(1.0-perct),   3e8*(1.0+perct))
             pars["N7"] <-  runif(1, 1e10*(1.0-perct),  1e10*(1.0+perct))
             pars["N8"] <-  runif(1, 1e11*(1.0-perct),  1e11*(1.0+perct))
             pars["N9"] <-  runif(1, 5e9*(1.0-perct),   5e9*(1.0+perct))
             ###
             #print(pars[c("N1", "N2", "N3", "N4", "N5", 
                          #"N6", "N7", "N8", "N9")])
             ###
             sarED <- try(sar(), silent=TRUE)
             if (class(sarED)!="try-error")  {
                 iter <- iter + 1
                 EDmat[iter,] <- sarED
                 ###
                 print(iter)
                 print(sarED)
             } # end if.
             ###
             if (iter>=nsim) break       
         } # end repeat.
    ) # 
    #dev.off()

    #pdf("Fig2.pdf")
    dbED(EDmat)
    #dev.off()


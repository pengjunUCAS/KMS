
    rm(list=ls())
    library(deSolve)
    library(numOSL)
    source("Bailey04_Model.R")

    ### A function for simulating the SAR-OSL protocol 
    ### according to the kinetic model by Bailey (2004).
    sar <- function() {
        ### Initialize the model.
        setInis()
      
        ### Burial dose and test dose.
        nDose <- 200 
        testDose <- 0.1*nDose

        ### A series of regenerative doses. 
        reDose <- c(1e-13, 0.4*nDose, 0.8*nDose, 1.2*nDose, 
                    1.6*nDose, 1e-13, 0.4*nDose)

        # Natural dose rate (1 Gy/ka).
        ndr <- 3.17098e-11  

        # Laboratory dose rate (0.1 Gy/s).
        ldr <- 0.1                     

        ### Calculate absolute photon flux.    
        waveLength <- 470  #(nm)
        stimPower <- 20    #(mW /cm^2)
        hv<- (h<-6.63e-34)*(v<-3e8/(waveLength*1e-9))
        photonFlux <- stimPower/1e3/hv

        ### 50 Ma at 20 deg.C with a dose rate of 1 Gy/ka.
        irradiate(temp=20, tim=50000/ndr, doseRate=ndr) 
        ### Relaxation period. 
        heatAt(temp=20, tim=60) 
       
         ### 6000 s illumination at 20 deg.C.
        stimOSL(temp=20, tim=6000, pValue=photonFlux, nChannel=6000) 

        ### Cycles of 10 ka irradiation with a dose rate of 1 Gy/ka 
        ### + 6000 s illumination at 20 deg.C.
        nCycle <- 2
        for (i in seq(nCycle)) { 
            irradiate(temp=20,tim=10/ndr, doseRate=ndr) 
            heatAt(temp=20, tim=60)
            stimOSL(temp=20, tim=6000, pValue=photonFlux, nChannel=6000) 
        } #end for.

        ### Burial dose - at 20 deg.C at
        ### a very low natural dose rate.
        irradiate(temp=20,tim=nDose/ndr, doseRate=ndr)
        heatAt(temp=20, tim=60)

        len <- length(reDose)
        LxTx<-sLxTx<-vector(length=len) 
       
        ###
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
             ### Change initial concentrations of electron and 
             ### hole randonly using a uniform distribution. 
             pars["N1"] <-  runif(1, 1.42e10*(1.0-perct), 1.42e10*(1.0+perct))
             pars["N2"] <-  runif(1, 1.5e9*(1.0-perct),   1.5e9*(1.0+perct))
             pars["N3"] <-  runif(1, 2.05e11*(1.0-perct), 2.05e11*(1.0+perct))
             pars["N4"] <-  runif(1, 7.04e10*(1.0-perct), 7.04e10*(1.0+perct))
             pars["N5"] <-  runif(1, 1.77e11*(1.0-perct), 1.77e11*(1.0+perct))
             pars["N6"] <-  runif(1, 2.53e11*(1.0-perct), 2.53e11*(1.0+perct))
             pars["N7"] <-  runif(1, 3.58e12*(1.0-perct), 3.58e12*(1.0+perct))
             pars["N8"] <-  runif(1, 1.28e13*(1.0-perct), 1.28e13*(1.0+perct))
             pars["N9"] <-  runif(1, 8.83e13*(1.0-perct), 8.83e13*(1.0+perct))
             pars["N10"] <- runif(1, 1.15e14*(1.0-perct), 1.15e14*(1.0+perct))
             pars["N11"] <- runif(1, 4.16e12*(1.0-perct), 4.16e12*(1.0+perct))
             pars["N12"] <- runif(1, 4.2e11*(1.0-perct),  4.2e11*(1.0+perct))
             ###
             #print(pars[c("N1", "N2", "N3", "N4", "N5", "N6", 
                          #"N7", "N8", "N9", "N10", "N11", "N12")])
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


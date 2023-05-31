##' @title LINKAGES main function for spinup with Monte Carlo Simulations
##' @author Ann Raiho
##' edited Hayden Gallo May 2023
##'
##' @param iplot        PEcAn settings object
##' @param nyear        number of years
##' @param nspec        number of species
##' @param fc           field capacity
##' @param dry          wilting point
##' @param bgs          beginning of growing season
##' @param egs          end of growing season
##' @param max.ind      maximum number of indiviuals allowed to grow in a plot
##' @param plat         latitude
##' @param temp.mat     monthly temperature matrix (nyear x 12)
##' @param precip.mat   monthly precipitation matrix (nyear x 12)
##' @param spp.params   species parameter matrix
##' @param fdat         underground parameter matrix
##' @param clat         climate correction factor matrix
##' @param basesc       initial humus weight
##' @param basesn       initial N in soil
##'
##' @description Main function for running all LINKAGES subroutines
##'
##' @return year vector of years
##' @return ag.biomass above ground biomass in kgC/m2
##' @return total.soil.carbon soil organic matter + leaf litter in kgC/m2
##' @return leaf.litter leaf litter in kgC/m2
##' @return ag.npp above ground net primary production in kgC/m2/second
##' @return hetero.resp heterotrophic respiration in kgC/m2/second
##' @return nee net ecosystem exchange in kgC/m2/second
##' @return et annual evapotransipiration kgC/m2/s
##' @return agb.pft = agb.pft above ground biomass by plant functional type (species in LINKAGES)
##' @return f.comp = f.comp fractional composition (sums to 1 each year)
##' @return ntrees.birth=ntrees.birth number of trees born of each species each year
##' @return ntrees.kill = ntrees.kill number of trees killed of each species each year
##' @return tstem=tstem stem density
##' @return tab=tab total aboveground biomass without unit conversion
##' @return fl=fl total leaf litter without unit conversion
##' @return totl=totl leaf litter nitrogen
##' @return tnap=tnap total net above ground production without unit conversion
##' @return avln=avln available nitrogen
##' @return cn=cn carbon to nitrogen ratio
##' @return sco2c=sco2c soil repiration without unit conversion
##' @return som=som soil organic matter
##' @return bar=bar biomass by species without unit conversion
##' @return aet.save=aet.save annual evapotranspiration without unit conversion
##' @return nogro.save=nogro.save matrix of trees not growing each year
##' @return dbh.save=dbh.save matrix of dbh increment of trees each year
##' @return iage.save=iage.save matrix of age of each tree each year
##'

linkages <- function(linkages.input, outdir, spinup_input){

  if(is.null(spinup_input)){

      load(linkages.input)

      max.ind <- 200

      #temp.mat <- matrix(temp.mat,nyear,12)
      #precip.mat <- matrix(precip.mat,nyear,12)

      #Storage
      gf.vec.save <- array(NA,dim=c(nspec,4,nyear,iplot))
      tstem = matrix(0,nyear,iplot) #number of stems
      area = matrix(0,nyear,iplot)
      water = matrix(0,nyear,iplot)
      tab = matrix(0,nyear,iplot) #total aboveground biomass
      abvgrnwood = matrix(0,nyear,iplot) #total aboveground woody biomass
      fl = matrix(0,nyear,iplot) #leaf litter
      totl = matrix(0,nyear,iplot) #leaf litter N
      tnap = matrix(0,nyear,iplot) #net aboveground production
      avln = matrix(0,nyear,iplot) #available nitrogen
      cn = matrix(0,nyear,iplot) #humus C:N ratio
      sco2c = matrix(0,nyear,iplot) #soil co2 evolution
      som = matrix(0,nyear,iplot) #soil organic matter
      aet.save = matrix(0,nyear,iplot)
      ncohrt.save = matrix(0,nyear,iplot)
      tyl.save = array(0,dim=c(20,nyear,iplot))
      ntrees.birth <- array(0,dim=c(nspec,nyear,iplot))
      ntrees.grow <- array(0,dim=c(nspec,nyear,iplot))
      ntrees.kill <- array(0,dim=c(nspec,nyear,iplot))
      # gf.vec.save <- array(0,dim=c(4,nyear,iplot))
      bar <- array(0,dim=c(nspec,nyear,iplot))
      algf.save.keep<- array(NA,dim=c(max.ind,nspec,nyear,iplot))
      nogro.save <- array(0,dim=c(max.ind,nyear,iplot))
      dbh.save <- array(0,dim=c(max.ind,nyear,iplot))
      iage.save <- array(0,dim=c(max.ind,nyear,iplot))
      awp.save <- array(0,dim=c(max.ind,nyear,iplot))
      npp.spp.save <- array(0,dim=c(nspec,nyear,iplot))

      for(k in 1:iplot){ #loop over plots

          plotin.out <- plotin(iplot = k, basesc = basesc, basesn = basesn, max.ind = max.ind,
                               nspec = nspec) # initializes storage matrices with zeros for each plot

          ncohrt <- unlist(plotin.out$ncohrt, use.names = FALSE)
          tyl <- unlist(plotin.out$tyl, use.names = FALSE)
          C.mat <- unlist(plotin.out$C.mat, use.names = FALSE)
          ntrees <- unlist(plotin.out$ntrees, use.names = FALSE)
          dbh <- unlist(plotin.out$dbh, use.names = FALSE)
          nogro <- unlist(plotin.out$nogro, use.names = FALSE)
          ksprt <- unlist(plotin.out$ksprt, use.names = FALSE)
          iage <- unlist(plotin.out$iage, use.names = FALSE)


        Rprof(interval=.005)
        for(i in 1:nyear){

          #calculates degree days for the year
          degd <- tempe(temp.vec = temp.mat[i,1:12])
          #degd = unlist(tempe.out$degd, use.names = FALSE)

          bgs <- (which(temp.mat[i,]>10)[1]-1)*30
          egs <- bgs + length(which(temp.mat[i,]>10))*30

          #calculates aet
          moist.out <- moist(kyr = i, temp.vec = temp.mat[i,1:12], precip.vec = precip.mat[i,1:12],
                             fc = fc, dry = dry, bgs = bgs, egs = egs, plat = plat, clat = clat)
          aet <- unlist(moist.out$aet, use.names = FALSE)
          fj <- unlist(moist.out$fj, use.names = FALSE)

          #decomposition subroutine
          decomp.out <- decomp(fdat = fdat, aet = aet,
                               ncohrt = ncohrt, fc = fc, dry = dry,
                               tyl = tyl, C.mat = C.mat)

          ff <- unlist(decomp.out$ff, use.names = FALSE)
          availn <- unlist(decomp.out$availn, use.names = FALSE)
          tyln <- unlist(decomp.out$tyln, use.names = FALSE)
          hcn <- unlist(decomp.out$hcn, use.names = FALSE)
          sco2 <- unlist(decomp.out$sco2, use.names = FALSE)
          ncohrt <- unlist(decomp.out$ncohrt, use.names = FALSE)
          C.mat <- unlist(decomp.out$C.mat, use.names = FALSE)

          #calculates "growth multipliers"
          gmult.out <- gmult(bgs = bgs, egs = egs, availn = availn,
                             degd = degd, dmin = spp.params$DMIN,
                             dmax = spp.params$DMAX, d3 = spp.params$D3, fj = fj,
                             cm1 = spp.params$CM1, cm3 = spp.params$CM3, cm2 = spp.params$CM2,
                             cm4 = spp.params$CM4, cm5 = spp.params$CM5, nspec = nspec)

          smgf <- unlist(gmult.out$smgf, use.names = FALSE) #soil moisture growth factor
          sngf <- unlist(gmult.out$sngf, use.names = FALSE) #soil nitrogen growth factor
          degdgf <- unlist(gmult.out$degdgf, use.names = FALSE) #degree day growth factor
          availn <- unlist(gmult.out$availn, use.names = FALSE) #available nitrogen

          #birth subroutine
          birth.out <- birth(nspec = nspec, ntrees = ntrees, frt = spp.params$FRT, iage = iage,
                             slta = spp.params$SLTA, sltb = spp.params$SLTB, spp.num = spp.params$Spp_Number,dbh = dbh,
                             fwt = spp.params$FWT,
                             degd = degd, dmin = spp.params$DMIN, dmax = spp.params$DMAX,
                             frost = spp.params$FROST, rt = temp.mat[i,1:12], itol = spp.params$ITOL,
                             mplant = spp.params$MPLANT, nogro = nogro,
                             ksprt = ksprt, sprtnd = spp.params$SPRTND, max.ind = max.ind, smgf=smgf,
                             degdgf = degdgf)

          if(is.null(unlist(birth.out$ntrees, use.names = FALSE))){
            ntrees[,i,k] <- rep(0,nspec)
            ntrees.birth[,i,k] <- rep(0,nspec)
          } else {
            ntrees.birth[,i,k] <- unlist(birth.out$ntrees, use.names = FALSE)
            ntrees <- unlist(birth.out$ntrees, use.names = FALSE)
          }

          dbh <- unlist(birth.out$dbh, use.names = FALSE)
          nogro <- unlist(birth.out$nogro, use.names = FALSE)
          ksprt <- unlist(birth.out$ksprt, use.names = FALSE)
          iage <- unlist(birth.out$iage, use.names = FALSE)

          #if(dbh[sum(ntrees)]==0) browser()

          #growth subroutine - increments dbh


          grow.out <- grow.opt(max.ind = max.ind, nspec = nspec, ntrees = ntrees, frt = spp.params$FRT, slta = spp.params$SLTA,
                               sltb = spp.params$SLTB, dbh = dbh, fwt = spp.params$FWT, b2 = spp.params$B2,
                               b3 = spp.params$B3, itol =spp.params$ITOL, g = spp.params$G, degdgf = degdgf,
                               smgf = smgf, sngf= sngf,frost = spp.params$FROST, rt = temp.mat[i,1:12], iage = iage,
                               nogro=nogro)


          if(is.null(unlist(grow.out$ntrees, use.names = FALSE))){
            ntrees <- rep(0,nspec)
            ntrees.grow[,i,k] <- rep(0,nspec)
          } else {
            ntrees.grow[,i,k] <- unlist(grow.out$ntrees, use.names = FALSE)
            ntrees <- unlist(grow.out$ntrees, use.names = FALSE)
          }
          dbh <- unlist(grow.out$dbh, use.names = FALSE)
          awp <- unlist(grow.out$awp, use.names = FALSE)
          nogro <- unlist(grow.out$nogro, use.names = FALSE)
          #browser()
          gf.vec.save[1:nspec,1:4,i,k] <- grow.out$gf.vec
          algf.save.keep[,,i,k] <- grow.out$algf.save
          # gf.vec.save[,i,k] <- gf.vec

          #if(dbh[sum(ntrees)]==0) browser()

          #kill subroutine
          kill.out<- kill(nspec = nspec, ntrees= ntrees,slta = spp.params$SLTA, sltb = spp.params$SLTB,
                          dbh = dbh, agemx = spp.params$AGEMX, ksprt = ksprt,
                          sprtmn = spp.params$SPRTMN, sprtmx = spp.params$SPRTMX, iage  = iage,
                          nogro  = nogro,tl = spp.params$TL,rtst = spp.params$RTST, fwt = spp.params$FWT,
                          max.ind = max.ind, frt = spp.params$FRT)

          ntrees <- unlist(kill.out$ntrees, use.names = FALSE)
          ntrees.kill[,i,k] <- unlist(kill.out$ntrees, use.names = FALSE)
          dbh <- unlist(kill.out$dbh, use.names = FALSE)
          nogro <- unlist(kill.out$nogro, use.names = FALSE)
          ksprt <- unlist(kill.out$ksprt, use.names = FALSE)
          iage <- unlist(kill.out$iage, use.names = FALSE)
          tyl <- unlist(kill.out$tyl, use.names = FALSE)
          tyl.save[,i,k] <- unlist(kill.out$tyl, use.names = FALSE)
          tyl[is.na(tyl)] <- 0

          #output subroutine
          output.out <- output(availn = availn, tyln = tyln, nspec = nspec, frt=spp.params$FRT,
                               iage = iage,slta = spp.params$SLTA, max.ind = max.ind,
                               sltb = spp.params$SLTB,dbh = dbh,fwt = spp.params$FWT,tyl = tyl,
                               ntrees=ntrees,awp=awp)

          #conversion factors
          DEFAULT.C <- 0.48  ## mass percent C of biomass
          PLOT.AREA <- 833 ## m^2
          toKG <- 1000 ## g in Kg
          yearSecs <- (3.15569 * 10^7) ## seconds in a year

          #save variables
          awp.save[1:length(grow.out$awp),i,k] = unlist(grow.out$awp, use.names = FALSE)
          tstem[i,k] = unlist(output.out$atot, use.names = FALSE) #number of stems
          tab[i,k] = unlist(output.out$tbar, use.names = FALSE) #total aboveground biomass
          abvgrnwood[i,k] = unlist(output.out$twbar, use.names = FALSE) #total aboveground biomass
          area[i,k] = unlist(output.out$area, use.names = FALSE)/10 #LAI
          water[i,k] = unlist(moist.out$water, use.names = FALSE) #soil moisture
          fl[i,k] = unlist(kill.out$tyl, use.names = FALSE)[17] #leaf litter
          totl[i,k] = unlist(output.out$tyln, use.names = FALSE) #leaf litter N
          tnap[i,k] = unlist(output.out$tynap, use.names = FALSE) #net aboveground production
          avln[i,k] = unlist(gmult.out$availn, use.names = FALSE) #available nitrogen
          cn[i,k] = unlist(decomp.out$hcn, use.names = FALSE) #humus C:N ratio
          sco2c[i,k] = unlist(decomp.out$sco2, use.names = FALSE) #soil co2 evolution
          som[i,k] = unlist(decomp.out$ff[19,2], use.names = FALSE) #soil organic matter
          bar[,i,k] = unlist(output.out$bar, use.names = FALSE) #species biomass
          aet.save[i,k] = aet #annual evapotranspiration
          nogro.save[,i,k] = unlist(kill.out$nogro, use.names = FALSE)
          dbh.save[,i,k] = unlist(kill.out$dbh, use.names = FALSE)
          iage.save[,i,k] = unlist(kill.out$iage, use.names = FALSE)
          ncohrt.save[i,k] = ncohrt
          npp.spp.save[,i,k] = (unlist(grow.out$npp.spp,use.names=FALSE) * (1 / PLOT.AREA) * (1 / yearSecs) * DEFAULT.C)

          print(paste("year = ",i))
        }
        Rprof(NULL)
        summaryRprof()


        print(paste("PLOT = ",k))
      }

      #unit conversions for variables of interest #need to recheck more carefully later
      year <- seq(1,nyear,1)
      ag.biomass <- (tab  * (1 / PLOT.AREA) * DEFAULT.C) # Above Ground Biomass in kgC/m2 #total aboveground biomass
      abvgroundwood.biomass <- (abvgrnwood  * (1 / PLOT.AREA) * DEFAULT.C) # Above Ground Biomass in kgC/m2 #total aboveground biomass
      total.soil.carbon <- (som + fl) * 907.185 * (1/10000) * DEFAULT.C # TotSoilCarb in kgC/m2
      leaf.litter <- fl * 907.185 * (1/10000) * DEFAULT.C # leaf litter in kgC/m2
      ag.npp <- (tnap * (1 / 10000) * (1 / yearSecs) * DEFAULT.C * 907.185) # GWBI = NPP in linkages
      hetero.resp <- (sco2c * (1/10000) * (1 / yearSecs) * 907.185) # HeteroResp in kgC/m^2/s
      nee <- ((ag.npp - hetero.resp))# NEE #possibly questionable
      et <- aet.save * (1 / yearSecs) # Evap in mm/s
      agb.pft <- (bar  * (1 / PLOT.AREA) * DEFAULT.C) #biomass by PFT
      if(nspec>1){
        f.comp <- t(t(bar[,,1]  * (1 / PLOT.AREA) * DEFAULT.C) / colSums((as.matrix(bar[,,1]) * (1 / PLOT.AREA) * DEFAULT.C))) #f composition
        f.comp[is.na(f.comp)]<-0 #look into prop.table()
      }else{
        f.comp <- matrix(1,nspec,nyear)
      }


      #NOT USED IN CURRENT PECAN OUTPUT #Add? SoilMoisture? LAI? StemDensity?
      #What about MIP stuff?
      #Can we get root biomass from C.mat?
      #tstem[i,k] <- unlist(output.out$atot, use.names = FALSE) #number of stems
      #totl[i,k] = unlist(output.out$tyln, use.names = FALSE) #leaf litter N
      #avln[i,k] = unlist(gmult.out$availn, use.names = FALSE) #available nitrogen
      #cn[i,k] = unlist(decomp.out$hcn, use.names = FALSE) #humus C:N ratio

      output.file <- file.path(outdir,"linkages.out.Rdata")
      sprintf("%s",output.file)

      save(year = year, ag.biomass = ag.biomass, abvgroundwood.biomass =  abvgroundwood.biomass, total.soil.carbon = total.soil.carbon,
           leaf.litter = leaf.litter, ag.npp = ag.npp, hetero.resp = hetero.resp,
           nee = nee, et = et, agb.pft = agb.pft, f.comp = f.comp,
           ntrees.birth = ntrees.birth, ntrees.kill = ntrees.kill, tstem = tstem,
           tab = tab, abvgrnwood=abvgrnwood,fl = fl,totl = totl,tnap = tnap,avln = avln,cn = cn,sco2c = sco2c,
           som = som,bar = bar,aet.save = aet.save,nogro.save = nogro.save,
           dbh.save = dbh.save, iage.save = iage.save, C.mat = C.mat, tyl = tyl,
           ncohrt = ncohrt, area = area, water = water, ksprt = ksprt, tyl.save = tyl.save,
           ff=ff, gf.vec.save = gf.vec.save, algf.save.keep = algf.save.keep,
           npp.spp.save=npp.spp.save, awp.save = awp.save, file = output.file)

      file.exists(output.file)

      output.check(linkages.input,outdir)

    }

  else{

  load(linkages.input) #normal linkages input
  load(spinup_input) #input from the spinup

  nyear = nyear + 1 # for indexing purposes

  #if(spinup == TRUE){
  #  load(spinup_input)
  #  max.ind <- length(dbh)
  #}else{
  max.ind <- 200
  #}

  #temp.mat <- matrix(temp.mat,nyear,12)
  #precip.mat <- matrix(precip.mat,nyear,12)

  ### Reassigning spinup_input values to new names as to not interfere with storage matrices to be created
  assign('gf.vec.save.spinup',gf.vec.save)
  assign('tstem.spinup',tstem)
  assign('area.spinup',area)
  assign('water.spinup',water)
  assign('tab.spinup',tab)
  assign('abvgrnwood.spinup',abvgrnwood)
  assign('fl.spinup',fl)
  assign('totl.spinup',totl)
  assign('tnap.spinup',tnap)
  assign('avln.spinup',avln)
  assign('cn.spinup',cn)
  assign('sco2c.spinup',sco2c)
  assign('som.spinup',som)
  assign('aet.save.spinup',aet.save)
  #assign('ncohrt.save.spinup',ncohrt.save)
  assign('ntrees.birth.spinup',ntrees.birth)
  # assign('ntrees.grow.spinup',ntrees.grow)
  assign('ntrees.kill.spinup',ntrees.kill)
  assign('tyl.save.spinup',tyl.save)
  assign('bar.spinup',bar)
  assign('algf.save.keep.spinup',algf.save.keep)
  assign('nogro.save.spinup',nogro.save)
  assign('dbh.save.spinup',dbh.save)
  assign('iage.save.spinup',iage.save)
  assign('awp.save.spinup',awp.save)
  assign('npp.spp.save.spinup',npp.spp.save)

  #Storage
  gf.vec.save <- array(NA,dim=c(nspec,4,nyear,iplot))
  tstem = matrix(0,nyear,iplot) #number of stems
  area = matrix(0,nyear,iplot)
  water = matrix(0,nyear,iplot)
  tab = matrix(0,nyear,iplot) #total aboveground biomass
  abvgrnwood = matrix(0,nyear,iplot) #total aboveground woody biomass
  fl = matrix(0,nyear,iplot) #leaf litter
  totl = matrix(0,nyear,iplot) #leaf litter N
  tnap = matrix(0,nyear,iplot) #net aboveground production
  avln = matrix(0,nyear,iplot) #available nitrogen
  cn = matrix(0,nyear,iplot) #humus C:N ratio
  sco2c = matrix(0,nyear,iplot) #soil co2 evolution
  som = matrix(0,nyear,iplot) #soil organic matter
  aet.save = matrix(0,nyear,iplot)
  ncohrt.save = matrix(0,nyear,iplot)
  tyl.save = array(0,dim=c(20,nyear,iplot))
  ntrees.birth <- array(0,dim=c(nspec,nyear,iplot))
  ntrees.grow <- array(0,dim=c(nspec,nyear,iplot))
  ntrees.kill <- array(0,dim=c(nspec,nyear,iplot))
  # gf.vec.save <- array(0,dim=c(4,nyear,iplot))
  bar <- array(0,dim=c(nspec,nyear,iplot))
  algf.save.keep<- array(NA,dim=c(max.ind,nspec,nyear,iplot))
  nogro.save <- array(0,dim=c(max.ind,nyear,iplot))
  dbh.save <- array(0,dim=c(max.ind,nyear,iplot))
  iage.save <- array(0,dim=c(max.ind,nyear,iplot))
  awp.save <- array(0,dim=c(max.ind,nyear,iplot))
  npp.spp.save <- array(0,dim=c(nspec,nyear,iplot))

  ### taking inputs from spinup_input and putting at index 1 in storage matrices/arrays
  spinup_len
  ### add to arrays


  gf.vec.save[,,1,] <- gf.vec.save.spinup[,,spinup_len,]
  tstem[1,] <- tstem.spinup[spinup_len,]
  area[1,] <- area.spinup[spinup_len,]
  water[1,] <- water.spinup[spinup_len,]
  tab[1,] <- tab.spinup[spinup_len,]
  abvgrnwood[1,] <- abvgrnwood.spinup[spinup_len,]
  fl[1,] <- fl.spinup[spinup_len,]
  totl[1,] <- totl.spinup[spinup_len,]
  tnap[1,] <- tnap.spinup[spinup_len,]
  avln[1,] <- avln.spinup[spinup_len,]
  cn[1,] <- cn.spinup[spinup_len,]
  sco2c[1,] <- sco2c.spinup[spinup_len,]
  som[1,] <- som.spinup[spinup_len,]
  aet.save[1,] <- aet.save.spinup[spinup_len,]
  #ncohrt.save[1,] <- ncohrt.save.spinup[spinup_len,]
  tyl.save[,1,] <- tyl.save.spinup[,spinup_len,]
  ntrees.birth[,1,] <- ntrees.birth.spinup[,spinup_len,]
  #ntrees.grow[,1,] <- ntrees.grow.spinup[,spinup_len,]
  ntrees.kill[,1,] <- ntrees.kill.spinup[,spinup_len,]
  bar[,1,] <- bar.spinup[,spinup_len,]
  algf.save.keep[,,1,] <- algf.save.keep.spinup[,,spinup_len,]
  nogro.save[,1,] <- nogro.save.spinup[,spinup_len,]
  dbh.save[,1,] <- dbh.save.spinup[,spinup_len,]
  iage.save[,1,] <- iage.save.spinup[,spinup_len,]
  awp.save[,1,] <- awp.save.spinup[,spinup_len,]
  npp.spp.save[,1,] <- npp.spp.save.spinup[,spinup_len,]


  for(k in 1:iplot){ #loop over plots

    #if(spinup == FALSE){

    #plotin.out <- plotin(iplot = k, basesc = basesc, basesn = basesn, max.ind = max.ind,
    #                     nspec = nspec) # initializes storage matrices with zeros for each plot

    #ncohrt <- unlist(plotin.out$ncohrt, use.names = FALSE)
    #tyl <- unlist(plotin.out$tyl, use.names = FALSE)
    #C.mat <- unlist(plotin.out$C.mat, use.names = FALSE)
    #ntrees <- unlist(plotin.out$ntrees, use.names = FALSE)
    #dbh <- unlist(plotin.out$dbh, use.names = FALSE)
    #nogro <- unlist(plotin.out$nogro, use.names = FALSE)
    #ksprt <- unlist(plotin.out$ksprt, use.names = FALSE)
    #iage <- unlist(plotin.out$iage, use.names = FALSE)

    # had to also initialize these values with the last entry in the spinup matrices

    tyl <- tyl.save.spinup[,spinup_len,]
    ntrees <- ntrees.kill.spinup[,spinup_len,]
    dbh <- dbh.save.spinup[,spinup_len,]
    nogro <- nogro.save.spinup[,spinup_len,]
    iage <- iage.save.spinup[,spinup_len,]


    #change loop to start at 2 and go to nyear
    #after that small change nothing else should be different from normal Linkages.R

    #### double check indexing with temp and precip matrix
    ### add in 1899 if it's never called push matrices forward one year to compensate for loop here

    Rprof(interval=.005)
    for(i in 2:nyear){

      #calculates degree days for the year
      degd <- tempe(temp.vec = temp.mat[i-1,1:12])
      #degd = unlist(tempe.out$degd, use.names = FALSE)

      bgs <- (which(temp.mat[i-1,]>10)[1]-1)*30
      egs <- bgs + length(which(temp.mat[i-1,]>10))*30

      #calculates aet
      moist.out <- moist(kyr = i, temp.vec = temp.mat[i-1,1:12], precip.vec = precip.mat[i-1,1:12],
                         fc = fc, dry = dry, bgs = bgs, egs = egs, plat = plat, clat = clat)
      aet <- unlist(moist.out$aet, use.names = FALSE)
      fj <- unlist(moist.out$fj, use.names = FALSE)

      #decomposition subroutine
      decomp.out <- decomp(fdat = fdat, aet = aet,
                           ncohrt = ncohrt, fc = fc, dry = dry,
                           tyl = tyl, C.mat = C.mat)

      ff <- unlist(decomp.out$ff, use.names = FALSE)
      availn <- unlist(decomp.out$availn, use.names = FALSE)
      tyln <- unlist(decomp.out$tyln, use.names = FALSE)
      hcn <- unlist(decomp.out$hcn, use.names = FALSE)
      sco2 <- unlist(decomp.out$sco2, use.names = FALSE)
      ncohrt <- unlist(decomp.out$ncohrt, use.names = FALSE)
      C.mat <- unlist(decomp.out$C.mat, use.names = FALSE)

      #calculates "growth multipliers"
      gmult.out <- gmult(bgs = bgs, egs = egs, availn = availn,
                         degd = degd, dmin = spp.params$DMIN,
                         dmax = spp.params$DMAX, d3 = spp.params$D3, fj = fj,
                         cm1 = spp.params$CM1, cm3 = spp.params$CM3, cm2 = spp.params$CM2,
                         cm4 = spp.params$CM4, cm5 = spp.params$CM5, nspec = nspec)

      smgf <- unlist(gmult.out$smgf, use.names = FALSE) #soil moisture growth factor
      sngf <- unlist(gmult.out$sngf, use.names = FALSE) #soil nitrogen growth factor
      degdgf <- unlist(gmult.out$degdgf, use.names = FALSE) #degree day growth factor
      availn <- unlist(gmult.out$availn, use.names = FALSE) #available nitrogen

      #birth subroutine
      birth.out <- birth(nspec = nspec, ntrees = ntrees, frt = spp.params$FRT, iage = iage,
                         slta = spp.params$SLTA, sltb = spp.params$SLTB, spp.num = spp.params$Spp_Number,dbh = dbh,
                         fwt = spp.params$FWT,
                         degd = degd, dmin = spp.params$DMIN, dmax = spp.params$DMAX,
                         frost = spp.params$FROST, rt = temp.mat[i-1,1:12], itol = spp.params$ITOL,
                         mplant = spp.params$MPLANT, nogro = nogro,
                         ksprt = ksprt, sprtnd = spp.params$SPRTND, max.ind = max.ind, smgf=smgf,
                         degdgf = degdgf)

      if(is.null(unlist(birth.out$ntrees, use.names = FALSE))){
        ntrees[,i,k] <- rep(0,nspec)
        ntrees.birth[,i,k] <- rep(0,nspec)
      } else {
        ntrees.birth[,i,k] <- unlist(birth.out$ntrees, use.names = FALSE)
        ntrees <- unlist(birth.out$ntrees, use.names = FALSE)
      }

      dbh <- unlist(birth.out$dbh, use.names = FALSE)
      nogro <- unlist(birth.out$nogro, use.names = FALSE)
      ksprt <- unlist(birth.out$ksprt, use.names = FALSE)
      iage <- unlist(birth.out$iage, use.names = FALSE)

      #if(dbh[sum(ntrees)]==0) browser()

      #growth subroutine - increments dbh


      grow.out <- grow.opt(max.ind = max.ind, nspec = nspec, ntrees = ntrees, frt = spp.params$FRT, slta = spp.params$SLTA,
                           sltb = spp.params$SLTB, dbh = dbh, fwt = spp.params$FWT, b2 = spp.params$B2,
                           b3 = spp.params$B3, itol =spp.params$ITOL, g = spp.params$G, degdgf = degdgf,
                           smgf = smgf, sngf= sngf,frost = spp.params$FROST, rt = temp.mat[i-1,1:12], iage = iage,
                           nogro=nogro)


      if(is.null(unlist(grow.out$ntrees, use.names = FALSE))){
        ntrees <- rep(0,nspec)
        ntrees.grow[,i,k] <- rep(0,nspec)
      } else {
        ntrees.grow[,i,k] <- unlist(grow.out$ntrees, use.names = FALSE)
        ntrees <- unlist(grow.out$ntrees, use.names = FALSE)
      }
      dbh <- unlist(grow.out$dbh, use.names = FALSE)
      awp <- unlist(grow.out$awp, use.names = FALSE)
      nogro <- unlist(grow.out$nogro, use.names = FALSE)
      #browser()
      gf.vec.save[1:nspec,1:4,i,k] <- grow.out$gf.vec
      algf.save.keep[,,i,k] <- grow.out$algf.save
      # gf.vec.save[,i,k] <- gf.vec

      #if(dbh[sum(ntrees)]==0) browser()

      #kill subroutine
      kill.out<- kill(nspec = nspec, ntrees= ntrees,slta = spp.params$SLTA, sltb = spp.params$SLTB,
                      dbh = dbh, agemx = spp.params$AGEMX, ksprt = ksprt,
                      sprtmn = spp.params$SPRTMN, sprtmx = spp.params$SPRTMX, iage  = iage,
                      nogro  = nogro,tl = spp.params$TL,rtst = spp.params$RTST, fwt = spp.params$FWT,
                      max.ind = max.ind, frt = spp.params$FRT)

      ntrees <- unlist(kill.out$ntrees, use.names = FALSE)
      ntrees.kill[,i,k] <- unlist(kill.out$ntrees, use.names = FALSE)
      dbh <- unlist(kill.out$dbh, use.names = FALSE)
      nogro <- unlist(kill.out$nogro, use.names = FALSE)
      ksprt <- unlist(kill.out$ksprt, use.names = FALSE)
      iage <- unlist(kill.out$iage, use.names = FALSE)
      tyl <- unlist(kill.out$tyl, use.names = FALSE)
      tyl.save[,i,k] <- unlist(kill.out$tyl, use.names = FALSE)
      tyl[is.na(tyl)] <- 0

      #output subroutine
      output.out <- output(availn = availn, tyln = tyln, nspec = nspec, frt=spp.params$FRT,
                           iage = iage,slta = spp.params$SLTA, max.ind = max.ind,
                           sltb = spp.params$SLTB,dbh = dbh,fwt = spp.params$FWT,tyl = tyl,
                           ntrees=ntrees,awp=awp)

      #conversion factors
      DEFAULT.C <- 0.48  ## mass percent C of biomass
      PLOT.AREA <- 833 ## m^2
      toKG <- 1000 ## g in Kg
      yearSecs <- (3.15569 * 10^7) ## seconds in a year

      #save variables
      awp.save[1:length(grow.out$awp),i,k] = unlist(grow.out$awp, use.names = FALSE)
      tstem[i,k] = unlist(output.out$atot, use.names = FALSE) #number of stems
      tab[i,k] = unlist(output.out$tbar, use.names = FALSE) #total aboveground biomass
      abvgrnwood[i,k] = unlist(output.out$twbar, use.names = FALSE) #total aboveground biomass
      area[i,k] = unlist(output.out$area, use.names = FALSE)/10 #LAI
      water[i,k] = unlist(moist.out$water, use.names = FALSE) #soil moisture
      fl[i,k] = unlist(kill.out$tyl, use.names = FALSE)[17] #leaf litter
      totl[i,k] = unlist(output.out$tyln, use.names = FALSE) #leaf litter N
      tnap[i,k] = unlist(output.out$tynap, use.names = FALSE) #net aboveground production
      avln[i,k] = unlist(gmult.out$availn, use.names = FALSE) #available nitrogen
      cn[i,k] = unlist(decomp.out$hcn, use.names = FALSE) #humus C:N ratio
      sco2c[i,k] = unlist(decomp.out$sco2, use.names = FALSE) #soil co2 evolution
      som[i,k] = unlist(decomp.out$ff[19,2], use.names = FALSE) #soil organic matter
      bar[,i,k] = unlist(output.out$bar, use.names = FALSE) #species biomass
      aet.save[i,k] = aet #annual evapotranspiration
      nogro.save[,i,k] = unlist(kill.out$nogro, use.names = FALSE)
      dbh.save[,i,k] = unlist(kill.out$dbh, use.names = FALSE)
      iage.save[,i,k] = unlist(kill.out$iage, use.names = FALSE)
      ncohrt.save[i,k] = ncohrt
      npp.spp.save[,i,k] = (unlist(grow.out$npp.spp,use.names=FALSE) * (1 / PLOT.AREA) * (1 / yearSecs) * DEFAULT.C)

      print(paste("year = ",i))
    }
    Rprof(NULL)
    summaryRprof()


    print(paste("PLOT = ",k))
  }

  #unit conversions for variables of interest #need to recheck more carefully later
  year <- seq(1,nyear,1)
  ag.biomass <- (tab  * (1 / PLOT.AREA) * DEFAULT.C) # Above Ground Biomass in kgC/m2 #total aboveground biomass
  abvgroundwood.biomass <- (abvgrnwood  * (1 / PLOT.AREA) * DEFAULT.C) # Above Ground Biomass in kgC/m2 #total aboveground biomass
  total.soil.carbon <- (som + fl) * 907.185 * (1/10000) * DEFAULT.C # TotSoilCarb in kgC/m2
  leaf.litter <- fl * 907.185 * (1/10000) * DEFAULT.C # leaf litter in kgC/m2
  ag.npp <- (tnap * (1 / 10000) * (1 / yearSecs) * DEFAULT.C * 907.185) # GWBI = NPP in linkages
  hetero.resp <- (sco2c * (1/10000) * (1 / yearSecs) * 907.185) # HeteroResp in kgC/m^2/s
  nee <- ((ag.npp - hetero.resp))# NEE #possibly questionable
  et <- aet.save * (1 / yearSecs) # Evap in mm/s
  agb.pft <- (bar  * (1 / PLOT.AREA) * DEFAULT.C) #biomass by PFT
  if(nspec>1){
    f.comp <- t(t(bar[,,1]  * (1 / PLOT.AREA) * DEFAULT.C) / colSums((as.matrix(bar[,,1]) * (1 / PLOT.AREA) * DEFAULT.C))) #f composition
    f.comp[is.na(f.comp)]<-0 #look into prop.table()
  }else{
    f.comp <- matrix(1,nspec,nyear)
  }


  #NOT USED IN CURRENT PECAN OUTPUT #Add? SoilMoisture? LAI? StemDensity?
  #What about MIP stuff?
  #Can we get root biomass from C.mat?
  #tstem[i,k] <- unlist(output.out$atot, use.names = FALSE) #number of stems
  #totl[i,k] = unlist(output.out$tyln, use.names = FALSE) #leaf litter N
  #avln[i,k] = unlist(gmult.out$availn, use.names = FALSE) #available nitrogen
  #cn[i,k] = unlist(decomp.out$hcn, use.names = FALSE) #humus C:N ratio

  output.file <- file.path(outdir,"linkages.out.Rdata")
  sprintf("%s",output.file)

  save(year = year, ag.biomass = ag.biomass, abvgroundwood.biomass =  abvgroundwood.biomass, total.soil.carbon = total.soil.carbon,
       leaf.litter = leaf.litter, ag.npp = ag.npp, hetero.resp = hetero.resp,
       nee = nee, et = et, agb.pft = agb.pft, f.comp = f.comp,
       ntrees.birth = ntrees.birth, ntrees.kill = ntrees.kill, tstem = tstem,
       tab = tab, abvgrnwood=abvgrnwood,fl = fl,totl = totl,tnap = tnap,avln = avln,cn = cn,sco2c = sco2c,
       som = som,bar = bar,aet.save = aet.save,nogro.save = nogro.save,
       dbh.save = dbh.save, iage.save = iage.save, C.mat = C.mat, tyl = tyl,
       ncohrt = ncohrt, area = area, water = water, ksprt = ksprt, tyl.save = tyl.save,
       ff=ff, gf.vec.save = gf.vec.save, algf.save.keep = algf.save.keep,
       npp.spp.save=npp.spp.save, awp.save = awp.save, file = output.file)

  file.exists(output.file)

  output.check(linkages.input,outdir)

  }
}

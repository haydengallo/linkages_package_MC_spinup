##' @title LINKAGES output function
##' @author Ann Raiho
##'
##' @param availn available nitrogen
##' @param tyln leaf litter N content from "decomp.r"
##' @param nspec number of species
##' @param frt,slta,sltb,fwt species specific parameters
##' @param iage age of each individual
##' @param dbh diameter of each individual
##' @param tyl total yearly litter
##' @param max.ind maximum number of individual
##' @param ntrees number of trees of each species
##' @param awp aboveground woody production
##'
##' @description    OUTPUT coverts model variables into ecosystem variables of interest
##'
##' @return atot total number of trees per ha
##' @return tbar total aboveground biomass
##' @return tyln leaf litter N content
##' @return tynap total aboveground production
##' @return availn available N
##' @return bar specieis biomass t/ha
##' @export
##'
output <- function(availn,tyln,nspec,frt,iage,slta,sltb,dbh,fwt,tyl,max.ind,ntrees,awp,bio_method,spp.num = spp.params$Spp_Number){

  #initialization
  area = 0 #leaf area
  folw = 0 #leaf biomass
  availn = availn*1000 #available nitrogen from "gmult.r"
  tbar = 0 #total aboveground biomass
  twbar = 0 #total aboveground wood biomass
  tawp = 0 #total aboveground woody production
  ntot = 0 #number of trees
  tyln = tyln*1000 #leaf litter N content from "decomp.r"

  bar = numeric(nspec)
  abvgrnwood = numeric(nspec)

  #calculate spp biomass, total biomass, total number of stems, leaf area, and total woody production
  nl = 1
  for(i in 1:nspec){

    # initialize
    if(bio_method == 'lambert'){
      y_wood = 0
      y_bark = 0
      y_foliage = 0
      y_branches = 0
      tot = 0
    }
    bar[i] = 0
    abvgrnwood[i] = 0
    if(ntrees[i]==0) next
    nu = nl + ntrees[i] - 1
    ret = frt[i]

    code = spp.num[i]

    # add coefficients
    if(bio_method == 'chojnacky'){
      if(code == 'ACRU'){b_0 = -2.047; b_1 = 2.3852} # acer rubrum, red maple
      if(code == 'ACSA3'){b_0 = -1.8011; b_1 = 2.3852} # acer saccharum, sugar maple
      if(code == 'BEAL2'){b_0 = -1.8096; b_1 = 2.348} # betula alleghaniensis, yellow birch
      if(code == 'BELE'){b_0 = -2.2652; b_1 = 2.5349} # betula lenta, sweet birch
      if(code == 'FAGR'){b_0 = -2.0705; b_1 = 2.441} # fagus grandifolia, american beech
      if(code == 'PIRU'){b_0 = -2.1364; b_1 = 2.3233} # picea rubens, red spruce
      if(code == 'PIST'){b_0 = -2.6177; b_1 = 2.4638} # pinus strobus, white pine
      if(code == 'QUAL'){b_0 = -2.0705; b_1 = 2.441} # quercus alba, white oak
      if(code == 'QUMO'){b_0 = -2.0705; b_1 = 2.441} # quercus montana, chestnut oak
      if(code == 'QURU'){b_0 = -2.0705; b_1 = 2.441} # quercus rubra, red oak
      if(code == 'QUVE'){b_0 = -2.0705; b_1 = 2.441} # quercus velutina, black oak
      if(code == 'THOC2'){b_0 = -1.9615; b_1 = 2.1063} # thuja occidentalis, northern white-cedar
      if(code == 'TSCA'){b_0 = -2.348; b_1 = 2.3876} # tsuga canadensis, eastern hemlock
    }

    if(bio_method == 'lambert'){
      # acer rubrum, red maple
      if(code == 'ACRU'){bwood1 = 0.1014; bwood2 = 2.3448; bbark1 = 0.0291; bbark2 = 2.0893;
                         bbranches1 = 0.0175; bbranches2 = 2.4846; bfoliage1 = 0.0515; bfoliage2 = 1.5198}
      # acer saccharum, sugar maple
      if(code == 'ACSA3'){bwood1 = 0.1315; bwood2 = 2.3129; bbark1 = 0.0631; bbark2 = 1.6241;
                          bbranches1 = 0.033; bbranches2 = 2.3741; bfoliage1 = 0.0393; bfoliage2 = 1.693}
      # betula alleghaniensis, yellow  birch
      if(code == 'BEAL2'){bwood1 = 0.1932; bwood2 = 2.1569; bbark1 = 0.0192; bbark2 = 2.2475;
                          bbranches1 = 0.0305; bbranches2 = 2.4044; bfoliage1 = 0.1119; bfoliage2 = 1.3973}
      # fagus grandifolia, american beech
      if(code == 'FAGR'){bwood1 = 0.1478; bwood2 = 2.2986; bbark1 = 0.012; bbark2 = 2.2388;
                         bbranches1 = 0.037; bbranches2 = 2.368; bfoliage1 = 0.0376; bfoliage2 = 1.6164}
      # picea rubens, red spruce
      if(code == 'PIRU'){bwood1 = 0.0989; bwood2 = 2.2814; bbark1 = 0.022; bbark2 = 2.0908;
                         bbranches1 = 0.0005; bbranches2 = 3.275; bfoliage1 = 0.0066; bfoliage2 = 2.4213}
      # pinus strobus, white pine
      if(code == 'PIST'){bwood1 = 0.0997; bwood2 = 2.2709; bbark1 = 0.0192; bbark2 = 2.2038;
                         bbranches1 = 0.0056; bbranches2 = 2.6011; bfoliage1 = 0.0284; bfoliage2 = 1.9375}
      # quercus alba, white oak
      if(code == 'QUAL'){bwood1 = 0.0762; bwood2 = 2.3335; bbark1 = 0.0338; bbark2 = 1.9845;
                         bbranches1 = 0.0113; bbranches2 = 2.6211; bfoliage1 = 0.0188; bfoliage2 = 1.7881}
      # quercus rubra, red oak
      if(code == 'QURU'){bwood1 = 0.1754; bwood2 = 2.1616; bbark1 = 0.0381; bbark2 = 2.0991;
                         bbranches1 = 0.0085; bbranches2 = 2.779; bfoliage1 = 0.0373; bfoliage2 = 1.674}
      # tsuga canadensis, eastern hemlock
      if(code == 'TSCA'){bwood1 = 0.0619; bwood2 = 2.3821; bbark1 = 0.0139; bbark2 = 2.3282;
                         bbranches1 = 0.0217; bbranches2 = 2.2653; bfoliage1 = 0.0776; bfoliage2 = 1.6995}
      # betula lenta, sweet birch (used values for hop-hornbeam, classification from chojnacky)
      if(code == 'BELE'){bwood1 = 0.1929; bwood2 = 1.9672; bbark1 = 0.0671; bbark2 = 1.5911;
                         bbraches1 = 0.0278; bbranches2 = 2.1336; bfoliage1 = 0.0293; bfoliage2 = 1.9502}
      # quercus montana, chestnut oak (used values for white oak, classification from chojnacky)
      if(code == 'QUMO'){bwood1 = 0.0762; bwood2 = 2.3335; bbark1 = 0.0338; bbark2 = 1.9845;
                         bbranches1 = 0.0113; bbranches2 = 2.6211; bfoliage1 = 0.0188; bfoliage2 = 1.7881}
      # quercus velutina, black oak (used values for white oak, classification from chojnacky)
      if(code == 'QUVE'){bwood1 = 0.0762; bwood2 = 2.3335; bbark1 = 0.0338; bbark2 = 1.9845;
                         bbranches1 = 0.0113; bbranches2 = 2.6211; bfoliage1 = 0.0188; bfoliage2 = 1.7881}
      # thuja occidentalis, northern white-cedar
      if(code == 'THOC2'){bwood1 = 0.0654; bwood2 = 2.2121; bbark1 = 0.0114; bbark2 = 2.1432;
                          bbranches1 = 0.0335; bbranches2 = 1.9367; bfoliage1 = 0.0499; bfoliage2 = 1.7278}
    }
    for(j in nl:nu){
      age = iage[j]

      if(bio_method == 'default'){
        #calculate leaf biomass (kg/tree)
        folw = ((slta[i]+sltb[i]*dbh[j])/2)^2 * 3.14 * fwt[i] * ret * .001
        #calculate species biomass (kg/plot)
        bar[i] = bar[i] + .1193 * dbh[j]^2.393 + folw
      }

      if(bio_method == 'lambert'){
        #Lambert allometry
        y_wood[j] = bwood1*(dbh[j]^bwood2)
        y_bark[j] = bbark1*(dbh[j]^bbark2)
        y_foliage[j] = bfoliage1*(dbh[j]^bfoliage2)
        y_branches[j] = bbranches1*(dbh[j]^bbranches2)
        tot[j] = y_wood[j] + y_bark[j] + y_foliage[j] + y_branches[j]
        bar[i] = bar[i] + tot[j]
      }

      if(bio_method == 'chojnacky'){
        #Chonjnacky allometry species biomass (kg/plot)
        bar[i] = bar[i] + (exp(b_0[i])*(dbh[j]^b_1[i]))
      }

      #calculate species biomass (kg/plot)
      #bar[i] = bar[i] + .1193 * dbh[j]^2.393 + folw

      if(dbh[j]>10) {
        abvgrnwood[i] = abvgrnwood[i] + .1193 * dbh[j]^2.393
      }
      if(is.na(bar[i])) bar[i] <- 0

      #calculate leaf area index
      area = area + 1.9283295 * 10^-4 * dbh[j]^2.129

      #calculate woody production (kg/plot)
      tawp = tawp + awp[j]
    }
    #calculate total aboveground biomass (kg/plot)
    tbar = tbar + bar[i]
    twbar = twbar + abvgrnwood[i]
    nl = nu + 1
    #calculate total number of trees per plot
    ntot = ntot + ntrees[i]
    if(ntot > max.ind) print("too many trees -- output")
  }

  #convert number of treees per plot to number per ha
  atot = ntot
  atot = atot*12

  #convert total aboveground biomass and woody production to t/ha
  #tbar = tbar * .012
  tawp = tawp * .012

  #calculate total aboveground production
  tynap = tawp + tyl[17]

  #convert spp biomass to t/ha
  #bar = bar * .012

  return(list(atot=atot,tbar=tbar,twbar=twbar,tyln=tyln,tynap=tynap,availn=availn,
              bar=bar,area=area))

}

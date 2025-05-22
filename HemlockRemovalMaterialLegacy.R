### Data from Audrey Barker Plotkin, Harvard Forest (HFR) ###
### input data here is archived in less-processed form as follows: ###
## https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-hfr.125.25 (Dead wood)
## https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-hfr.106.35 (saplings; 106-05)
## https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-hfr.126.19 (trees; biomass calculated using allometries cited in Barker Plotkin et al. 2024)
## see more information & background in Barker Plotkin et al. 2024 (DOI:10.1002/eap.2957)

### predictor variable: volume or mass of dead wood (the material legacy). 
## Can be a total or broken down into standing dead wood, CWD (downed wood >7.5cm diameter), FWD (downed wood 0.6-7.5 cm diameter)
## input file = DeadByPlot.csv

### response variable: number or biomass of trees (trees are stems >5 cm dbh), and/or number of tree saplings (at least 1.3 m tall but < 5cm dbh)
## Can be all trees together, or broken down into just HEMLOCK and/or just BIRCH (B. lenta + B. alleghaniensis)
## input file = TreeSapDensSpecies.csv

library(tidyverse)
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis/Datasets/Harvard Forest') #FILL IN YOUR OWN WORKING DIRECTORY!

# Dead wood sampled in 2005, 2007, 2009, 2011, 2013, 2015, 2017, 2021 (next 2025).
# Dead wood amounts inferred for 2004, 2014, and 2019 (to line up with tree census dates)
predictorsAll<-read.csv('DeadByPlot.csv', header=TRUE)
str(predictorsAll)
predictorsAll<-predictorsAll %>% mutate_at(c(2,3,5), as.factor)

# Trees (>5 cm dbh) measured in 2004, 2009, 2014, 2019, 2024
# Sapling trees (>1.3 m tall but <5 cm dbh) measured 2004, 2007, 2009, 2011, 2013, 2015, 2017, 2019, 2021, 2023
responseAll<-read.csv('TreeSapDensSpecies.csv', header=TRUE)
str(responseAll)
responseAll$species<-ifelse(is.na(responseAll$species), 'none', responseAll$species)
responseAll<-responseAll %>% mutate_at(c(2,3,5,8), as.factor)

#### Data Wrangling ####
levels(responseAll$species)
# group trees/saplings as hemlock (initial dominant), birch (post-disturbance dominant), and other
responseAll$spgroup<-ifelse(responseAll$species == 'TSCA', 'hemlock', ifelse(responseAll$species %in% c('BELE', 'BEAL'), 'birch', 'other'))

responseL<-responseAll %>% 
  group_by (plot, trt, block, year, spgroup, stratum) %>%
  summarise (dens_ha = sum(dens.ha),
             biomass_gm2=sum(biomass.gm2, na.rm=TRUE))

responseW<-responseL %>%
  pivot_wider(names_from = spgroup, values_from = c(dens_ha, biomass_gm2))
responseW[is.na(responseW)] <- 0 #not exactly since biomass is NA for saplings
responseW$dens_ha_total<-responseW$dens_ha_hemlock +responseW$dens_ha_birch + responseW$dens_ha_other 
responseW$biomass_gm2_total<-responseW$biomass_gm2_hemlock +responseW$biomass_gm2_birch + responseW$biomass_gm2_other 


predictors<-predictorsAll[, c(1:7)] #the notes mess up the pivot
predictorsW<-predictors %>%
  pivot_wider(names_from = group, values_from = c(volm3ha, massgm2))

#use this data frame for looking at matched response & predictor at each time-point
pred.resp <- predictorsW %>% inner_join( responseW, 
                                           by=c('plot','year', 'trt', 'block'))

#another data frame to see how response rate is affected by INITIAL input of legacy material (let's pick 2007 since that is when the girdled trees would all have died)
predictor07<-subset(predictorsW, predictorsW$year == 2007)
predictor07<-predictor07[,c(1:3, 5:12)]
pred.respINIT<- predictor07 %>% inner_join( responseW, 
                                             by=c('plot', 'trt', 'block'))


# okay, these are workable if messy data frames

#### exploring the wrangled data ####

# Q1: does dead wood amount in some form (the legacy) relate to sapling density or tree biomass (the resilience)?
# it probably doesn't make sense to include the non-manipulated treatments since the question is response to disturbance

pred.resp$trt <- ordered(pred.resp$trt, levels = c("hemlock", "girdled", "logged", "hardwood"))

ggplot(subset(pred.resp, stratum=='Tree'), aes(x=massgm2_TotalDead, y=biomass_gm2_total, colour=trt)) + 
  geom_point(aes(shape=block), size=4) +
  theme_bw(base_size = 18) + 
  scale_colour_manual(values=c("#440154FF" ,"#31688EFF","#35B779FF","#FDE725FF"), name="treatment")+
  labs(x="Deadwood mass (gm2)", y="Tree Biomass", title="")

ggplot(subset(pred.resp, stratum=='Sapling'), aes(x=massgm2_TotalDead, y=dens_ha_total, colour=trt)) + 
  geom_point(aes(shape=block), size=4) +
  theme_bw(base_size = 18) + 
  scale_colour_manual(values=c("#440154FF" ,"#31688EFF","#35B779FF","#FDE725FF"), name="treatment")+
  labs(x="Deadwood mass (gm2)", y="Sapling Density", title="")


# Q2: does the initial input of material legacy (dead wood amount in 2009 since the girdled trees weren't all dead in 2009) affect the RATE of recovery (e.g. gain in tree biomass or sapling number)
# for this question, the non-manipulated treatments MIGHT be relevant to include?
pred.respINIT$trt <- ordered(pred.respINIT$trt, levels = c("hemlock", "girdled", "logged", "hardwood"))

# can also look at biomass of just hemlock or birch; or density (stem numbers) for all trees, birch, or hemlock
# to me, the total biomass seems the most promising
pred.respINIT$time<-pred.respINIT$year-2005

ggplot(subset(pred.respINIT, stratum=='Tree'), aes(x=year, y=biomass_gm2_total, colour=trt)) + 
  geom_point(aes(shape=block), size=4) +
  theme_bw(base_size = 18) + 
  scale_colour_manual(values=c("#440154FF" ,"#31688EFF","#35B779FF","#FDE725FF"), name="treatment")+
  labs(x="Year", y="Total tree biomass", title="")

# no biomass for saplings, but can look at density of all trees, or just hemlocks or birches
ggplot(subset(pred.respINIT, stratum=='Sapling'), aes(x=year, y=dens_ha_total, colour=trt)) + 
  geom_point(aes(shape=block), size=4) +
  theme_bw(base_size = 18) + 
  scale_colour_manual(values=c("#440154FF" ,"#31688EFF","#35B779FF","#FDE725FF"), name="treatment")+
  labs(x="Year", y="Sapling Density", title="")

# if this seems promising, then what you'd do next is figure out the rate of change from these and then plot the rate vs initial deadwood
# note that some/most of these are not linear changes over time (biomass )
# simple-minded lm
# look at each plot separately

GRbiomasstree<-lm(biomass_gm2_total~time, data=subset(pred.respINIT, pred.respINIT$time>0 & pred.respINIT$trt=='girdled' & pred.respINIT$stratum=='Tree' & pred.respINIT$block=='ridge'))
summary(GRbiomasstree)
# slope sig 229
GVbiomasstree<-lm(biomass_gm2_total~time, data=subset(pred.respINIT, pred.respINIT$time>0 & pred.respINIT$trt=='girdled' & pred.respINIT$stratum=='Tree' & pred.respINIT$block=='valley'))
summary(GVbiomasstree)
# slope sig 229

LRbiomasstree<-lm(biomass_gm2_total~time, data=subset(pred.respINIT, pred.respINIT$time>0 & pred.respINIT$trt=='logged' & pred.respINIT$stratum=='Tree' & pred.respINIT$block=='ridge'))
summary(LRbiomasstree)
#slope sig 331

LVbiomasstree<-lm(biomass_gm2_total~time, data=subset(pred.respINIT, pred.respINIT$time>0 & pred.respINIT$trt=='logged' & pred.respINIT$stratum=='Tree' & pred.respINIT$block=='valley'))
summary(LVbiomasstree)
#slope sig 357

trt<-c('girdled', 'logged', 'logged', 'girdled')
block<-c('valley', 'valley', 'ridge', 'ridge')
initdead<-subset(predictor07, predictor07$trt %in% c('girdled', 'logged'))
Treebiomassslope<-c(229, 357, 331, 229)
trythis<-data.frame(trt, block, initdead, Treebiomassslope)

ggplot(trythis, aes(x=massgm2_TotalDead, y=Treebiomassslope, colour=trt)) + 
  geom_point(size=4) +
  theme_bw(base_size = 18) + 
  scale_colour_manual(values=c("#31688EFF","#35B779FF"), name="treatment")+
  labs(x="Initial deadwood mass input (LEGACY)", y="Tree biomass recovery (RESILIENCE)", title="")
#great --- but I think the rate of recovery may not really be driven by the quantity of deadwood. . . 


### Another Option -- just the years in which trees were censused (2004, 2009, 2014, 2019, 2024) ####
#fudge time since disturbance so that 2013 and 2014 are the same, and 2023 and 2024 are the same to include saplings
predictorsSUB<-subset(predictorsAll, predictorsAll$year %in% c(2004, 2005, 2009, 2014, 2019)) #we will be measuring dead wood in fall 2025 so don't have an update to match with the 2024 tree data yet
responsesSUB<-subset(responseW, responseW$year %in% c(2004, 2005, 2009, 2013, 2014, 2019, 2023, 2024))  
year<-c(2004, 2005, 2009, 2013, 2014, 2019, 2023, 2024)
timefudge<-c(-1, 0.5, 4, 9, 9, 14, 19, 19)
yeartime<-data.frame(year, timefudge)

responseSUB<-merge(responsesSUB, yeartime, all.y=FALSE)
predictorSUB<-merge(predictorsSUB, yeartime, all.y=FALSE)

pred.respSUB<- predictorsSUB %>% inner_join( responseSUB, 
                                             by=c('plot','timefudge', 'trt', 'block'))

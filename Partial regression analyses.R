#### Partial regression analyses for observational datasets ####

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis/Datasets/Partial regression/')

#### Virginia Coast Reserve ####

# Read in Virginia Coast Reserve dataset
vcr.part_reg <- read_csv("vcr.part_reg.csv")

# First plot relationship between adult density and juvenile density
ggplot(vcr.part_reg, aes(x = adult.mean, y = juv.mean)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab('Adult density') +
  ylab('Juvenile density') +
  theme(aspect.ratio = 1)

# Model this relationship, then extract residuals
vcr.lm.demo_live <- lm(juv.mean ~ adult.mean, data = vcr.part_reg)
summary(vcr.lm.demo_live)

vcr.part_reg$vcr.resid.demo_live <- resid(vcr.lm.demo_live)

# Next plot relationship between adult density and dead density
ggplot(vcr.part_reg, aes(x = adult.mean, y = dead.mean)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab('Adult density') +
  ylab('Dead density') +
  theme(aspect.ratio = 1)

# Model this relationship and extract residuals
vcr.lm.dead_live <- lm(dead.mean ~ adult.mean, data = vcr.part_reg)
summary(vcr.lm.dead_live)

vcr.part_reg$vcr.resid.dead_live <- resid(vcr.lm.dead_live)

# Now plot residuals of the live-demographic relationship against residuals of the live-dead relationship
ggplot(vcr.part_reg, aes(x = vcr.resid.dead_live, y = vcr.resid.demo_live)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme(aspect.ratio = 1) +
  ggtitle('Virginia Coast Reserve (VCR)')+  
  xlab('Resid. dead density\n(controlling for adult density)') +
  ylab('Resid. juvenile density\n(controlling for adult density)') +
  theme_minimal()

# Finally, model this relationship
vcr.lm.demo_dead.resid <- lm(vcr.resid.demo_live ~ vcr.resid.dead_live, data = vcr.part_reg)
summary(vcr.lm.demo_dead.resid)


#### Andrews forest ####

# Read in H. J. Andrews dataset
hja.part_reg <- read_csv("hja.part_reg.csv")

# First plot relationship between live tree basal area and growth of live trees
ggplot(hja.part_reg, aes(x = baph0_spp, y = tree_growth_ind)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab('Live tree basal area') +
  ylab('Tree growth rate') +
  theme(aspect.ratio = 1)

# Model this relationship, then extract residuals
hja.lm.demo_live <- lm(tree_growth_ind ~ baph0_spp, data = hja.part_reg)
summary(hja.lm.demo_live)

hja.part_reg$hja.resid.demo_live <- resid(hja.lm.demo_live)

# Next plot relationship between live tree basal area and mass of dead wood
ggplot(hja.part_reg, aes(x = baph0_spp, y = dw_mass_ha)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab('Live tree basal area') +
  ylab('Mass of deadwood') +
  theme(aspect.ratio = 1)

# Model this relationship and extract residuals
hja.lm.dead_live <- lm(dw_mass_ha ~ baph0_spp, data = hja.part_reg)
summary(hja.lm.dead_live)

hja.part_reg$hja.resid.dead_live <- resid(hja.lm.dead_live)

# Now plot residuals of the live-demographic relationship against residuals of the live-dead relationship
ggplot(hja.part_reg, aes(x = hja.resid.dead_live, y = hja.resid.demo_live)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme(aspect.ratio = 1) +
  ggtitle('Andrews Forest (AND)')+  
  xlab('Resid. deadwood mass\n(controlling for live basal area)') +
  ylab('Resid. tree growth\n(controlling for ive basal area)') +
  theme_minimal()

# Finally, model this relationship
hja.lm.demo_dead.resid <- lm(hja.resid.demo_live ~ hja.resid.dead_live, data = hja.part_reg)
summary(hja.lm.demo_dead.resid)

#### Florida Coastal Everglades ####

# Read in Florida Coastal Everglades dataset
fce.part_reg <- read_csv("fce.part_reg.csv")

# First plot relationship between root biomass and root production
ggplot(fce.part_reg, aes(x = root_biomass.mean, y = root_prod.mean)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab('Mean root biomass') +
  ylab('Mean root production') +
  theme(aspect.ratio = 1)

# Model this relationship, then extract residuals
fce.lm.demo_live <- lm(root_prod.mean ~ root_biomass.mean, data = fce.part_reg)
summary(fce.lm.demo_live)

fce.part_reg$fce.resid.demo_live <- resid(fce.lm.demo_live)

# Next plot relationship between live tree basal area and mass of dead wood
ggplot(fce.part_reg, aes(x = root_biomass.mean, y = mean_litter)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab('Mean root biomass') +
  ylab('Mean litter mass') +
  theme(aspect.ratio = 1)

# Model this relationship and extract residuals
fce.lm.dead_live <- lm(mean_litter ~ root_biomass.mean, data = fce.part_reg)
summary(fce.lm.dead_live)

fce.part_reg$fce.resid.dead_live <- resid(fce.lm.dead_live)

# Now plot residuals of the live-demographic relationship against residuals of the live-dead relationship
ggplot(fce.part_reg, aes(x = fce.resid.dead_live, y = fce.resid.demo_live)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme(aspect.ratio = 1) +
  ggtitle('Florida Coastal Everglades (FCE)')+  
  xlab('Resid. litter mass\n(controlling for root biomass)') +
  ylab('Resid. root production\n(controlling for root biomass)') +
  theme_minimal()

# Finally, model this relationship
fce.lm.demo_dead.resid <- lm(fce.resid.demo_live ~ fce.resid.dead_live, data = fce.part_reg)
summary(fce.lm.demo_dead.resid)

#### Bonanza Creek ####

# Read in Bonanza Creek seed rain dataset
bnz.part_reg <- read_csv("bnz.part_reg.csv")

# First plot relationship between pre-fire total basal area of trees and seed rain density
ggplot(bnz.part_reg, aes(x = bs_ba, y = total_m2)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab('Pre-fire tree basal area') +
  ylab('Seed rain density') +
  theme(aspect.ratio = 1)

# Model this relationship, then extract residuals
bnz.lm.demo_live <- lm(total_m2 ~ bs_ba, data = bnz.part_reg)
summary(bnz.lm.demo_live)

bnz.part_reg$bnz.resid.demo_live <- resid(bnz.lm.demo_live)

# Next plot relationship between live tree basal area and mass of dead wood
# Force intercept to be zero?
ggplot(bnz.part_reg, aes(x = bs_ba, y = bs_stg_ba)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlab('Pre-fire tree basal area') +
  ylab('Standing burned tree basal area') +
  theme(aspect.ratio = 1)

# Model this relationship and extract residuals
bnz.lm.dead_live <- lm(bs_stg_ba ~ bs_ba, data = bnz.part_reg)
summary(bnz.lm.dead_live)

bnz.part_reg$bnz.resid.dead_live <- resid(bnz.lm.dead_live)

# Now plot residuals of the live-demographic relationship against residuals of the live-dead relationship
ggplot(bnz.part_reg, aes(x = bnz.resid.dead_live, y = bnz.resid.demo_live)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme(aspect.ratio = 1) +
  ggtitle('Bonanza Creek (BNZ)')+  
  xlab('Resid. Standing burned tree basal area\n(controlling for pre-fire tree basal area)') +
  ylab('Resid. Seed rain density\n(controlling for pre-fire tree basal area)') +
  theme_minimal()

# Finally, model this relationship
bnz.lm.demo_dead.resid <- lm(bnz.resid.demo_live ~ bnz.resid.dead_live, data = bnz.part_reg)
summary(bnz.lm.demo_dead.resid)

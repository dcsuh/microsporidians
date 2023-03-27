


#load packages
library(tidyverse)
library(magrittr)
library(patchwork)

#read data
dat <- read_csv("/Users/dcsuh/Documents/GitHub/microsporidians/data/GAponds.master.summary.csv")
yield <- read_csv("/Users/dcsuh/Documents/GitHub/microsporidians/workshop/microsporidians.csv")

dat %<>% mutate(total_abund = dpa.ncount + dam.ncount + dla.ncount + sim.ncount + dia.ncount + cer.ncount + bos.ncount) #make total abundance

yield %<>% select(Species, Spores) # get spore yield for each species

yield %<>% mutate(species = case_when(Species == "Am" ~ "dam",
                                      Species == "Amb" ~ "dam",
                                      Species == "Bos" ~ "bos",
                                      Species == "Cerio" ~ "cer",
                                      Species == "Dia" ~ "dia",
                                      Species == "La" ~ "dla",
                                      Species == "Lav" ~ "dla",
                                      Species == "Pa" ~ "dpa",
                                      Species == "Simo" ~ "sim"))

yield %<>% group_by(species) %>% summarize(mean_spores = mean(Spores, na.rm=T),
                                           se_spores = sd(Spores, na.rm=T)/sqrt(sum(!is.na(Spores))),
                                           n = sum(!is.na(Spores)))

# get data for each sampling event
dat1 <- dat %>% select(lakeday, year, lake_id, jday, jday_cum, 
                       dpa.ncount, dam.ncount, dla.ncount, sim.ncount, dia.ncount, cer.ncount, bos.ncount,
                       dpa.ninfect, dam.ninfect, dla.ninfect, sim.ninfect, dia.ninfect, cer.ninfect, bos.ninfect, 
                       dpa.dens, dam.dens, dla.dens, sim.dens, dia.dens, cer.dens, bos.dens, total_abund) %>% distinct() 

#calculate relative abundance and prevalence
dat1 %<>% mutate(dpa.prev = dpa.ninfect/dpa.ncount,
                 dam.prev = dam.ninfect/dam.ncount,
                 dla.prev = dla.ninfect/dla.ncount,
                 sim.prev = sim.ninfect/sim.ncount,
                 dia.prev = dia.ninfect/dia.ncount,
                 cer.prev = cer.ninfect/cer.ncount,
                 bos.prev = bos.ninfect/bos.ncount,
                 dpa.ra = dpa.ncount/total_abund,
                 dam.ra = dam.ncount/total_abund,
                 dla.ra = dla.ncount/total_abund,
                 sim.ra = sim.ncount/total_abund,
                 dia.ra = dia.ncount/total_abund,
                 cer.ra = cer.ncount/total_abund,
                 bos.ra = bos.ncount/total_abund)

#site x event x species
dens <- dat1 %>% pivot_longer(., cols = ends_with("dens"), values_to = "dens") %>% mutate(species = gsub(".dens","",name)) %>% select(lakeday, lake_id, jday_cum, dens, species)

abund <- dat1 %>% pivot_longer(., cols = ends_with("ncount"), values_to = "AB") %>% mutate(species = gsub(".ncount","",name)) %>% select(lakeday, lake_id, jday_cum, AB, species)

rel_abund <- dat1 %>% pivot_longer(., cols = ends_with("ra"), values_to = "RA") %>% mutate(species = gsub(".ra","",name)) %>% select(lakeday, lake_id, jday_cum, RA, species)

prev <- dat1 %>% pivot_longer(., cols = ends_with("prev"), values_to = "prev") %>% mutate(species = gsub(".prev","",name)) %>% select(lakeday, lake_id, jday_cum, prev, species)



#merge
lake_int <- left_join(yield, rel_abund)
lake_int %<>% left_join(., abund)
lake_int %<>% left_join(., dens)
lake_int %<>% left_join(., prev)

#calculate # of spores contributed by each species at each site at each sampling event
lake_int %<>% mutate(rel_spores = ifelse(!is.na(prev),dens*prev*mean_spores, 0))

#calculate total # of spores produced by infected hosts of all species at each site at each sampling event
total_lake_int <- lake_int %>% group_by(lakeday, lake_id, jday_cum) %>% summarize(tot_spores = sum(rel_spores, na.rm = T))

#merge
lake_int %<>% left_join(., total_lake_int)
#dat1 %<>% left_join(.,total_lake_int)

# #look at relationship between total host abundance (infectd + uninfected) and total spores at each site at each sampling event
# dat1 %>% ggplot(.,aes(x=log(total_abund), y=log(tot_spores))) + geom_point() + geom_smooth(method = "lm")
# 
# #relative contribution of each species at each site over all time
# lake_int %>% ggplot(.,aes(x = lake_id, y=rel_spores, fill = species)) + geom_bar(stat = "identity", position = "stack")
# 
# #relative contribution of each species at each site at each sampling event
# lake_int %>% ggplot(., aes(x=rel_spores, y=tot_spores, color=species, size = AB*prev)) + geom_point()


# lake_int %>% filter(!lake_id %in% c("NN1", "Chapman", "Memorial", "Oglethorpe", "NN3", "Sister 1")) %>%
#   filter(!is.na(lake_id)) %>%
#   filter(tot_spores > 0) %>% 
#   ggplot(., aes(x = as.factor(jday_cum), y = log(rel_spores), fill = species)) + 
#   geom_bar(stat = "identity", position = "stack") + 
#   facet_wrap(vars(lake_id),
#              ncol=1)

# 
# lake_int %>% filter(!lake_id %in% c("NN1", "Chapman", "Memorial", "Oglethorpe", "NN3", "Sister 1", "Catfish")) %>%
#   filter(!is.na(lake_id)) %>%
#   filter(jday_cum %in% c(90:260)) %>% 
#   filter(!is.na(tot_spores)) %>% 
#   ggplot(., aes(x = as.factor(jday_cum), y = rel_spores, fill = species)) + 
#   geom_bar(stat = "identity", position = "stack") + 
#   facet_wrap(vars(lake_id),
#              ncol=1)

# 
# lake_int %>% filter(!lake_id %in% c("NN1", "Chapman", "Memorial", "Oglethorpe")) %>%
#   filter(jday_cum %in% c(420:500)) %>% 
#   filter(!is.na(tot_spores)) %>% 
#   ggplot(., aes(x = as.factor(jday_cum), y = rel_spores, fill = species)) + 
#   geom_bar(stat = "identity", position = "stack") + 
#   facet_wrap(vars(lake_id),
#              ncol=1)


#lake_int %>% ggplot(.,aes(x = AB*prev, y = tot_spores, color = species)) + geom_point()
#lake_int %>% ggplot(.,aes(x = mean_spores, y = rel_spores, color = species)) + geom_point()

yield %>% arrange(desc(mean_spores)) %>% mutate(species=factor(species, levels=unique(species))) %>%
  filter(!is.na(mean_spores)) %>% 
  filter(mean_spores > 0) %>%
  ggplot(.,aes(x=species, y=mean_spores)) + geom_point() + geom_linerange(aes(ymin = mean_spores - se_spores, 
                                                                              ymax = mean_spores + se_spores))




p1 <- lake_int %>% arrange(desc(tot_spores)) %>% mutate(lakeday=factor(lakeday, levels=unique(lakeday))) %>%
  filter(species != "bos") %>% 
  filter(tot_spores > 100000000) %>% 
  filter(!lake_id %in% c("NN1", "Chapman", "Memorial", "Oglethorpe")) %>% 
  ggplot(., aes(x = lakeday, y = rel_spores, fill = species)) + 
  geom_bar(stat = "identity", position = "stack") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) + labs(y = "# of Spores")

p2 <- lake_int %>% arrange(desc(tot_spores)) %>% mutate(lakeday=factor(lakeday, levels=unique(lakeday))) %>%
  filter(species != "bos") %>% 
  filter(tot_spores > 100000000) %>% 
  filter(!lake_id %in% c("NN1", "Chapman", "Memorial", "Oglethorpe")) %>% 
  ggplot(., aes(x = lakeday, y = rel_spores/mean_spores, fill = species)) + 
  geom_bar(stat = "identity", position = "stack") + 
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) + guides(fill = "none") + labs(y = "# of infected individuals")

# lake_int %>% arrange(desc(tot_spores)) %>% mutate(lakeday=factor(lakeday, levels=unique(lakeday))) %>%
#   filter(tot_spores > 0 & tot_spores < 10000000) %>% 
#   filter(!lake_id %in% c("NN1", "Chapman", "Memorial", "Oglethorpe")) %>% 
#   group_by(lakeday) %>% slice_max(rel_spores) %>% 
#   ggplot(., aes(x = lakeday, y = AB*prev, color = species, fill = mean_spores)) + 
#   scale_fill_viridis_c() + 
#   geom_bar(stat = "identity", position = "stack") + 
#   theme(axis.text.x = element_text(angle = 70, hjust = 1)) + labs(y = "# of infected individuals")


p1 / p2






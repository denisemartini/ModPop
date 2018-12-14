tbl=read.table("~/Desktop/biallelic_common_snps.8.Q")
pops=read.table("~/Desktop/admixture/pops.txt")
pops2=read.table("~/Desktop/pops2.txt")
colnames(pops) <- c("Individual", "Population", "Island")
colnames(pops2) <- c("Individual", "Population", "Island")

full_tbl <- data.frame(pops, tbl)

library(reshape2)
m_tbl = melt(full_tbl, id.vars=c("Individual", "Population", "Island"), 
            variable.name="Ancestry", value.name="Fraction")

library(dplyr)

m_tbl %>% mutate(Population = factor(Population, levels = c('LittleBarrier', 'Pureora', 
                                                            'Zealandia', 'Kapiti', 
                                                            'Nelson', 'Westland', 
                                                            'Fiordland', 'Codfish'))) %>% arrange(Population) -> ord_tbl

library(ggplot2)

ggplot(ord_tbl, aes(x=Individual, y=Fraction, fill=Ancestry, order=Population)) +
  geom_bar(stat="identity", position="stack", width=1, colour="white") +
  facet_grid(. ~ Population, drop=TRUE, space="free", scales="free", switch = "x") +
  scale_fill_manual(values=brewer.pal(8, "RdYlGn")) +
  theme(legend.position = "None") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.spacing.x=grid:::unit(0.5, "lines")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank())


tbl2=read.table("~/Desktop/excluded_pops.2.Q")
full_tbl2 <- data.frame(pops2, tbl2)
m_tbl2 = melt(full_tbl2, id.vars=c("Individual", "Population", "Island"), 
             variable.name="Ancestry", value.name="Fraction")
m_tbl2 %>% mutate(Population = factor(Population, levels = c('LittleBarrier', 'Pureora', 
                                                            'Zealandia', 'Kapiti', 
                                                            'Nelson', 'Westland', 
                                                            'Fiordland', 'Codfish'))) %>% arrange(Population) -> ord_tbl2
ggplot(ord_tbl2, aes(x=Individual, y=Fraction, fill=Ancestry, order=Population)) +
  geom_bar(stat="identity", position="stack", width=1, colour="white") +
  facet_grid(. ~ Population, drop=TRUE, space="free", scales="free", switch = "x") +
  scale_fill_manual(values=brewer.pal(4, "RdYlGn")) +
  theme(legend.position = "None") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.spacing.x=grid:::unit(0.5, "lines")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank())


tbl3=read.table("~/Desktop/excluded_pops.3.Q")
full_tbl3 <- data.frame(pops2, tbl3)
m_tbl3 = melt(full_tbl3, id.vars=c("Individual", "Population", "Island"), 
              variable.name="Ancestry", value.name="Fraction")
m_tbl3 %>% mutate(Population = factor(Population, levels = c('LittleBarrier', 'Pureora', 
                                                             'Zealandia', 'Kapiti', 
                                                             'Nelson', 'Westland', 
                                                             'Fiordland', 'Codfish'))) %>% arrange(Population) -> ord_tbl3
ggplot(ord_tbl3, aes(x=Individual, y=Fraction, fill=Ancestry, order=Population)) +
  geom_bar(stat="identity", position="stack", width=1, colour="white") +
  facet_grid(. ~ Population, drop=TRUE, space="free", scales="free", switch = "x") +
  scale_fill_manual(values=brewer.pal(4, "RdYlGn")) +
  theme(legend.position = "None") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.spacing.x=grid:::unit(0.5, "lines")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank())

tbl4=read.table("~/Desktop/admixture/biallelic_common_snps.4.Q")
full_tbl4 <- data.frame(pops, tbl4)
m_tbl4 = melt(full_tbl4, id.vars=c("Individual", "Population", "Island"), 
              variable.name="Ancestry", value.name="Fraction")
m_tbl4 %>% mutate(Population = factor(Population, levels = c('LittleBarrier', 'Pureora', 
                                                             'Zealandia', 'Kapiti', 
                                                             'Nelson', 'Westland', 
                                                             'Fiordland', 'Codfish'))) %>% arrange(Population) -> ord_tbl4
ggplot(ord_tbl4, aes(x=Individual, y=Fraction, fill=Ancestry, order=Population)) +
  geom_bar(stat="identity", position="stack", width=1, colour="white") +
  facet_grid(. ~ Population, drop=TRUE, space="free", scales="free", switch = "x") +
  scale_fill_manual(values=brewer.pal(4, "RdYlGn")) +
  theme(legend.position = "None") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.spacing.x=grid:::unit(0.5, "lines")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.ticks.x=element_blank())

tbl5=read.table("~/Desktop/excluded_pops.5.Q")
full_tbl5 <- data.frame(pops2, tbl5)
m_tbl5 = melt(full_tbl5, id.vars=c("Individual", "Population", "Island"), 
              variable.name="Ancestry", value.name="Fraction")
m_tbl5 %>% mutate(Population = factor(Population, levels = c('LittleBarrier', 'Pureora', 
                                                             'Zealandia', 'Kapiti', 
                                                             'Nelson', 'Westland', 
                                                             'Fiordland', 'Codfish'))) %>% arrange(Population) -> ord_tbl5
ggplot(ord_tbl5, aes(x=Individual, y=Fraction, fill=Ancestry, order=Population)) +
  geom_bar(stat="identity", position="stack", width=1, colour="white") +
  facet_grid(. ~ Population, drop=TRUE, space="free", scales="free", switch = "x") +
  scale_fill_manual(values=brewer.pal(6, "RdYlGn")) +
  theme(legend.position = "None") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.spacing.x=grid:::unit(0.5, "lines")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank())


tbl6=read.table("~/Desktop/excluded_pops.6.Q")
full_tbl6 <- data.frame(pops2, tbl6)
m_tbl6 = melt(full_tbl6, id.vars=c("Individual", "Population", "Island"), 
              variable.name="Ancestry", value.name="Fraction")
m_tbl6 %>% mutate(Population = factor(Population, levels = c('LittleBarrier', 'Pureora', 
                                                             'Zealandia', 'Kapiti', 
                                                             'Nelson', 'Westland', 
                                                             'Fiordland', 'Codfish'))) %>% arrange(Population) -> ord_tbl6
ggplot(ord_tbl6, aes(x=Individual, y=Fraction, fill=Ancestry, order=Population)) +
  geom_bar(stat="identity", position="stack", width=1, colour="white") +
  facet_grid(. ~ Population, drop=TRUE, space="free", scales="free", switch = "x") +
  scale_fill_manual(values=brewer.pal(6, "RdYlGn")) +
  theme(legend.position = "None") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.spacing.x=grid:::unit(0.5, "lines")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank())


tbl7=read.table("~/Desktop/biallelic_common_snps.7.Q")
full_tbl7 <- data.frame(pops, tbl7)
m_tbl7 = melt(full_tbl7, id.vars=c("Individual", "Population", "Island"), 
              variable.name="Ancestry", value.name="Fraction")
m_tbl7 %>% mutate(Population = factor(Population, levels = c('LittleBarrier', 'Pureora', 
                                                             'Zealandia', 'Kapiti', 
                                                             'Nelson', 'Westland', 
                                                             'Fiordland', 'Codfish'))) %>% arrange(Population) -> ord_tbl7
ggplot(ord_tbl7, aes(x=Individual, y=Fraction, fill=Ancestry, order=Population)) +
  geom_bar(stat="identity", position="stack", width=1, colour="white") +
  facet_grid(. ~ Population, drop=TRUE, space="free", scales="free", switch = "x") +
  scale_fill_manual(values=brewer.pal(8, "RdYlGn")) +
  theme(legend.position = "None") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.spacing.x=grid:::unit(0.5, "lines")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank())


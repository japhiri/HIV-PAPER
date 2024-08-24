#installing and loading required packages----
pacman::p_load(char = c('lubridate',"gtsummary", 'tidyverse', "dplyr", "here", "rio", "scales", "boot", 
                        "magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl", 
                        "ggsignif", "ggpubr", "ggeasy", "cowplot","ggExtra", "PupillometryR","hrbrthemes", "ggstance",
                        "survival","survminer","sysfonts","showtext","nlme"))
# Preparing input samples ----
# read in Micro.csv and MPO NETs .csv data file
Micro <- read_csv("Data/Main_Files_Thesis/Micro.csv")
NETSMPO <- read_csv("Data/Main_Files_Thesis/MPO_NETs.csv")
Innate_Like_2024 <- read_csv("Data/Main_Files_Thesis/Innate-Like_2024.csv") %>%
  select(-Sample) %>%
  filter(`Sample Type`=="Nasal") %>%
  select(-`Sample Type`,-`Sample Quality`)
# List all .csv files in the folder
Maincsv_files <- list.files(path = "Data/Main_Files_thesis", pattern = "_Main.csv", full.names = TRUE)
Maincsv_files
# read the CSV files
Maindata_list <- lapply(Maincsv_files, read_csv)
# Combine all the main data
Maindata <- do.call(rbind, Maindata_list) %>%
  select(-Sample) %>%
  mutate(NER = `Neutrophil count`/`Epithelial count`,
         MER = `Monocyte count`/`Epithelial count`,
         LER = `CD45+ count`/`Epithelial count`) %>%
  mutate(NCD45ER = `Neutrophil count`/`LER`,
         MCD45ER = `Monocyte count`/`LER`,
         CD11bER = `Neutrophil CD11b++`/`Epithelial count`) %>%
  filter(`Sample Type`=="Nasal",
         `Sample Quality`=="Good")
# List all .csv files in the folder containing NeutrophilCD10CD63
NeutCD10CD63csv_files <- list.files(path = "Data/Main_Files_thesis", pattern = "_NeutrophilCD10CD63.csv", full.names = TRUE)
NeutCD10CD63csv_files
# read the CSV files
NeutCD10CD63data_list <- lapply(NeutCD10CD63csv_files, read_csv)
# Combine all the NeutCD10CD63 data
NeutCD10CD63 <- do.call(rbind, NeutCD10CD63data_list) %>%
  select(-Sample,-`Sample Type`,-`Sample Quality`,-`CD3 Staining`,-`CD11b Staining`,-`Panel`)
# List all .csv files in the folder containing MonocyteCD10CD63
MonoCD10CD63csv_files <- list.files(path = "Data/Main_Files_thesis", pattern = "_MonocyteCD10CD63.csv", full.names = TRUE)
MonoCD10CD63csv_files
# read the CSV files
MonoCD10CD63data_list <- lapply(MonoCD10CD63csv_files, read_csv)
# Combine all the MonoCD10CD63 data
MonoCD10CD63 <- do.call(rbind, MonoCD10CD63data_list) %>%
  select(-Sample,-`Sample Type`,-`Sample Quality`,-`CD3 Staining`,-`CD11b Staining`,`Panel`)
# Merge all the files
Main_Myeloid <- Maindata %>%
  merge(NeutCD10CD63, by=c("LAB ID"),all=T) %>%
  merge(MonoCD10CD63, by=c("LAB ID"),all=T) %>%
  merge(NETSMPO, by=c("LAB ID"),all=T) %>%
  merge(Innate_Like_2024, by=c("LAB ID"),all = T) %>%
  merge(Micro,by=c("LAB ID"),all=F)

# Save the Main_results  
write.csv(Main_Myeloid,'Results/Main_Myeloid.csv',row.names = F)
Main_Myeloid <- read.csv('Results/Main_Myeloid.csv')

# FIGURE 1: NASAL IMMUNE CELL LANSCAPE DURING HIV INFECTION----
# Immune cell abundance during HIV infection
Immune <- Main_Myeloid %>%
  filter(`HIV Status`=="HIV-",
         `Epithelial count`>100,
         `CD45+ count`>200,
         `Visit`=='Week 1',
         `Sample Type`=="Nasal",
         `CD3 Staining`=='Stained'
  ) %>%
  pivot_longer(cols = c("Neutrophil proportion","Tcell proportion","Monocyte proportion"),
               names_to = "Immune cells",
               values_to = "Proportion of immune cells") %>%
  group_by(`Sample Type`,`HIV Status`) %>%
  ggplot(aes(x=factor(`Immune cells`,
                      levels=c("Monocyte proportion","Tcell proportion","Neutrophil proportion")),
             y=`Proportion of immune cells`,
             fill=factor(`Immune cells`,
                         levels=c("Monocyte proportion","Tcell proportion","Neutrophil proportion"))))+
  geom_boxplot(width=0.5,notch = T)+
  #geom_jitter(position = position_jitter(width = 0.1),size=4,aes(shape = `HIV Status`, color=`HIV Status`)) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.2)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               #shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p.format",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_fill_manual(values = c('#009193','#941100','#005493'))+
  labs(x='',
       y='Proportion of Immune cells',
       fill='Immune cells')+
  scale_y_continuous(limits = c(0,120), breaks = c(20,40,60,80,100,120))+
  scale_x_discrete(labels=c("Monocyte proportion"="Monocytes",
                            "Tcell proportion"="T cells",
                            "Neutrophil proportion"="Neutrophils"))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30,face = "bold"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

Immune

# Neutrophil to epithelial cell ratio during HIV infection
NER <- Main_Myeloid %>%
  filter(`Sample Quality`=='Good',
         `Epithelial count`>100,
         `CD45+ count`>200,
         `Sample Type`=='Nasal',
         Visit=="Week 1",
         #`Carriage Status`=="Carriage Negative"
         ) %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`NER`),
             fill=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = T)+
  #geom_jitter(position = position_jitter(width = 0.1),size=4,aes(shape = `HIV Status`, color=`HIV Status`)) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.2)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               #shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p.format",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  labs(x='',
       y='Neutrophil Epithelial Ratio',
       fill='HIV Status')+
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                limits = c(10^-2.5,10^2))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="ART<3m",
                            "PLHIV ART >1y"="ART>1y"))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30,face = "bold"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

NER

# Monocyte to epithelial cell ratio during HIV infection 
MER <- Main_Myeloid %>%
  filter(`Sample Quality`=='Good',
         `Epithelial count`>100,
         `CD45+ count`>200,
         `Sample Type`=='Nasal',
         Visit=="Week 1",
         #`Carriage Status`=="Carriage Negative"
  ) %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`MER`),
             fill=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = T)+
  #geom_jitter(position = position_jitter(width = 0.1),size=4,aes(shape = `HIV Status`, color=`HIV Status`)) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.2)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               #shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p.format",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  labs(x='',
       y='Monocyte to Epithelial Ratio',
       fill='HIV Status')+
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                limits = c(10^-4,10^1))+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="ART<3m",
                            "PLHIV ART >1y"="ART>1y"))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30,face = "bold"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

MER

# CD4+ T cells during HIV infection
CD4 <- Main_Myeloid %>%
  filter(`CD4+`<80) %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD4+`),
             fill=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = T)+
  #geom_jitter(position = position_jitter(width = 0.1),size=0.5) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.2)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               #shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p.format",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  labs(x='',
       y='Proportion of CD4+ T cells',
       fill='HIV Status')+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="ART<3m",
                            "PLHIV ART >1y"="ART>1y"))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30,face = "bold"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

CD4

# CD8+ T cells during HIV infection
CD8 <- Main_Myeloid %>%
  filter(`CD8+`>20) %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD8+`),
             fill=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = T)+
  #geom_jitter(position = position_jitter(width = 0.1),size=0.5) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.2)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               #shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p.format",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  labs(x='',
       y='Proportion of CD8+ T cells',
       fill='HIV Status')+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="ART<3m",
                            "PLHIV ART >1y"="ART>1y"))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30,face = "bold"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

CD8

# CD8+ T cells during HIV infection
`DN-` <- Main_Myeloid %>%
  #filter(`CD8+`>20) %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`DN-`),
             fill=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = T)+
  #geom_jitter(position = position_jitter(width = 0.1),size=0.5) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.2)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               #shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p.format",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  labs(x='',
       y='Proportion of DN- T cells',
       fill='HIV Status')+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="ART<3m",
                            "PLHIV ART >1y"="ART>1y"))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30,face = "bold"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

`DN-`

# CD3+ Mait cells during HIV infection
 `CD3+Mait`<- Main_Myeloid %>%
  #filter(`CD8+`>20) %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD3+Mait`),
             fill=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = T)+
  #geom_jitter(position = position_jitter(width = 0.1),size=0.5) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.2)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               #shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p.format",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  labs(x='',
       y='Proportion of CD3+ MAIT cells',
       fill='HIV Status')+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="ART<3m",
                            "PLHIV ART >1y"="ART>1y"))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30,face = "bold"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

`CD3+Mait`

# CD3+ Gamma Delta T cells during HIV infection
`CD3+yD`<- Main_Myeloid %>%
  #filter(`CD8+`>20) %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD3+TCRgd+`),
             fill=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = T)+
  #geom_jitter(position = position_jitter(width = 0.1),size=0.5) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.2)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               #shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p.format",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  labs(x='',
       y='Proportion of CD3+ Gamma delta T cells',
       fill='HIV Status')+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="ART<3m",
                            "PLHIV ART >1y"="ART>1y"))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30,face = "bold"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

`CD3+yD`

# CD3+ Natural Killer T cells during HIV infection
`CD3+NK`<- Main_Myeloid %>%
  #filter(`CD8+`>20) %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD3+CD56+`),
             fill=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = T)+
  #geom_jitter(position = position_jitter(width = 0.1),size=0.5) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.2)+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               #shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p.format",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  labs(x='',
       y='Proportion of CD3+ Natural killer T cells',
       fill='HIV Status')+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="ART<3m",
                            "PLHIV ART >1y"="ART>1y"))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30,face = "bold"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

`CD3+NK`
# SAVE Figure 1
Fig1 <- (Immune|NER|MER|CD4|CD8|`DN-`|`CD3+Mait`|`CD3+yD`|`CD3+NK`)+
  plot_layout(ncol = 4,nrow = 3,widths = c(1,1,1,1))+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30,face = "bold"))

# Save Fig1
ggsave(filename = here("Figures/Figure 1","Fig1.pdf"),
       plot = Fig1,
       width = 30,height = 32,units = "in",dpi = 700)
# Save Fig1
ggsave(filename = here("Figures/Figure 1","Fig1.png"),
       plot = Fig1,
       width = 30,height = 32,units = "in",dpi = 700)

# Figure 2: Nasal Inflammation in PLHIV characterised by high abundance of activate neutrophils ----



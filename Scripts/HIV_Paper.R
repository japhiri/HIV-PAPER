#installing and loading required packages----
pacman::p_load(char = c('lubridate',"gtsummary", 'tidyverse', "dplyr", "here", "rio", "scales", "boot", "corrplot",
                        "magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl","tidyheatmaps","Hmisc",
                        "ggsignif", "ggpubr", "ggeasy", "cowplot","ggExtra", "PupillometryR","hrbrthemes", "ggstance",
                        "survival","survminer","sysfonts","showtext","nlme","pheatmap","ComplexHeatmap","glue","reshape2"))
# Preparing input samples ----
# read in Micro.csv and MPO NETs .csv data file
Clinical_Data <- read_csv("Data/Main_Files_Thesis/Clinical_Data.csv")
Lab_Data <- read_csv("Data/Main_Files_Thesis/Lab_Data.csv")
Micro <- read_csv("Data/Main_Files_Thesis/Micro.csv")
NETSMPO <- read_csv("Data/Main_Files_Thesis/MPO_NETs.csv")
Cytokines_averages <- read_csv("Data/Main_Files_Thesis/Cytokines_averages.csv") 
Innate_Like_2024 <- read_csv("Data/Main_Files_Thesis/Innate-Like_2024.csv") %>%
  dplyr::select(-Sample) %>%
  dplyr::filter(`Sample Type`=="Nasal") %>%
  dplyr::select(-`Sample Type`,-`Sample Quality`)
# List all .csv files in the folder
CD3_files <- list.files(path = "Data/Main_Files_thesis", pattern = "_HIVPaper_CD3.csv", full.names = TRUE)
CD3_files
# read the CSV files
CD3_files <- lapply(CD3_files, read_csv)
# Combine all the main data
CD3_files <- do.call(rbind, CD3_files) %>%
  dplyr::select(-Sample)
# List all .csv files in the folder
Maincsv_files <- list.files(path = "Data/Main_Files_thesis", pattern = "_Main.csv", full.names = TRUE)
Maincsv_files
# read the CSV files
Maindata_list <- lapply(Maincsv_files, read_csv)
# Combine all the main data
Maindata <- do.call(rbind, Maindata_list) %>%
  dplyr::select(-Sample) %>%
  dplyr::mutate(NER = `Neutrophil count`/`Epithelial count`,
         MER = `Monocyte count`/`Epithelial count`,
         LER = `CD45+ count`/`Epithelial count`) %>%
  dplyr::mutate(NCD45ER = `Neutrophil count`/`LER`,
         MCD45ER = `Monocyte count`/`LER`,
         CD11bER = `Neutrophil CD11b++`/`Epithelial count`) %>%
  dplyr::filter(`Sample Type`=="Nasal",
         `Sample Quality`=="Good") %>%
  dplyr::select(-`Sample Quality`,-`Sample Type`,-Panel)
Maindata$`LAB ID` <- gsub("CHU106","CUH106", Maindata$`LAB ID`)
Maindata$`LAB ID` <- gsub("CUF12BN1","CUF12B", Maindata$`LAB ID`)
# List all .csv files in the folder containing NeutrophilCD10CD63
NeutCD10CD63csv_files <- list.files(path = "Data/Main_Files_thesis", pattern = "_NeutrophilCD10CD63.csv", full.names = TRUE)
NeutCD10CD63csv_files
# read the CSV files
NeutCD10CD63data_list <- lapply(NeutCD10CD63csv_files, read_csv)
# Combine all the NeutCD10CD63 data
NeutCD10CD63 <- do.call(rbind, NeutCD10CD63data_list) %>%
  dplyr::select(-Sample,-`Sample Type`,-`Sample Quality`,-`CD3 Staining`,-`CD11b Staining`,-`Panel`)
# List all .csv files in the folder containing MonocyteCD10CD63
MonoCD10CD63csv_files <- list.files(path = "Data/Main_Files_thesis", pattern = "_MonocyteCD10CD63.csv", full.names = TRUE)
MonoCD10CD63csv_files
# read the CSV files
MonoCD10CD63data_list <- lapply(MonoCD10CD63csv_files, read_csv)
# Combine all the MonoCD10CD63 data
MonoCD10CD63 <- do.call(rbind, MonoCD10CD63data_list) %>%
  dplyr::select(-Sample,-`Sample Type`,-`Sample Quality`,-`CD3 Staining`,-`CD11b Staining`,-`Panel`)
# Merge all the files
Masterfile <- Maindata %>%
  merge(NeutCD10CD63, by=c("LAB ID"),all=T) %>%
  merge(MonoCD10CD63, by=c("LAB ID"),all=T) %>%
  merge(NETSMPO, by=c("LAB ID"),all=T) %>%
  merge(Cytokines_averages, by=c("LAB ID"),all=T) %>%
  merge(Innate_Like_2024, by=c("LAB ID"),all = T) %>%
  merge(Lab_Data, by=c("LAB ID"),all = T) %>%
  merge(Clinical_Data, by=c("LAB ID"),all = T) %>%
  merge(Micro,by=c("LAB ID"),all=T) %>%
  merge(CD3_files, by=c("LAB ID"),all = T)
Masterfile <- Masterfile[!duplicated(Masterfile$`LAB ID`), ]

# Save the Main_results  
write.csv(Masterfile,'Results/Masterfile.csv',row.names = F)
#Masterfile <- read_csv('Results/Masterfile.csv')

HIV_status_colors <- c('#A9A9A9','#941100','#005493')
Masterfile %>%
  dplyr::filter(`Epithelial count` > 100,
                Visit=="Week 1") %>%
  ggplot(aes(`HIV Status`,`Epithelial count`))+
  geom_boxplot()+
  geom_pwc()

# FIGURE 1: NASAL IMMUNE CELL LANSCAPE DURING HIV INFECTION----
# Immune cell abundance during HIV infection
group_counts <- Masterfile %>%
  dplyr::filter(`HIV Status` == "HIV-",
         `Epithelial count` > 100,
         `CD45+ count` > 200,
         `Visit` == 'Week 1',
         `CD3 Staining` == 'Stained') %>%
 dplyr::mutate(TER = `T cell count` / `Epithelial count`) %>%
  dplyr::select(NER,MER,TER) %>%
  pivot_longer(cols = c("NER", "TER", "MER"),
               names_to = "Immune cells",
               values_to = "Immune cell to epithelial cell ratio") %>%
  dplyr::group_by(`Immune cells`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")


Immune <- Masterfile %>%
  filter(`HIV Status`=="HIV-",
         `Epithelial count`>100,
         `CD45+ count`>200,
         `Visit`=='Week 1',
         `CD3 Staining`=='Stained'
  ) %>%
  mutate(TER=`T cell count`/`Epithelial count`) %>%
  pivot_longer(cols = c("NER","TER","MER"),
               names_to = "Immune cells",
               values_to = "Immune cell to epithelial cell ratio") %>%
  ggplot(aes(x=factor(`Immune cells`,
                      levels=c("NER","MER","TER")),
             y=`Immune cell to epithelial cell ratio`,
             shape=factor(`Immune cells`,
                         levels=c("NER","MER","TER")),
             color=factor(`Immune cells`,
                          levels=c("NER","MER","TER")),
             fill=factor(`Immune cells`,
                          levels=c("NER","MER","TER"))))+
  #geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=10,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 'wilcox.test',
           label = "p = {p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F,
           vjust = -0.2)+
  scale_color_manual(values = c('#000','#000','#000'))+
  scale_fill_manual(values = c('#A9A9A9','#A9A9A9','#A9A9A9'))+
  scale_shape_manual(values = c(21, 24, 22))+
  labs(x='',
       y="Immune cell to epithelial cell ratio",
       fill='Immune cells')+
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                limits = c(10^-3,10^1.5), breaks=c(10^-3,10^-2,10^-1,10^0,10^1,10^2))+
  scale_x_discrete(labels=c("MER"="Monocytes",
                            "TER"="T cells",
                            "NER"="Neutrophils"))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 20, face = "bold"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

Immune

ggsave(Immune,filename="Figures/Figure 1/immune.png",
       width = 7,height = 9,dpi = 1080)


# Neutrophil to epithelial cell ratio during HIV infection
HIV_status_colors <- c('#A9A9A9','#941100','#005493')
group_counts <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         `Visit`=="Week 1") %>%
  dplyr::select(NER,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")

NER <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         `Visit`=="Week 1") %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`NER`),
             color=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p = {p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

NER

ggsave(NER,filename="Figures/Figure 1/NER.png",
       width = 7,height = 9,dpi = 1080)

# Monocyte to epithelial cell ratio during HIV infection 
group_counts <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         `Visit`=="Week 1") %>%
  select(MER,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")

MER <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         Visit=="Week 1") %>%
  group_by(`HIV Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`MER`),
             color=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p = {p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y='Monocyte to Epithelial Ratio',
       fill='HIV Status')+
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                limits = c(10^-4,10^1), breaks=c(10^-4,10^-3,10^-2,10^-1,10^0,10^1))+
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

MER

ggsave(MER,filename="Figures/Figure 1/MER.png",
       width = 7,height = 9,dpi = 1080)


# T cell to epithelial cell ratio during HIV infection 

# Monocyte to epithelial cell ratio during HIV infection 
group_counts <- Masterfile %>%
  dplyr::filter(`Epithelial count`>100,
                `CD45+ count`>200,
                Visit=="Week 1",
                `CD3 Staining`=="Stained") %>%
  dplyr::mutate(TER=`T cell count`/`Epithelial count`) %>%
  dplyr::select(TER,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
group_counts

CD3 <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         Visit=="Week 1",
         `CD3 Staining`=="Stained"
         #`Carriage Status`=="Carriage Negative"
  ) %>%
  mutate(TER=`T cell count`/`Epithelial count`) %>%
  filter(TER<1) %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`TER`),
             color=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't_test',
           label = "p = {p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y='T cell to Epithelial Ratio',
       fill='HIV Status')+
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1))+
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

CD3

ggsave(CD3,filename="Figures/Figure 1/CD3.png",
       width = 7,height = 9,dpi = 1080)

# SAVE Figure 1
Fig1c_f <- (Immune|NER|MER|CD3)+
  plot_layout(ncol = 4,nrow = 3,widths = c(1,1,1,1))+
  plot_annotation(tag_levels = list(c("C","D","E","F"))) &
  theme(plot.tag = element_text(size = 40,face = "bold"))

# Save Fig1
ggsave(filename = here("Figures/Figure 1","Fig1c_f.pdf"),
       plot = Fig1c_f,
       width = 30,height = 32,units = "in",dpi = 700)

# Save Fig1
ggsave(filename = here("Figures/Figure 1","Fig1c_f.png"),
       plot = Fig1c_f,
       width = 30,height = 32,units = "in",dpi = 700)


# CD4+ T cells during HIV infection
group_counts <- Masterfile %>%
  dplyr::filter(`CD4+`<80,
                `HIV Status`!='NA',
                Visit=="Week 1") %>%
  dplyr::select(`CD4+`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
group_counts

CD4 <- Masterfile %>%
  filter(`CD4+`<80,
         `HIV Status`!='NA',
         Visit=="Week 1") %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD4+`),
             color=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color="black")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p={p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Proportion of CD4"^"+"~"T cells"),
       fill='HIV Status')+
  scale_y_continuous(limits = c(0,60),breaks = c(0,10,20,30,40,50,60))+
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

CD4

ggsave(CD4,filename="Figures/Figure 1/CD4.png",
       width = 7,height = 9,dpi = 1080)

# CD8+ T cells during HIV infection
group_counts <- Masterfile %>%
  dplyr::filter(`CD8+`>20,
                Visit=="Week 1",
                `HIV Status`!="NA") %>%
  dplyr::select(`CD8+`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
group_counts

CD8 <- Masterfile %>%
  filter(`CD8+`>20,
         Visit=="Week 1",
         `HIV Status`!="NA") %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD8+`),
             color=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p={p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Proportion of CD8"^"+"~"T cells"),
       fill='HIV Status')+
  scale_y_continuous(limits = c(20,120),breaks = c(20,40,60,80,100,120))+
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

CD8


ggsave(CD8,filename="Figures/Figure 1/CD8.png",
       width = 7,height = 9,dpi = 1080)


# CD8+ T cells during HIV infection
group_counts <- Masterfile %>%
  dplyr::filter(`CD8+`>20,
                Visit=="Week 1",
                `HIV Status`!="NA") %>%
  dplyr::select(`DN-`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
group_counts

`DN-` <- Masterfile %>%
  filter(`CD8+`>20,
         Visit=="Week 1",
         `HIV Status`!="NA") %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`DN-`),
             color=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color="black")+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p={p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Proportion of DN"^"-"~"T cells"),
       fill='HIV Status')+
  scale_y_continuous(limits = c(0,90),breaks = c(0,20,40,60,80))+
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

`DN-`

ggsave(`DN-`,filename="Figures/Figure 1/DN-.png",
       width = 7,height = 9,dpi = 1080)

# Save figure
Fig1h_j <- (CD4|CD8|`DN-`)+
  plot_layout(ncol = 4,nrow = 3,widths = c(1,1,1,1))+
  plot_annotation(tag_levels = list(c("H","I","J"))) &
  theme(plot.tag = element_text(size = 40,face = "bold"))

# Save Fig1
ggsave(filename = here("Figures/Figure 1","Fig1h_j.pdf"),
       plot = Fig1h_j,
       width = 30,height = 32,units = "in",dpi = 700)
# Save Fig1
ggsave(filename = here("Figures/Figure 1","Fig1h_j.png"),
       plot = Fig1h_j,
       width = 30,height = 32,units = "in",dpi = 700)


# CD3+ Mait cells during HIV infection
group_counts <- Masterfile %>%
  dplyr::filter(Visit=="Week 1",
                `HIV Status`!="NA") %>%
  dplyr::select(`CD3+Mait`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  na.omit() %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
group_counts

 `CD3+Mait`<- Masterfile %>%
  filter(`HIV Status`!='NA',
         Visit=="Week 1") %>%
  group_by(`HIV Status`,`Carriage Status`) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD3+Mait`),
             color=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
   geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
   geom_jitter(position = position_jitter(width = 0.3),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15, color='black')+
   stat_summary(geom = "point",
                fun = median,
                color='black',
                size=12,
                shape=95,
                position = position_dodge(width = 0.75))+
   geom_pwc(method = 't.test',
            label = "p={p.format}",
            label.size = 8,
            tip.length = 0.01,
            hide.ns = F)+
   scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Proportion of CD3"^"+"~"MAIT cells"),
       fill='HIV Status')+
   scale_y_continuous(limits = c(0,40),breaks = c(0,10,20,30,40))+
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
         axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
         #axis.text.x = element_blank(),
         axis.text.y = element_text(size = 15,face = "bold"),
         axis.title = element_text(size = 30))

`CD3+Mait`

ggsave(`CD3+Mait`,filename="Figures/Figure 1/CD3+Mait.png",
       width = 7,height = 9,dpi = 1080)

# CD3+ Gamma Delta T cells during HIV infection
group_counts <- Masterfile %>%
  dplyr::filter(Visit=="Week 1",
                `HIV Status`!="NA") %>%
  dplyr::select(`CD3+TCRgd+`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  na.omit() %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
group_counts

`CD3+yD`<- Masterfile %>%
  filter(`HIV Status`!='NA',
         Visit=="Week 1") %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD3+TCRgd+`),
             color=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15, color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p={p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Proportion of TCR"~gamma~delta~""^"+"~"T cells"),
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

`CD3+yD`

ggsave(`CD3+yD`,filename="Figures/Figure 1/CD3+yD.png",
       width = 7,height = 9,dpi = 1080)

# CD3+ Natural Killer T cells during HIV infection
group_counts <- Masterfile %>%
  dplyr::filter(Visit=="Week 1",
                `HIV Status`!="NA") %>%
  dplyr::select(`CD3+CD56+`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  na.omit() %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
group_counts

`CD3+NK`<- Masterfile %>%
  filter(`HIV Status`!='NA',
         Visit=="Week 1") %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`CD3+CD56+`),
             color=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p={p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Proportion of CD3"^"+"~"Natural Killer T cells"),
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

`CD3+NK`

ggsave(`CD3+NK`,filename="Figures/Figure 1/CD3+NK.png",
       width = 7,height = 9,dpi = 1080)

# Save figure
Fig1L_N <- (`CD3+Mait`|`CD3+yD`|`CD3+NK`)+
  plot_layout(ncol = 4,nrow = 3,widths = c(1,1,1,1))+
  plot_annotation(tag_levels = list(c("l","M","N"))) &
  theme(plot.tag = element_text(size = 40,face = "bold"))

# Save Fig1
ggsave(filename = here("Figures/Figure 1","Fig1L_N.pdf"),
       plot = Fig1L_N,
       width = 30,height = 32,units = "in",dpi = 700)
# Save Fig1
ggsave(filename = here("Figures/Figure 1","Fig1L_N.png"),
       plot = Fig1L_N,
       width = 30,height = 32,units = "in",dpi = 700)

# SAVE Figure 1
Fig1 <- (Immune|NER|CD3|MER|CD4|CD8|`DN-`|`CD3+Mait`|`CD3+yD`|`CD3+NK`)+
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
# Relationship between Neutrophil CD10+CD11b+CD62L-CD63- and neutrophil CD11b++
group_counts <- Masterfile %>%
  dplyr::filter(Visit=="Week 1",
                `HIV Status`!="NA") %>%
  dplyr::select(`CD3+CD56+`,`HIV Status`) %>%
  dplyr::group_by(`HIV Status`) %>%
  na.omit() %>%
  dplyr::summarize(n = dplyr::n(), .groups = "drop")
group_counts
`CD11b++Correlation` <- Masterfile %>%
  filter(`CD11b Staining`=="Stained",
         #`Visit`=="Week 1",
         `HIV Status`=="HIV-",
         `Neutrophil CD10+CD11b+CD62L-CD63-`>20,
         `Neutrophil CD11b++`>20 &
           `Neutrophil CD11b++`<90
         ) %>%
  ggplot(aes(`Neutrophil CD10+CD11b+CD62L-CD63-`,`Neutrophil CD11b++`))+
  geom_point(size=7)+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman",size=12)+
  labs(x= expression(atop("Proportion of CD10"^"+"~"CD11b"^"++"~"CD63"^"-"~"CD62L"^"-", 
      "Neutrophils"
    )
  ),
       y=expression("Proportion of CD11b"^"++"~" Neutrophils"))+
  scale_y_continuous(limits = c(50,75),breaks = c(50,60,70))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30,face = 'bold'),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size =30,face = 'bold'))
`CD11b++Correlation`
ggsave(`CD11b++Correlation`,filename="Figures/Figure 2/CD11b++Correlation.png",
       width = 10,height = 12,dpi = 1080)

# Save CD11b++ correlation
Fig2c <- (`CD11b++Correlation`)+
  plot_layout(ncol = 3,nrow = 1,widths = c(1,1,1,1))+
  #plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30,face = "bold"))

# Save Fig2
ggsave(filename = here("Figures/Figure 2","Fig2c.pdf"),
       plot = Fig2c,
       width = 30,height = 13,units = "in",dpi = 700)
# Save Fig2
ggsave(filename = here("Figures/Figure 2","Fig2c.png"),
       plot = Fig2c,
       width = 30,height = 13,units = "in",dpi = 700)

# Proportion of CD11b++ neutrophils in nasal lining fluid during HIV infection
`NeutCD11b++` <- Masterfile %>%
  filter(`CD11b Staining`=="Stained",
         `HIV Status`!="NA",
         `Visit`=="Week 1",
         #`Neutrophil CD11b++`<60,
         `Neutrophil CD11b++`>5) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('PLHIV ART <3m', 'HIV-','PLHIV ART >1y')),
             y=as.numeric(`Neutrophil CD11b++`),
             color=factor(`HIV Status`,
                         levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p={p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#941100','#A9A9A9','#005493'))+
  labs(x='',
       y=expression("Proportion of CD11b"^"++"~" Neutrophils"),
       fill='HIV Status')+
  scale_x_discrete(labels=c("HIV-"="ART<3m",
                            "PLHIV ART <3m"="HIV-",
                            "PLHIV ART >1y"="ART>1y"))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

`NeutCD11b++`

ggsave(`NeutCD11b++`,filename="Figures/Figure 2/NeutCD11b++.png",
       width = 7,height = 9,dpi = 1080)


# Proportion of CD62L+ neutrophils in nasal lining fluid during HIV infection
`NeutCD62L+` <- Masterfile %>%
  filter(`CD11b Staining`=="Stained",
         `Visit`=="Week 1",
         `Neutrophil CD11b++`<70) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`Neutrophil CD62L proportion`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p={p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Proportion of CD62L"^"+"~" Neutrophils"),,
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
`NeutCD62L+`

ggsave(`NeutCD62L+`,filename="Figures/Figure 2/NeutCD62L+.png",
       width = 7,height = 9,dpi = 1080)

# Median Flourescence intensity (MFI) of CD11b during HIV infection
NeutCD11bMFI <- Masterfile %>%
  filter(`CD11b Staining`=="Stained",
         `HIV Status`!="NA",
         Visit=="Week 1",
         `Neutrophil CD11b++`<70) %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`Neutrophil CD11b+ MFI`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p={p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y='Neutrophil CD11b Median Flourescence Intensity',
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
NeutCD11bMFI

ggsave(NeutCD11bMFI,filename="Figures/Figure 2/NeutCD11bMFI.png",
       width = 7,height = 9,dpi = 1080)

# Median fluorescence intensity (MFI) of CD62L neutrohils during HIV infection
NeutCD62LMFI <- Masterfile %>%
  filter(`Visit`=="Week 1",
         `Neutrophil CD62L+ MFI`<5000,
         `HIV Status`!="NA") %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`Neutrophil CD62L+ MFI`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.3),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=12,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p={p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y='Neutrophil CD62L Median \n Flourescence Intensity',
       fill='HIV Status')+
  scale_x_discrete(labels=c("HIV-"="HIV-",
                            "PLHIV ART <3m"="ART<3m",
                            "PLHIV ART >1y"="ART>1y"))+
  #scale_y_continuous(labels = label_scientific())+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
NeutCD62LMFI

ggsave(NeutCD62LMFI,filename="Figures/Figure 2/NeutCD62LMFI.png",
       width = 7,height = 9,dpi = 1080)

# Relationship between the proportion of CD11b+ neutrophils and the proportion of CD66b++ neutrophils
Masterfile %>%
  filter(`CD11b Staining`=="Stained",
         `Visit`=="Week 1",
         #`HIV Status`=="HIV-"
         ) %>%
  ggplot(aes(`Neutrophil CD11b proportion`,`Neutrophil CD66b++`))+
  geom_point(size=7)+
  geom_smooth(method = "lm")+
  stat_cor(method = 'spearman')+
  labs(x='Proportion of CD11b+ Neutrophils',
       y='Proportion of CD66b++ Neutrophils',
       fill='HIV Status')+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25),
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))


Masterfile %>%
  filter(`CD11b Staining`=="Stained",
         `Visit`=="Week 1",
         `HIV Status`=="HIV-") %>%
  ggplot(aes(x=as.numeric(`Neutrophil CD11b proportion`),
             y=as.numeric(`Neutrophil CD62L proportion`)))+
  geom_point(size=7)+
  geom_smooth(method = "lm")+
  stat_cor(method = 'spearman')+
  labs(x='Proportion of CD11b+ Neutrophils',
       y='Proportion of CD62L Neutrophils',
       fill='HIV Status')+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))



# Save CD11b++ correlation
Fig2ef <- (`NeutCD11b++`|NeutCD62LMFI)+
  plot_layout(ncol = 4,nrow = 3,widths = c(1,1,1,1))+
  #plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30,face = "bold"))

# Save Fig2
ggsave(filename = here("Figures/Figure 2","Fig2ef.pdf"),
       plot = Fig2ef,
       width = 30,height = 32,units = "in",dpi = 700)
# Save Fig2
ggsave(filename = here("Figures/Figure 2","Fig2ef.png"),
       plot = Fig2ef,
       width = 30,height = 32,units = "in",dpi = 700)

# Concentration of myeloperoxidase in nasal lining fluid during HIV infection
MPO <- Masterfile %>%
  filter(`HIV Status`!="NA") %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`Myeloperoxidase (pg/mL)`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.1),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15, color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p.format",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y='Concentration of Myeloperoxidase (pg/ML)',
       fill='HIV Status')+
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                limits = c(10^4,10^6.5))+
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

MPO

ggsave(MPO,filename="Figures/Figure 1/MPO.png",
       width = 7,height = 9,dpi = 1080)

# Concentration of Neutrophil extracellular traps (NETs) in nasal lining fluid during HIV infection
NETs <- Masterfile %>%
  filter(`HIV Status`!="NA") %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`NETs Concentration`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F,outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.1),size=7) +
  geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15, color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p.format",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y='Concentration of Neutrophil \n extracellular traps (pg/ML)',
       fill='HIV Status')+
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                limits = c(10^-0.5,10^3.5))+
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))

NETs

ggsave(NETs,filename="Figures/Figure 1/NETs.png",
       width = 7,height = 9,dpi = 1080)


# SAVE Figure 2
Fig2 <- (MPO|NETs)+
  plot_layout(ncol = 4,nrow = 3,widths = c(1,1,1,1))+
  #plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30,face = "bold"))

# Save Fig2
ggsave(filename = here("Figures/Figure 2","Fig2.pdf"),
       plot = Fig2,
       width = 30,height = 32,units = "in",dpi = 700)
# Save Fig2
ggsave(filename = here("Figures/Figure 2","Fig2.png"),
       plot = Fig2,
       width = 30,height = 32,units = "in",dpi = 700)

# Relationship between NER and Myeloeroxidase count
NERMPO <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         `HIV Status`!="NA") %>%
  ggplot(aes(`Myeloperoxidase (pg/mL)`,NER)) +
  geom_point(size=7, aes(color=`HIV Status`))+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman",size=12)+
  labs(y= 'Neutrophil to Epithelial cell ratio',
       x='Concentration of \n Myeloperoxidase (pg/ML)')+
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                limits = c(10^4,10^6), breaks=c(10^4,10^5,10^6))+
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                limits = c(10^-2,10^1), breaks=c(10^-2,10^-1,10^0,10^1))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
NERMPO
ggsave(NERMPO,filename="Figures/Figure 1/NERMPO.png",
       width = 7,height = 9,dpi = 1080)

# Relationship between NER and NETs concentration
NERNETs <- Masterfile %>%
  filter(`Epithelial count`>100,
         `CD45+ count`>200,
         `HIV Status`!="NA") %>%
  ggplot(aes(`NETs Concentration`,NER)) +
  geom_point(size=7)+
  geom_smooth(method = "lm")+
  stat_cor(method = "spearman",size=12)+
  labs(y= 'Neutrophil to Epithelial cell ratio',
       x='Concentration of NETs (pg/ML)')+
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                limits = c(10^-0.5,10^2.5))+
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                limits = c(10^-2,10^1), breaks=c(10^-2,10^-1,10^0,10^1))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background = element_rect(color='white'),
        legend.title = element_text(size=30,face = 'bold'),
        legend.text = element_text(size=25), 
        #axis.line = element_line(),
        #strip.text.x = element_text(size = 40,face = 'bold'),
        axis.text.x = element_text(size = 30),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
NERNETs
ggsave(NERNETs,filename="Figures/Figure 1/NERNETs.png",
       width = 7,height = 9,dpi = 1080)

# Correlogram using ggplot
Nasal_cytokines <- Masterfile %>%
  dplyr::filter(Visit=="Week 1",
         `CD11b Staining`=="Stained",
         `CD45+ count`>500,
         `Epithelial count`>100,
         ) %>%
  dplyr::select(`IFN-y`,`IL-10`,`IL-12p70`,`IL-13`,`IL-1B`,`IL-2`,`IL-4`,`IL-6`,`IL-8`,`TNF-a`,
         `HIV Status`) %>%
  na.omit()

# Convert only numeric columns, excludingg the specified colum (HIV Status)  
char_column <- "HIV Status"
Nasal_cytokines <- Nasal_cytokines %>%
  mutate(across(-all_of(char_column), ~ as.numeric(as.character(.))))
# Remove rows with any "0" values across numeric colums
Nasal_cytokines <- Nasal_cytokines %>%
  filter(if_all(-all_of(char_column), ~ . != 0))

# Step 1: Split the data by Antigen and calculate correlations and p-values for each group
cor_data <- Nasal_cytokines %>%
  dplyr::select(`IFN-y`,`IL-10`,`IL-12p70`,`IL-13`,`IL-1B`,`IL-2`,`IL-4`,`IL-6`,`IL-8`,`TNF-a`,
         `HIV Status`) %>%  # Select relevant columns
  group_by(`HIV Status`) %>%
  nest() %>%
  dplyr::mutate(cor_matrix = map(data, ~ {
    # Ensure only numeric columns are selected for correlation calculation
    numeric_data <- dplyr::select(.x, where(is.numeric))
    
    # Initialize matrices for correlation and p-values
    cor_matrix <- cor(numeric_data, method = "spearman")
    p_matrix <- matrix(NA, ncol = ncol(cor_matrix), nrow = nrow(cor_matrix))
    
    # Loop over combinations of variables to calculate p-values
    for (row in 1:nrow(cor_matrix)) {
      for (col in 1:ncol(cor_matrix)) {
        if (row != col) {
          test <- cor.test(numeric_data[[row]], numeric_data[[col]], method = "spearman")
          p_matrix[row, col] <- test$p.value
        } else {
          p_matrix[row, col] <- NA
        }
      }
    }
    
    list(correlation = cor_matrix, p_value = p_matrix)
  }))

# Step 2: Convert correlation matrices to long format with p-values
cor_long <- data.frame()

for (i in seq_len(nrow(cor_data))) {
  `HIV Status` <- cor_data$`HIV Status`[i]
  cor_matrix <- cor_data$cor_matrix[[i]]$correlation  # Extract the correlation matrix
  p_matrix <- cor_data$cor_matrix[[i]]$p_value  # Extract the p-value matrix
  
  if (!is.null(cor_matrix)) {  # Check if correlation matrix is not NULL
    # Convert the matrix to long format
    cor_df <- as.data.frame(as.table(cor_matrix))
    p_df <- as.data.frame(as.table(p_matrix))
    
    cor_df$`HIV Status` <- `HIV Status`  # Add Antigen column
    cor_df$p_value <- p_df$Freq  # Add p-value column
    
    # Combine the correlation data
    cor_long <- rbind(cor_long, cor_df)
  }
}

# Rename columns to make them suitable for plotting
colnames(cor_long) <- c("Variable1", "Variable2", "Correlation", "HIV Status", "p_value")

#cor_long_filtered <- cor_long %>%
#filter(as.numeric(factor(Variable1)) > as.numeric(factor(Variable2)))

cor_long <- cor_long %>%
  mutate(significance = case_when(
    p_value <= 0.05  ~ "*",     # One asterisk for p-value between 0.06 and 0.1
    p_value >= 0.01 & p_value < 0.05  ~ "**",    # Two asterisks for p-value between 0.03 and 0.05
    p_value >= 0.001 & p_value < 0.01  ~ "***",   # Three asterisks for p-value between 0.01 and 0.03
    p_value < 0.001                    ~ "****",  # Four asterisks for p-value less than 0.01
    TRUE                              ~ ""       # No asterisks for p-values greater than 0.1
  ))
# Step 3: Plot the faceted correlation matrix with renamed variables and significance markers
cor_matrix_plot <- ggplot(cor_long, aes(x = Variable1, y = Variable2, fill = Correlation)) +
  geom_tile(color = "white") +
  
  # Add correlation values inside the tiles
  geom_text(aes(label = round(Correlation, 2)), size = 3, color = "black") +
  
  # Add asterisks for significant correlations based on the significance column
  geom_text(aes(label = significance), size = 3, color = "black", vjust = -0.5) +
  
  scale_size(range = c(3, 10)) +  # Adjust the range for circle sizes
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Spearman\nCorrelation") +
  #labs(title = expression("Nasal cytokines"))+
  
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1,face = "bold"),
        axis.text.y = element_text(size = 8,face = "bold"),
        strip.text.x = element_text(size = 10,face = 'bold'),
        axis.title = element_blank(),
        plot.title = element_text(size = 15,face = "bold",hjust = 0.5)) +
  coord_fixed() +
  
  facet_wrap(~`HIV Status`)  # Facet by Antigen

# Print the plot
print(cor_matrix_plot)

ggsave(cor_matrix_plot,filename="Figures/Figure 3/Nasal correlogram.png",
       width = 14,height = 12,dpi = 1080)
ggsave(cor_matrix_plot,filename="Figures/Figure 3/Nasal correlogram.pdf",
       width = 14,height = 12,dpi = 1080)


# Only one half of the plot 
# Filter the data to keep only one half of the correlation matrix
cor_long <- cor_long %>%
  dplyr::filter(as.character(Variable1) < as.character(Variable2))

# Continue with the plot creation
cor_matrix_plot <- ggplot(cor_long, aes(x = Variable1, y = Variable2, fill = Correlation)) +
  geom_tile(color = "white") +
  
  # Add correlation values inside the tiles
  geom_text(aes(label = round(Correlation, 2)), size = 3, color = "black") +
  
  # Add asterisks for significant correlations based on the significance column
  geom_text(aes(label = significance), size = 3, color = "black", vjust = -0.5) +
  
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Spearman\nCorrelation") +
  
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 8, face = "bold"),
        strip.text.x = element_text(size = 10, face = 'bold'),
        axis.title = element_blank(),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
  coord_fixed() +
  
  facet_wrap(~`HIV Status`)
cor_matrix_plot

# Correlation plot using cor function
Nasal_cytokines_split <- split(Nasal_cytokines,
                               Nasal_cytokines$`HIV Status`)

cor_matrices <- lapply(Nasal_cytokines_split,function(df){
  cor(df %>% dplyr::select(-`HIV Status`))
})

par(mfrow = c(1, length(cor_matrices)))
col_grad <- colorRampPalette(c("blue","white","red"))
for (status in names(cor_matrices)) {
  corrplot(cor_matrices[[status]],
           method="circle",
           col = col_grad(200),
           type = "lower",
           title=paste("HIV Status:",status),
           mar=c(0,0,2,0),
           cl.lim=c(-1,1),
           addgrid.col = "black",
           tl.cex = 0.8,
           number.cex = 0.7)
}

par(mfrow = c(1, 1))




library(Hmisc)
library(corrplot)
library(dplyr)

# Split the data by HIV Status
Nasal_cytokines_split <- split(Nasal_cytokines,
                               Nasal_cytokines$`HIV Status`)

# Calculate correlation matrices and p-values
cor_results <- lapply(Nasal_cytokines_split, function(df) {
  data <- df %>% dplyr::select(-`HIV Status`)
  rcorr(as.matrix(data)) # Calculates correlations and p-values
})

# Set up plot layout: one plot for each HIV status
num_plots <- length(cor_results)
par(mfrow = c(1,3))  # Adjust for the number of statuses
col_grad <- colorRampPalette(c("blue", "white", "red"))

# Plot each correlation matrix with p-values
for (status in names(cor_results)) {
  cor_mat <- cor_results[[status]]$r   # Correlation coefficients
  p_mat <- cor_results[[status]]$P  # P-values
}

# Plot each correlation matrix with p-values
for (status in names(cor_results)) {
  cor_mat <- cor_results[[status]]$r   # Correlation coefficients
  p_mat <- cor_results[[status]]$P    # P-values
  
  # Add correlation plot
  corrplot(cor_mat,
           method = "circle",
           type = "lower",
           col = col_grad(200),
           title = paste("HIV Status:", status),
           mar = c(0, 0, 2, 0),
           cl.lim = c(-1, 1),
           addgrid.col = "black",
           tl.cex = 0.8,
           number.cex = 0.7,
           p.mat = p_mat,            # Add p-value matrix
           sig.level = 0.05,        # Highlight significant correlations
           insig = "label_sig")     # Show p-values as labels
}

# Reset plotting parameters
par(mfrow = c(1, 1))

png("Figures/corr_plot_Neg.png", width = 800,height = 800,dpi = 1080)
# Add correlation plot
corrplot(cor_results$`HIV-`$r,
         method = "circle",
         type = "lower",
         col = col_grad(200),
         title = paste("HIV-"),
         mar = c(0, 0, 2, 0),
         cl.lim = c(-1, 1),
         addgrid.col = "black",
         tl.cex = 1.5,
         number.cex = 0.7,
         p.mat = cor_Neg$P,            # Add p-value matrix
         sig.level = 0.05,        # Highlight significant correlations
         insig = "label_sig")     # Show p-values as label
dev.off()

# Add correlation plot
png("Figures/corr_plot_3M.png", width = 800,height = 800)
corrplot(cor_results$`PLHIV ART <3m`$r,
         method = "circle",
         type = "lower",
         col = col_grad(200),
         title = paste("ART<3m"),
         mar = c(0, 0, 2, 0),
         cl.lim = c(-1, 1),
         addgrid.col = "black",
         tl.cex = 1.5,
         number.cex = 0.7,
         p.mat = cor_results$`PLHIV ART <3m`$P,            # Add p-value matrix
         sig.level = 0.05,        # Highlight significant correlations
         insig = "label_sig")     # Show p-values as labels
dev.off()


# Add correlation plot
png("Figures/corr_plot_1Y.png", width = 800,height = 800)
corrplot(cor_results$`PLHIV ART >1y`$r,
         method = "circle",
         type = "lower",
         col = col_grad(200),
         title = paste("ART>1y"),
         mar = c(0, 0, 2, 0),
         cl.lim = c(-1, 1),
         addgrid.col = "black",
         tl.cex = 1.5,
         number.cex = 0.7,
         p.mat = cor_results$`PLHIV ART >1y`$P,            # Add p-value matrix
         sig.level = 0.05,        # Highlight significant correlations
         insig = "label_sig")     # Show p-values as labels
dev.off()

# Using ggplot
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare the data for ggplot2
cor_results_long <- lapply(names(cor_results), function(status) {
  cor_matrix <- cor_results[[status]]$cor
  p_matrix <- cor_results[[status]]$p
  
  # Create a long-format data frame for correlations
  cor_long <- as.data.frame(as.table(cor_matrix))
  colnames(cor_long) <- c("Var1", "Var2", "Correlation")
  
  # Add p-values
  cor_long$p_value <- as.vector(p_matrix)
  
  # Add significance markers
  cor_long <- cor_long %>%
    mutate(significance = case_when(
      p_value < 0.001 ~ "****",
      p_value < 0.01  ~ "***",
      p_value < 0.05  ~ "**",
      p_value < 0.1   ~ "*",
      TRUE            ~ ""
    ))
  
  # Add HIV status
  cor_long$HIV_Status <- status
  return(cor_long)
}) %>%
  bind_rows()

# Plot with ggplot2
Cytokines_cor <- ggplot(cor_results_long, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Correlation, 2)), size = 3, color = "black") +
  geom_text(aes(label = significance), size = 3, color = "black", vjust = -0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1), name = "Correlation") +
  facet_wrap(~HIV_Status) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12,face = 'bold'),
        axis.text.y = element_text(size = 12,face = 'bold'),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = 10)) +
  coord_fixed()
Cytokines_cor

ggsave(Cytokines_cor,filename="Figures/Figure 3/Nasal correlogram_new.png",
       width = 14,height = 12,dpi = 1080)
ggsave(Cytokines_cor,filename="Figures/Figure 3/Nasal correlogram_new.pdf",
       width = 14,height = 12,dpi = 1080)


Masterfile %>%
  dplyr::filter(Visit=="Week 1") %>%
  ggplot(aes(x=`Neutrophil proportion`,
             #x=`NER`,
             y=`IL-13`))+
  geom_point(size=7)+
  geom_smooth(method = "lm")+
  scale_x_log10()+
  scale_y_log10()+
  stat_cor(method = "spearman")+
  theme_bw()+
  theme(panel.grid = element_blank())

# CD11b Data -----
Neut_CD11b_data <- read_csv("Data/Main_Files_Thesis/NewMF_2024_Week1_NeutrophilCD11bCD63.csv") %>%
  merge(Micro, by=c("LAB ID"), all=FALSE) %>%
  dplyr::select(-Samples) %>%
  dplyr::mutate(`Neut CD11b+`=`Neut CD11b+CD63-`+`Neut CD11b+CD63+`,
                `Neut CD63+`=`Neut CD11b-CD63+`+`Neut CD11b+CD63+`)

HIV_status_colors <- c('#A9A9A9','#941100','#005493')

# Proportion of CD63+ Neutrophils
`CD63+_Neutrophils` <- Neut_CD11b_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`Neut CD63+`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p = {p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Proportion of CD63"^"+"~"Neutrophils"),
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
`CD63+_Neutrophils`

# Proportion of CD11b+ Neutrophils
`CD11b+_Neutrophils` <- Neut_CD11b_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`Neut CD11b+`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p = {p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Proportion of CD11b"^"+"~"Neutrophils"),
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
`CD11b+_Neutrophils`


# Proportion of CD11b+CD63- Neutrophils
`CD11b+CD63-_Neutrophils` <- Neut_CD11b_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`Neut CD11b+CD63-`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p = {p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Proportion of CD63"^"+"~"Neutrophils"),
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
`CD11b+CD63-_Neutrophils`

# Proportion of CD11b-CD63+ Neutrophils
`CD11b-CD63+_Neutrophils` <- Neut_CD11b_data %>%
  ggplot(aes(x=factor(`HIV Status`,levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y')),
             y=as.numeric(`Neut CD11b-CD63+`),
             color=factor(`HIV Status`,
                          levels=c('HIV-','PLHIV ART <3m', 'PLHIV ART >1y'))))+
  geom_boxplot(width=0.5,notch = F, outlier.shape = NA)+
  geom_jitter(position = position_jitter(width = 0.2),size=7) +
  #geom_errorbar(stat = 'summary',fun.data='mean_cl_normal',position = 'dodge',width=0.15,color='black')+
  stat_summary(geom = "point",
               fun = median,
               color='black',
               size=7,
               shape=95,
               position = position_dodge(width = 0.75))+
  geom_pwc(method = 't.test',
           label = "p = {p.format}",
           label.size = 8,
           tip.length = 0.01,
           hide.ns = F)+
  scale_color_manual(values = c('#A9A9A9','#941100','#005493'))+
  labs(x='',
       y=expression("Proportion of CD63"^"+"~"Neutrophils"),
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
        axis.text.x = element_text(size = 20, face = "bold",color = HIV_status_colors),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title = element_text(size = 30))
`CD11b-CD63+_Neutrophils`



# SINGLE CELL RNA SEQUENCING ------
# Load required packages
pacman::p_load(char = c("lubridate","gtsummary", "tidyverse", "dplyr", "here", "rio", "scales", "boot", "Matrix","ggpubr",
                        "magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl","DropletUtils","SeuratWrappers",
                        "Seurat","rjson","R2HTML","DT","cowplot","RCurl","glmGamPoi","DESeq2","ggrepel","rPanglaoDB","ComplexHeatmap",
                        "EnhancedVolcano","RColorBrewer","circlize","rmarkdown","biomaRt","biomartr","clusterProfiler","multinichenetr",
                        "AnnotationDbi","org.Hs.eg.db","CEMiTool","enrichplot","pathview","scmap","SingleR","S4Vectors","TMB","muscat",
                        "SingleCellExperiment","apeglm","edgeR","purrr","tibble","png","RColorBrewer","scran","microViz","CommPath",
                        "reshape2","scater","Azimuth","scCATCH","CellChat","SoupX","knitr","DoubletFinder","ggpubr","viridis","GSVA",
                        "ggmin","cluster","foreach","doParallel","BPCells","ggimage","ggbeeswarm","grid","data.table","scriabin",
                        "clusterExperiment","destiny","gam","corrplot","ggthemes","base64enc","Biobase","CATALYST","dittoSeq","viridis",
                        "DelayedArray","DelayedMatrixStats","limma","lme4","batchelor","HDF5Array","terra","ggrastr","nloptr","ggsignif",
                        'lubridate',"gtsummary", 'tidyverse', "dplyr", "here", "rio", "scales", "boot", 
                        "magrittr",  "mvtnorm", "zoo", "patchwork", "mgcv", "PropCIs", "writexl", "presto",
                        "ggsignif", "ggpubr", "ggeasy", "cowplot","ggExtra", "PupillometryR","hrbrthemes", "ggstance",
                        "survival","survminer","sysfonts","showtext","nlme",'glue'))

# install.packages("devtools")
devtools::install_github("saeyslab/nichenetr")
devtools::install_github("saeyslab/multinichenetr")
devtools::install_github('immunogenomics/presto')

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(here)


# Load all cell cluster
load("data/Single_Cell_Data/all_merged_subset_labelled.RData")
DimPlot(all_merged_subset_labelled,
        reduction = "umap.harmony",label = T,repel = T)+NoLegend()+
  labs(x="UMAP1",y="UMAP2")+
  theme_minimal()
# Renname idents
all_merged_subset_labelled <- RenameIdents(object = all_merged_subset_labelled,
                                           "Secretory cells" = "Secretory cells",
                                           "Goblet cells" = "Goblet cells",
                                           "CD3+ T cells" = "CD3+ T cells",
                                           "Phagocytes" = "Mono/Mac",
                                           "B cells" = "B cells",
                                           "Developing ciliated cells" = "Ciliated cells",
                                           "Neutrophils" = "Neutrophils",
                                           "Squamous cells" = 'Squamous cells',
                                           "FOXJ1++ Ciliated cells" = "Ciliated cells",
                                           "Stressed cells" = "Stressed cells",
                                           "B cells" = "B cells",
                                           "Club cells" = "Club cells",
                                           "Deuterosomal cells" = "Ciliated cells",
                                           "BEST++ Cilia++ Ciliated cells" = "Ciliated cells",
                                           "Ionocytes" = "Ionocytes",
                                           "Dendritic cells" = "Dendritic cells",
                                           "Neurons" = "Neurons")

all_merged_subset_labelled <- subset(all_merged_subset_labelled,
                                     idents = c("Stressed cells","Neurons","Club cells"),
                                     invert = T)
all_merged_subset_labelled$Clusters <- paste0(all_merged_subset_labelled@active.ident)
DimPlot(all_merged_subset_labelled,
        reduction = "umap.harmony",
        label = T,repel = T,
        #split.by = "HIV_Status"
)+NoLegend()+
  theme_minimal()+
  labs(x="UMAP-1",y="UMAP-2")+
  theme(#legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank())

# Define a consistent color palette for 17 clusters
Clusters <- unique(all_merged_subset_labelled$Clusters)
num_clusters <- 11
base_palette <- brewer.pal(min(num_clusters,11),'Paired')
#additional_colors <- colorRampPalette(brewer.pal(8,'Set3'))(num_clusters-12)
#color_palette <- c(base_palette,additional_colors)
#names(color_palette) <- Clusters
names(base_palette) <- Clusters
# Plot UMAP with consistent colors
# Figure 3a (Main UMAP)
main_umap <- DimPlot(all_merged_subset_labelled,
                     reduction = 'umap.harmony',
                     label = TRUE,
                     repel = TRUE,
                     label.size = 5)+
  scale_color_manual(values = base_palette)+
  #NoLegend()+
  labs(x='UMAP-1',
       y='UMAP-2')+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())
main_umap <- as.ggplot(main_umap)
main_umap
# Save Figure 3a
ggsave(main_umap,filename="Figures/Figure 3/Fig3a.png",
       width = 12,height = 9,dpi = 1080,units = "in")

ggsave(main_umap,filename="Figures/Figure 3/Fig3a.pdf",
       width = 12,height = 9,dpi = 1080, units = "in")


# Figure 3b (Pie chart)
# Extract cluster identities
cluster_info <- Idents(all_merged_subset_labelled)
clusters <- unique(cluster_info)
# Calculate proportion for each cluster
cluster_counts <- table(cluster_info)
cluster_df <- as.data.frame(cluster_counts)
colnames(cluster_df) <- c('Cluster','Count')
# Compute percentages
cluster_df$fraction <- cluster_df$Count / sum(cluster_df$Count)

# Compute the cumulative percentages (top of each rectangle)
cluster_df$ymax <- cumsum(cluster_df$fraction)

# Compute the bottom of each rectangle
cluster_df$ymin <- c(0, head(cluster_df$ymax, n=-1))

# Compute label position
cluster_df$labelPosition <- (cluster_df$ymax + cluster_df$ymin) / 2

# Compute a good label
cluster_df$label <- paste0(cluster_df$Cluster, "\n value: ", cluster_df$Count)

# Make the plot
all_data_pie <- ggplot(cluster_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Cluster
)) +
  geom_rect() +
  geom_point(aes(y = (ymax + ymin) / 2, x = 3.5, fill = Cluster), 
             shape = 21, size = 0, stroke = 0)+
  #geom_label(x=4,
  #aes(y=labelPosition, label=label), 
  #size=3) +
  scale_fill_manual(values = base_palette)+
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void()+
  theme(legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = 'bold'))+
  guides(fill=guide_legend(override.aes = list(shape=21, size=5,color=NA)))
all_data_pie

# Save Figure 3b
ggsave(all_data_pie,filename="Figures/Figure 3/Fig3b.png",
       width = 20,height = 7,dpi = 1080,units = "in")

ggsave(all_data_pie,filename="Figures/Figure 3/Fig3b.pdf",
       width = 20,height = 7,dpi = 1080, units = "in")


# Figure 3c (Heatmap)
# Heatmap for differential expressed gened for each cluster 
all_merged_subset_labelled <- FindVariableFeatures(all_merged_subset_labelled,
                                                   selection.method = 'vst',
                                                   verbose = TRUE)
genes <- VariableFeatures(all_merged_subset_labelled)
toplot <- SeuratExtend::CalcStats(all_merged_subset_labelled,
                                  features = genes,
                                  method = 'zscore',
                                  order = 'p',
                                  n = 5)
cluster_heatmap <-SeuratExtend::Heatmap(toplot,
                                        lab_fill='zscore',
                                        plot.margin = margin(l = 30),
                                        angle = 90)    
cluster_heatmap <- as.ggplot(cluster_heatmap)
cluster_heatmap
# Save Figure 3c
ggsave(cluster_heatmap,filename="Figures/Figure 3/Fig3c.png",
       width = 4,height = 8,dpi = 1080,units = "in")

ggsave(cluster_heatmap,filename="Figures/Figure 3/Fig3c.pdf",
       width = 4,height = 8,dpi = 1080, units = "in")
# ggplot 2 heatmap
gg_toplot <- as.data.frame(toplot) %>%
  mutate(feature = rownames(.)) %>%
  pivot_longer(cols = c(1:11),
               names_to = 'Cluster',
               values_to = 'zscore') %>%
  ggplot(aes(x=Cluster,
             y=feature,
             fill=zscore))+
    geom_tile(color='white')+
    scale_fill_gradient(low = "grey",high = "darkred")+
  theme_bw()+
  theme(plot.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 5),
        axis.title = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.7,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 5))+
  guides(fill = guide_colorbar(
    title = "zScore",
    title.theme = element_text(size = 9),
    barwidth = 0.8,
    barsize = 1))
gg_toplot

# Figure 3d (Main Bargraph)
# Calculate cluster proportion 
cluster_df <- as.data.frame(table(all_merged_subset_labelled$Clusters,
                                  all_merged_subset_labelled$HIV_Status))
colnames(cluster_df)<-c('Cluster','HIV Status','Count')

bargraph_main <- cluster_df %>%
  ggplot(aes(x=`HIV Status`,
             y=Count, fill=factor(Cluster,levels=c("Secretory cells","Goblet cells",
                                                   "CD3+ T cells","Mono/Mac","B cells",
                                                   "Ciliated cells","Neutrophils","Squamous cells",
                                                   "Ionocytes","Dendritic cells","Basal cells"))))+
  geom_col(position = "fill")+
  scale_fill_manual(values = base_palette)+
  scale_x_discrete(labels = c("HIV-"="HIV-",
                              "HIV+ ART<3 Months"="PLHIV on ART<3m",
                              "HIV+ ART>1 Year"="PLHIV on ART>1y"))+
  labs(x='',
       y='Frequency of cells')+
  theme_classic()+
  theme(legend.position = "right",
        legend.key.size = unit(0.3,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.text.y = element_text(angle = 0,hjust = 1,face = 'bold'))+
  coord_flip()
bargraph_main

# Save Figure 3d
ggsave(bargraph_main,filename="Figures/Figure 3/bargraph_main.png",
       width = 7,height = 2.3,dpi = 1080,units = "in")

ggsave(bargraph_main,filename="Figures/Figure 3/bargraph_main.pdf",
       width = 7,height = 2.3,dpi = 1080, units = "in")


# Perform GO analysis for each cluster
# Load DEG markers for each cluster
Markers <- FindAllMarkers(
  all_merged_subset_labelled,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = TRUE,
  min.cells.feature = 3,
  min.cells.group = 3)

write.csv(Markers,'scRNAseq_Results/Markers.csv',row.names = T)
#Markers <- read_csv('scRNAseq_Results/Markers.csv',row.name=T)
#column_to_rownames(Markers)

# Filter for significant DEGs (adjusted p-value <0.05)
significant_markers <- Markers %>%
  dplyr::filter(p_val_adj <= 0.05,abs(avg_log2FC) > 0.25)



# Split DEGs by cluster
DEGs_by_cluster <- split(significant_markers$gene,
                         significant_markers$cluster)


# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 5)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

GO_Terms_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
GO_Terms_Heatmap
GO_Terms_Heatmap <- as.ggplot(GO_Terms_Heatmap)

# Save Figure 3e
ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# ggplot 2 heatmap
cluster_order <- as.data.frame(colnames(heatmap_matrix))
colnames(cluster_order) <- "Cluster"

feature_order <- as.data.frame(rownames(heatmap_matrix))
colnames(feature_order) <- "features"

heatmap_matrix_new <- as.data.frame(heatmap_matrix)
heatmap_matrix_new$feature <- rownames(heatmap_matrix_new)
heatmap_matrix_new$wrapped_feature <- str_wrap(heatmap_matrix_new$feature, width = 50)

gg_heatmap <- as.data.frame(heatmap_matrix_new) %>%
  #mutate(feature = rownames(.)) %>%
  #str_wrap("feature", width = 15) %>%
  pivot_longer(cols = c(1:11),
               names_to = 'Cluster',
               values_to = 'zscore') %>%
  ggplot(aes(x=factor(Cluster,levels = cluster_order$Cluster),
             y=factor(feature,levels = feature_order$features),
             fill=zscore))+
  geom_tile(color='white')+
  scale_y_discrete(labels = heatmap_matrix_new$wrapped_feature)+
  scale_fill_gradient(low = "#FDFFFF",high = "darkred")+
  #scale_x_discrete(cluster_order)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 9),
        axis.title = element_blank(),
        #panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.7,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 5))+
  guides(fill = guide_colorbar(
    title = "zScore",
    title.theme = element_text(size = 9),
    barwidth = 0.8,
    barsize = 1))
gg_heatmap

# Save Figure 3e
ggsave(gg_heatmap,filename="Figures/Figure 3/gg_heatmap.png",
       width = 6,height = 11,dpi = 1080,units = "in")

ggsave(gg_heatmap,filename="Figures/Figure 3/gg_heatmap.pdf",
       width = 6,height = 11,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
GO_Terms_Heatmap <- SeuratExtend::Heatmap(heatmap_matrix,
                                          lab_fill='zscore',
                                          plot.margin = margin(l = 30),
                                          y_text_position = "right",
                                          angle = 90,
                                          hjust = 0.5,
                                          color_scheme = color_palette,
                                          vjust = 0.5)
GO_Terms_Heatmap <- as.ggplot(GO_Terms_Heatmap)
GO_Terms_Heatmap
# Save Figure 3e
ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e2.png",
       width = 11.5,height = 12,dpi = 1080,units = "in")

ggsave(GO_Terms_Heatmap,filename="Figures/Figure 3/Fig3e2.pdf",
       width = 11.5,height = 12,dpi = 1080, units = "in")

# Top10 highly significant genes in all the clusters
Top5_Markers <- Markers %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  dplyr::filter(avg_log2FC>1.0) %>%
  dplyr::filter(pct.1>0.2) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  slice_head(n=5)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("B cells","CD3+ T cells","Neutrophils","Dendritic cells",
                  "Mono/Mac","Secretory cells","Goblet cells","Ciliated cells",
                  "Squamous cells","Basal cells","Ionocytes")
all_merged_subset_labelled <- Seurat::SetIdent(all_merged_subset_labelled,value = factor(Idents(all_merged_subset_labelled),
                                                           levels = levels_order))
Top5_Markers_plot <- Seurat::DotPlot(all_merged_subset_labelled,
                                     features = unique(Top5_Markers$gene),
                                     cols = c("blue", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 7#,face = 'bold'
                                   ),
        axis.text.y = element_text(size = 9#, face = 'bold'
                                   ),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 1))
Top5_Markers_plot
ggsave(Top5_Markers_plot,filename="Figures/Figure 3/Fig3Dotplot.png",
       width = 10,height = 3.5,dpi = 1080,units = "in")

ggsave(Top5_Markers_plot,filename="Figures/Figure 3/Fig3Dotplot.pdf",
       width = 10,height = 3.5,dpi = 1080, units = "in")


# KEGG pathway and KEGG module analysis
significant_markers <- Markers %>%
  dplyr::filter(p_val_adj <= 0.05,abs(avg_log2FC) > 0.25)

gene_symbols <- unique(significant_markers$gene)

mapped_genes <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
  )
head(mapped_genes)


significant_markers <- merge(
  significant_markers,
  mapped_genes,
  by.x = "gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

head(significant_markers)

# Split DEGs by cluster
DEGs_by_cluster_KEGG <- split(significant_markers$ENTREZID,
                         significant_markers$cluster)


# Perform actual GO analysis
run_KEGG_analysis <- function(ENTREZID, cluster_name){
  enrichKEGG(
    gene = ENTREZID,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05
  )
}
# Perform GO analysis for all clusters
KEGG_results <- lapply(names(DEGs_by_cluster_KEGG),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_KEGG_analysis(DEGs_by_cluster_KEGG[[cluster]], cluster)
})

names(KEGG_results) <- names(DEGs_by_cluster_KEGG)
subsetted_KEGG_sucategories <- c("Cell growth and death","Cell motility","Cellular community - eukaryotes",
                                 "Immune disease","Immune system","Infectious disease: bacterial","Infectious disease: parasitic",
                                 "Infectious disease: viral","Signal transduction","Signalling molecules and interaction","Transport and catabolism") 
# Extract the @result slot and process data
zscore_data <- lapply(KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(KEGG_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(KEGG_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(KEGG_results)

# Subset the go_results object
KEGG_results_top20 <- lapply(KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 10)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
KEGG_results_top20 <- Filter(Negate(is.null), KEGG_results_top20)

# Check the first cluster's top 20 for example
head(KEGG_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(KEGG_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

KEGG_Terms_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = T,
  cluster_columns = F)
KEGG_Terms_Heatmap
KEGG_Terms_Heatmap <- as.ggplot(KEGG_Terms_Heatmap)

# Save Figure 3e
ggsave(KEGG_Terms_Heatmap,filename="Figures/Figure 3/KEGG_Terms_Heatmap.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(KEGG_Terms_Heatmap,filename="Figures/Figure 3/KEGG_Terms_Heatmap.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# ggplot 2 heatmap
cluster_order <- as.data.frame(colnames(heatmap_matrix))
colnames(cluster_order) <- "Cluster"

feature_order <- as.data.frame(rownames(heatmap_matrix))
colnames(feature_order) <- "features"

heatmap_matrix_new <- as.data.frame(heatmap_matrix)
heatmap_matrix_new$feature <- rownames(heatmap_matrix_new)
heatmap_matrix_new$wrapped_feature <- str_wrap(heatmap_matrix_new$feature, width = 50)

gg_KEGG_heatmap <- as.data.frame(heatmap_matrix_new) %>%
  pivot_longer(cols = c(1:11),
               names_to = 'Cluster',
               values_to = 'zscore') %>%
  ggplot(aes(x=factor(Cluster,levels = cluster_order$Cluster),
             y=factor(feature,levels = feature_order$features),
             fill=zscore))+
  geom_tile(color='white')+
  scale_y_discrete(labels = heatmap_matrix_new$wrapped_feature)+
  scale_fill_gradient(low = "#FDFFFF",high = "darkred")+
  #scale_x_discrete(cluster_order)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 9),
        axis.title = element_blank(),
        #panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.7,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 5))+
  guides(fill = guide_colorbar(
    title = "zScore",
    title.theme = element_text(size = 9),
    barwidth = 0.8,
    barsize = 1))
gg_KEGG_heatmap

# Save Figure 3e
ggsave(gg_KEGG_heatmap,filename="Figures/Figure 3/gg_KEGG_heatmap.png",
       width = 6,height = 11,dpi = 1080,units = "in")

ggsave(gg_KEGG_heatmap,filename="Figures/Figure 3/gg_KEGG_heatmap.pdf",
       width = 6,height = 11,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
KEGG_heatmap_SE <- SeuratExtend::Heatmap(heatmap_matrix,
                                          lab_fill='zscore',
                                          plot.margin = margin(l = 30),
                                          y_text_position = "right",
                                          angle = 90,
                                          hjust = 0.5,
                                          color_scheme = color_palette,
                                          vjust = 0.5)
KEGG_heatmap_SE <- as.ggplot(KEGG_heatmap_SE)
KEGG_heatmap_SE
# Save Figure 3e
ggsave(KEGG_heatmap_SE,filename="Figures/Figure 3/KEGG_heatmap_SE.png",
       width = 8,height = 16,dpi = 1080,units = "in")

ggsave(KEGG_heatmap_SE,filename="Figures/Figure 3/KEGG_heatmap_SE.pdf",
       width = 8,height = 16,dpi = 1080,units = "in")



# CEll TO CELL COMMUNICATION NETWORK USING MULTINICHENETR-----
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(nichenetr)
library(Seurat)

# Load the data
all_merged_subset_labelled_new <- readRDS("Data/Single_Cell_Data/all_merged_subset_labelled_new.rds")
# Load the organism
organism="human"
# Make gene names syntactically valid using make.names to avoid loss of genes and load the organism name
lr_network = readRDS("NicheNet/lr_network_human_allInfo_30112033.rds")%>%
  mutate(ligand = nichenetr::convert_alias_to_symbols(ligand, organism = organism),
         receptor = nichenetr::convert_alias_to_symbols(receptor, organism = organism)) %>% 
  dplyr::mutate(ligand = make.names(ligand), receptor = make.names(receptor)) %>% 
  dplyr::distinct(ligand, receptor)

ligand_target_matrix = readRDS("NicheNet/ligand_target_matrix_nsga2r_final.rds")
colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
  nichenetr::convert_alias_to_symbols(organism = organism) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
  nichenetr::convert_alias_to_symbols(organism = organism) %>% make.names()

lr_network = lr_network %>% dplyr::filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]

# Prepare the seurat object for conversion to single cell experiment object
all_merged_subset_labelled_new$Clusters <- paste0(all_merged_subset_labelled_new@active.ident)
all_merged_subset_labelled_new$HIV_Status <- paste0(all_merged_subset_labelled_new$HIV_Status)
all_merged_subset_labelled_new$sample <- paste0(all_merged_subset_labelled_new$sample)
# Make names syntactically valid
all_merged_subset_labelled_new$HIV_Status <- make.names(all_merged_subset_labelled_new$HIV_Status)
all_merged_subset_labelled_new$sample <- make.names(all_merged_subset_labelled_new$sample)
all_merged_subset_labelled_new$Clusters <- make.names(all_merged_subset_labelled_new$Clusters)
# Convert to a single cell experiment
sce <- Seurat::as.SingleCellExperiment(all_merged_subset_labelled_new)
SingleCellExperiment::colData(sce)
# Prepare cell-cell communication analysis
sample_id="sample"
group_id="HIV_Status"
celltype_id="Clusters"
covariates=NA
batches=NA
# Define contracts
contrasts_oi = c("'HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2','HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2','HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2'")
contrast_tbl = tibble(contrast =
                        c("HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2","HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2","HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2"),
                      group = c("HIV..ART.3.Months","HIV..ART.1.Year","HIV."))

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in%
            c(senders_oi,receivers_oi)]

conditions_keep = c("HIV..ART.3.Months","HIV..ART.1.Year","HIV.")
sce = sce[,SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep]
# Running Multinichenet core analysis
# Cell-type filtering
min_cells = 3
abundance_info =multinichenetr::get_abundance_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  batches = batches)

abundance_info$abund_plot_sample

abundance_df_summarized = abundance_info$abundance_data %>%
  dplyr::mutate(keep = as.logical(keep)) %>%
  dplyr::group_by(group_id,celltype_id) %>%
  dplyr::summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>%
  dplyr::filter(samples_present == 0) %>% 
  dplyr::pull(celltype_id) %>%
  unique()

celltypes_present_one_condition = abundance_df_summarized %>%
  dplyr::filter(samples_present >= 2) %>%
  dplyr::pull(celltype_id) %>%
  unique()

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>%
  unique() %>%
  length()

absent_celltypes <- abundance_df_summarized %>%
  dplyr::filter(samples_present < 2) %>% 
  dplyr::group_by(celltype_id) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>%
  dplyr::filter(n == total_nr_conditions) %>%
  pull(celltype_id)

analyse_condition_specific_celltypes = TRUE
if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% dplyr::setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% dplyr::setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    dplyr::setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    dplyr::setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

min_sample_prop = 0.50
fraction_cutoff = 0.05

frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% dplyr::pull(gene) %>% unique() 
sce = sce[genes_oi, ]

# Pseudobulk expression calculation
abundance_expression_info = multinichenetr::process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)

abundance_expression_info$celltype_info$pb_df %>% head()
abundance_expression_info$celltype_info$pb_df_group %>% head()
abundance_expression_info$sender_receiver_info$pb_df %>% head()
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()

# Differential expression (DE) analysis
DE_info = multinichenetr::get_DE_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
# Check DE results
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals


empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 

# Combine DE information for ligand-senders and receptors-receivers
sender_receiver_de = multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network)
sender_receiver_de %>% head(20)

# Ligand activity prediction: use DE analysis output to predict the activity
# of ligands in receiver cell types and infer their potential targer genes
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = FALSE 
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

# In case i want to use adjusted p values
geneset_assessment_adjustedPval = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj = TRUE, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment_adjustedPval

# Perform the ligand activity analysis and ligand-target inference
top_n_target = 250

verbose = TRUE
cores_system = 1
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 

# Running the ligand activity prediction
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  multinichenetr::get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores)))

ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

# Prioritization: Rank cell-cell communication patters through multi-criteria prioritization

ligand_activity_down = FALSE
sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(multinichenetr::generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down))
# Check output tables
prioritization_tables$group_prioritization_tbl %>% head(20)

# Calculate cross-samples expression correlation between ligand-receptor pairs and target genes
lr_target_prior_cor = multinichenetr::lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj)

# Save the output of multinichenetr
Path = "scRNAseq_Results/"
multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor) 
multinichenet_output = multinichenetr::make_lite_output(multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(Path, "multinichenet_output.rds"))
  
}

multinichenet_output <- readRDS("scRNAseq_Results/multinichenet_output.rds")

# Visualization of differential cell-cell interactions
prioritized_tbl_oi_all = multinichenetr::get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  top_n = 250, 
  rank_per_group = T)

prioritized_tbl_oi = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)



# Save the output of Circos_list
HIV_neg_circos <- circos_list$HIV.
HIV_1Y_circos <- circos_list$HIV..ART.1.Year
HIV_3M_circos <- circos_list$HIV..ART.3.Months

library(gridGraphics)
recordedplot_to_grob <- function(recorded_plot) {
  gridGraphics::grid.echo(recorded_plot)
  grid::grid.grab()
}

HIV_neg_circos <- recordedplot_to_grob(HIV_neg_circos)
HIV_3M_circos <- recordedplot_to_grob(HIV_3M_circos)
HIV_1Y_circos <- recordedplot_to_grob(HIV_1Y_circos)


library(ggplot2)

# Example: Convert grob to ggplot
ggplot_from_grob <- function(grob) {
  ggplot() +
    annotation_custom(grob) +
    theme_void()  # Remove axes, grid lines, etc.
}

# Convert your grob
HIV_neg_circos <- ggplot_from_grob(HIV_neg_circos)
HIV_3M_circos <- ggplot_from_grob(HIV_3M_circos)
HIV_1Y_circos <- ggplot_from_grob(HIV_1Y_circos)

# Combine Multiple ggplot Objects
library(patchwork)
combined_plot <- HIV_neg_circos + HIV_3M_circos + HIV_1Y_circos + plot_layout(ncol = 3)
print(combined_plot)

# Save Fig4
ggsave(combined_plot,filename="Figures/Figure 4/CElltoCell.png",
       width = 30,height = 8,dpi = 1080)
ggsave(combined_plot,filename="Figures/Figure 4/CElltoCell.pdf",
       width = 40,height = 10,dpi = 1080)

# Filtering the plots
filt_prioritized_tbl_oi <- prioritized_tbl_oi %>%
  dplyr::filter(receiver=="Neutrophils",sender!="CD3..T.cells",sender!="B.cells",sender!="Mono.Mac")
Neut_circos_list = make_circos_group_comparison(filt_prioritized_tbl_oi, colors_sender, colors_receiver)
Neg <- Neut_circos_list$HIV.
`3M` <- Neut_circos_list$HIV..ART.3.Months
`1Y` <- Neut_circos_list$HIV..ART.1.Year

library(gridGraphics)
recordedplot_to_grob <- function(recorded_plot) {
  gridGraphics::grid.echo(recorded_plot)
  grid::grid.grab()
}

Neg <- recordedplot_to_grob(Neg)
`3M` <- recordedplot_to_grob(`3M`)
`1Y` <- recordedplot_to_grob(`1Y`)

library(ggplot2)

# Example: Convert grob to ggplot
ggplot_from_grob <- function(grob) {
  ggplot() +
    annotation_custom(grob) +
    theme_void()  # Remove axes, grid lines, etc.
}

# Convert your grob
Neg <- ggplot_from_grob(Neg)
`3M` <- ggplot_from_grob(`3M`)
`1Y` <- ggplot_from_grob(`1Y`)
legend <- Neut_circos_list$legend

# Save Fig4
ggsave(Neg,filename="Figures/Figure 4/Neut_Neg_commun.png",
       width = 9,height = 8,dpi = 1080)
ggsave(`3M`,filename="Figures/Figure 4/Neut_3M_commun.png",
       width = 9,height = 8,dpi = 1080)
ggsave(`1Y`,filename="Figures/Figure 4/Neut_1Y_commun.png",
       width = 9,height = 8,dpi = 1080)

# Interpretable bubble plots
prioritized_tbl_oi_3M__50 <- filt_prioritized_tbl_oi %>%
  dplyr::filter(group=="HIV..ART.3.Months")

prioritized_tbl_oi_1Y__50 <- filt_prioritized_tbl_oi %>%
  dplyr::filter(group=="HIV..ART.1.Year")

plot_oi_3M <- multinichenetr::make_sample_lr_prod_activity_plots_Omnipath(
  multinichenet_output$prioritization_tables,
  prioritized_tbl_oi_3M__50 %>%dplyr::inner_join(lr_network)
)
````

# Further analysis of IMMUNE CELLS -----
Immune_cells <- subset(all_merged_subset_labelled,
                       idents = c("CD3+ T cells","Mono/Mac",
                                  "B cells","Neutrophils","Dendritic cells"),
                       invert = FALSE)
num_cells_immune_cells <- ncol(Immune_cells)
num_genes_immune_cells <- nrow(Immune_cells)
Immune_cells <- Seurat::NormalizeData(Immune_cells)
Immune_cells <- Seurat::FindVariableFeatures(Immune_cells)
Immune_cells <- Seurat::ScaleData(Immune_cells)
Immune_cells <- Seurat::RunPCA(Immune_cells, features = VariableFeatures(Immune_cells))
Immune_cells <- harmony::RunHarmony(
  object = Immune_cells, 
  group.by.vars = "sample",   # Column in metadata defining the grouping variable
  assay = "RNA",
  dims.use = 1:30# The assay to use (e.g., RNA)
)
Immune_cells <- Seurat::RunUMAP(Immune_cells, reduction = "harmony", dims = 1:30)
Immune_cells <- Seurat::FindNeighbors(Immune_cells, reduction = "harmony", dims = 1:30)
Immune_cells <- Seurat::FindClusters(Immune_cells, resolution = 0.5)
Main_umap_Immune <- Seurat::DimPlot(Immune_cells, reduction = "umap", group.by = "RNA_snn_res.0.5",label = T,repel = F)+NoLegend()+
  labs(x="UMAP-1",y="UMAP-2")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank())
Main_umap_Immune
Seurat::DotPlot(Immune_cells,
                features = unique(rownames(toplot)))+
  theme(axis.text.x = element_text(angle=90,hjust = 1, size = 7),
        axis.title = element_blank())

# Heatmap for differential expressed gened for each cluster 
#Immune_cells <- FindVariableFeatures(Immune_cells,selection.method = 'vst',verbose = TRUE)
Immune_genes <- VariableFeatures(Immune_cells)
Immune_toplot <- SeuratExtend::CalcStats(Immune_cells,
                                         features = genes,
                                         method = 'zscore',
                                         order = 'p',
                                         n = 5)
Immune_cluster_heatmap <-SeuratExtend::Heatmap(Immune_toplot,
                                               lab_fill='zscore',
                                               plot.margin = margin(l = 30),
                                               angle = 90)    
Immune_cluster_heatmap <- as.ggplot(Immune_cluster_heatmap)
Immune_cluster_heatmap
# Save Figure 3c
ggsave(Immune_cluster_heatmap,filename="Figures/Figure 3/Immune_cluster_heatmap.png",
       width = 4,height = 8,dpi = 1080,units = "in")

ggsave(Immune_cluster_heatmap,filename="Figures/Figure 3/Immune_cluster_heatmap.pdf",
       width = 4,height = 8,dpi = 1080, units = "in")

FeaturePlot(Immune_cells,
            features = c("CD3E","CD8A","TRAV10","KLRD1","KLRC2","KLRB1","CD4","NEAT1"),
            reduction = "umap",
            label = T)

VlnPlot(Immune_cells,
        features = c("CD3E","CD8A","TRAV10","KLRD1","KLRC2","KLRB1","CD4","NEAT1"))

Immune_cells <- RenameIdents(object = Immune_cells,
                             "0" = "Epithelial cells",
                             "1" = "CD8+ T cells",
                             "2" = "Neutrophils",
                             "3" = "Plasma B cells",
                             "4" = "Mac",
                             "5" = "Club cells",
                             "6" = "NK T cells",
                             "7" = "DCs",
                             "8" = "B cells",
                             "9" = "Mono",
                             "10" = "Mast cells",
                             "11" = "Epithelial cells")
# Remove the non immune cells
Immune_cells <- subset(Immune_cells,
                       idents = c("Epithelial cells","Club cells"),
                       invert = T)

# Dimplot of subsetted immune cells
Immune_cells$Clusters <- paste0(Immune_cells@active.ident)
Clusters <- unique(Immune_cells$Clusters)
num_clusters <- 9
base_palette <- brewer.pal(min(num_clusters,9),'Paired')
#additional_colors <- colorRampPalette(brewer.pal(8,'Set3'))(num_clusters-12)
#color_palette <- c(base_palette,additional_colors)
#names(color_palette) <- Clusters
names(base_palette) <- Clusters
# Plot UMAP with consistent colors
# Figure 3a (Main UMAP)
Main_umap_Immune <- Seurat::DimPlot(Immune_cells, 
                                    reduction = "umap",
                                    label = TRUE,
                                    repel = TRUE,
                                    label.size = 5)+
  NoLegend()+
  labs(x="UMAP-1",y="UMAP-2")+
  scale_color_manual(values = base_palette)+
  theme_bw()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())
Main_umap_Immune <- as.ggplot(Main_umap_Immune)
Main_umap_Immune
# Save Figure 3a
ggsave(Main_umap_Immune,filename="Figures/Figure 3/Main_umap_Immune.png",
       width = 8,height = 7,dpi = 1080,units = "in")

ggsave(Main_umap_Immune,filename="Figures/Figure 3/Main_umap_Immune.pdf",
       width = 8,height = 7,dpi = 1080, units = "in")

# Save the immune seurat object
save(Immune_cells, file = "Data/Single_Cell_Data/Immune_cells.RData")
# Load DEG markers for each cluster
Immune_Markers <- FindAllMarkers(
  Immune_cells,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = TRUE,
  min.cells.feature = 3,
  min.cells.group = 3)

write.csv(Immune_Markers,'scRNAseq_Results/Immune_Markers.csv',row.names = T)

## KEGG pathway and KEGG module analysis
Immune_sig_markers <- Immune_Markers %>%
  dplyr::filter(p_val_adj <= 0.05,abs(avg_log2FC) > 0.25)

gene_symbols <- unique(Immune_sig_markers$gene)

mapped_genes <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
head(mapped_genes)


Immune_sig_markers <- merge(
  Immune_sig_markers,
  mapped_genes,
  by.x = "gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

head(Immune_sig_markers)

# Split DEGs by cluster
Immune_DEGs_by_cluster_KEGG <- split(Immune_sig_markers$ENTREZID,
                                     Immune_sig_markers$cluster)


# Perform actual GO analysis
run_KEGG_analysis <- function(ENTREZID, cluster_name){
  enrichKEGG(
    gene = ENTREZID,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05
  )
}
# Perform KEGG analysis for all clusters
Immune_KEGG_results <- lapply(names(Immune_DEGs_by_cluster_KEGG),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_KEGG_analysis(Immune_DEGs_by_cluster_KEGG[[cluster]], cluster)
})

names(Immune_KEGG_results) <- names(Immune_DEGs_by_cluster_KEGG)
subsetted_KEGG_sucategories <- c("Cell growth and death","Cell motility","Cellular community - eukaryotes",
                                 "Immune disease","Immune system","Infectious disease: bacterial","Infectious disease: parasitic",
                                 "Infectious disease: viral","Signal transduction","Signalling molecules and interaction","Transport and catabolism") 
# Extract the @result slot and process data
zscore_data <- lapply(Immune_KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(Immune_KEGG_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(Immune_KEGG_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(Immune_KEGG_results)

# Subset the go_results object
KEGG_results_top20 <- lapply(Immune_KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 10)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
KEGG_results_top20 <- Filter(Negate(is.null), KEGG_results_top20)

# Check the first cluster's top 20 for example
head(KEGG_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(KEGG_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Immune_KEGG_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Immune_KEGG_Heatmap <- as.ggplot(Immune_KEGG_Heatmap)
Immune_KEGG_Heatmap

# Save Figure 3e
ggsave(Immune_KEGG_Heatmap,filename="Figures/Figure 3/Immune_KEGG_Heatmap.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Immune_KEGG_Heatmap,filename="Figures/Figure 3/Immune_KEGG_Heatmap.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# ggplot 2 heatmap
cluster_order <- as.data.frame(colnames(heatmap_matrix))
colnames(cluster_order) <- "Cluster"

feature_order <- as.data.frame(rownames(heatmap_matrix))
colnames(feature_order) <- "features"

heatmap_matrix_new <- as.data.frame(heatmap_matrix)
heatmap_matrix_new$feature <- rownames(heatmap_matrix_new)
heatmap_matrix_new$wrapped_feature <- str_wrap(heatmap_matrix_new$feature, width = 50)

Immune_gg_KEGG_heatmap <- as.data.frame(heatmap_matrix_new) %>%
  pivot_longer(cols = c(1:9),
               names_to = 'Cluster',
               values_to = 'zscore') %>%
  ggplot(aes(x=factor(Cluster,levels = cluster_order$Cluster),
             y=factor(feature,levels = feature_order$features),
             fill=zscore))+
  geom_tile(color='white')+
  scale_y_discrete(labels = heatmap_matrix_new$wrapped_feature)+
  scale_fill_gradient(low = "#FDFFFF",high = "darkred")+
  #scale_x_discrete(cluster_order)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 9),
        axis.title = element_blank(),
        #panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.7,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 5))+
  guides(fill = guide_colorbar(
    title = "zScore",
    title.theme = element_text(size = 9),
    barwidth = 0.8,
    barsize = 1))
Immune_gg_KEGG_heatmap

# Save Figure 3e
ggsave(Immune_gg_KEGG_heatmap,filename="Figures/Figure 3/Immune_gg_KEGG_heatmap.png",
       width = 6,height = 11,dpi = 1080,units = "in")

ggsave(Immune_gg_KEGG_heatmap,filename="Figures/Figure 3/Immune_gg_KEGG_heatmap.pdf",
       width = 6,height = 11,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Immune_KEGG_heatmap_SE <- SeuratExtend::Heatmap(heatmap_matrix,
                                         lab_fill='zscore',
                                         plot.margin = margin(l = 30),
                                         y_text_position = "right",
                                         angle = 90,
                                         hjust = 0.5,
                                         color_scheme = color_palette,
                                         vjust = 0.5)
Immune_KEGG_heatmap_SE <- as.ggplot(Immune_KEGG_heatmap_SE)
Immune_KEGG_heatmap_SE
# Save Figure 3e
ggsave(Immune_KEGG_heatmap_SE,filename="Figures/Figure 3/Immune_KEGG_heatmap_SE.png",
       width = 8,height = 16,dpi = 1080,units = "in")

ggsave(Immune_KEGG_heatmap_SE,filename="Figures/Figure 3/Immune_KEGG_heatmap_SE.pdf",
       width = 8,height = 16,dpi = 1080,units = "in")

# CEll TO CELL COMMUNICATION NETWORK AMONG IMMUNE CELLS USING MULTINICHENETR-----
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(nichenetr)
library(Seurat)

load("Data/Single_Cell_Data/Immune_cells.RData")
organism="human"
# Make gene names syntactically valid using make.names to avoid loss of genes and load the organism name
lr_network = readRDS("NicheNet/lr_network_human_allInfo_30112033.rds")%>%
  mutate(ligand = nichenetr::convert_alias_to_symbols(ligand, organism = organism),
         receptor = nichenetr::convert_alias_to_symbols(receptor, organism = organism)) %>% 
  dplyr::mutate(ligand = make.names(ligand), receptor = make.names(receptor)) %>% 
  dplyr::distinct(ligand, receptor)

ligand_target_matrix = readRDS("NicheNet/ligand_target_matrix_nsga2r_final.rds")
colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% 
  nichenetr::convert_alias_to_symbols(organism = organism) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% 
  nichenetr::convert_alias_to_symbols(organism = organism) %>% make.names()

lr_network = lr_network %>% dplyr::filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix = ligand_target_matrix[, lr_network$ligand %>% unique()]

# Prepare the seurat object for conversion to single cell experiment object
Immune_cells$Clusters <- paste0(Immune_cells@active.ident)
Immune_cells$HIV_Status <- paste0(Immune_cells$HIV_Status)
Immune_cells$sample <- paste0(Immune_cells$sample)
# Make names syntactically valid
Immune_cells$HIV_Status <- make.names(Immune_cells$HIV_Status)
Immune_cells$sample <- make.names(Immune_cells$sample)
Immune_cells$Clusters <- make.names(Immune_cells$Clusters)
# Convert to a single cell experiment
sce <- Seurat::as.SingleCellExperiment(Immune_cells)
SingleCellExperiment::colData(sce)
# Prepare cell-cell communication analysis
sample_id="sample"
group_id="HIV_Status"
celltype_id="Clusters"
covariates=NA
batches=NA
# Define contracts
contrasts_oi = c("'HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2','HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2','HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2'")
contrast_tbl = tibble(contrast =
                        c("HIV..ART.3.Months-(HIV..ART.1.Year+HIV.)/2","HIV..ART.1.Year-(HIV..ART.3.Months+HIV.)/2","HIV.-(HIV..ART.1.Year+HIV..ART.3.Months)/2"),
                      group = c("HIV..ART.3.Months","HIV..ART.1.Year","HIV."))

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in%
            c(senders_oi,receivers_oi)]

conditions_keep = c("HIV..ART.3.Months","HIV..ART.1.Year","HIV.")
sce = sce[,SummarizedExperiment::colData(sce)[,group_id] %in% 
            conditions_keep]
# Running Multinichenet core analysis
# Cell-type filtering
min_cells = 3
abundance_info =multinichenetr::get_abundance_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  batches = batches)

abundance_info$abund_plot_sample

abundance_df_summarized = abundance_info$abundance_data %>%
  dplyr::mutate(keep = as.logical(keep)) %>%
  dplyr::group_by(group_id,celltype_id) %>%
  dplyr::summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>%
  dplyr::filter(samples_present == 0) %>% 
  dplyr::pull(celltype_id) %>%
  unique()

celltypes_present_one_condition = abundance_df_summarized %>%
  dplyr::filter(samples_present >= 2) %>%
  dplyr::pull(celltype_id) %>%
  unique()

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>%
  unique() %>%
  length()

absent_celltypes <- abundance_df_summarized %>%
  dplyr::filter(samples_present < 2) %>% 
  dplyr::group_by(celltype_id) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>%
  dplyr::filter(n == total_nr_conditions) %>%
  pull(celltype_id)

analyse_condition_specific_celltypes = TRUE
if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% dplyr::setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% dplyr::setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    dplyr::setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    dplyr::setdiff(union(absent_celltypes, condition_specific_celltypes))
}

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
]

min_sample_prop = 0.50
fraction_cutoff = 0.05

frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

genes_oi = frq_list$expressed_df %>% 
  filter(expressed == TRUE) %>% dplyr::pull(gene) %>% unique() 
sce = sce[genes_oi, ]

# Pseudobulk expression calculation
abundance_expression_info = multinichenetr::process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, 
  receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)

abundance_expression_info$celltype_info$pb_df %>% head()
abundance_expression_info$celltype_info$pb_df_group %>% head()
abundance_expression_info$sender_receiver_info$pb_df %>% head()
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()

# Differential expression (DE) analysis
DE_info = multinichenetr::get_DE_info(
  sce = sce, 
  sample_id = sample_id, 
  group_id = group_id, 
  celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)
# Check DE results
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals


empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 

# Combine DE information for ligand-senders and receptors-receivers
sender_receiver_de = multinichenetr::combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network)
sender_receiver_de %>% head(20)

# Ligand activity prediction: use DE analysis output to predict the activity
# of ligands in receiver cell types and infer their potential targer genes
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = FALSE 
geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment

# In case i want to use adjusted p values
geneset_assessment_adjustedPval = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj = TRUE, p_val_threshold
  ) %>% 
  bind_rows() 
geneset_assessment_adjustedPval

# Perform the ligand activity analysis and ligand-target inference
top_n_target = 250

verbose = TRUE
cores_system = 1
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 

# Running the ligand activity prediction
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  multinichenetr::get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores)))

ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

# Prioritization: Rank cell-cell communication patters through multi-criteria prioritization

ligand_activity_down = FALSE
sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(multinichenetr::generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  scenario = "regular", # all prioritization criteria will be weighted equally
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender,
  ligand_activity_down = ligand_activity_down))
# Check output tables
prioritization_tables$group_prioritization_tbl %>% head(20)

# Calculate cross-samples expression correlation between ligand-receptor pairs and target genes
lr_target_prior_cor = multinichenetr::lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj)

# Save the output of multinichenetr
Path = "scRNAseq_Results/"
Immune_multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor) 
Immune_multinichenet_output = multinichenetr::make_lite_output(Immune_multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(Immune_multinichenet_output, paste0(Path, "Immune_multinichenet_output.rds"))
  
}

Immune_multinichenet_output <- readRDS("scRNAseq_Results/Immune_multinichenet_output.rds")

# Visualization of differential cell-cell interactions
prioritized_tbl_oi_all = multinichenetr::get_top_n_lr_pairs(
  Immune_multinichenet_output$prioritization_tables, 
  top_n = 10, 
  rank_per_group = T)

prioritized_tbl_oi = 
  Immune_multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% 
  left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)



# Save the output of Circos_list
HIV_neg_circos <- circos_list$HIV.
HIV_1Y_circos <- circos_list$HIV..ART.1.Year
HIV_3M_circos <- circos_list$HIV..ART.3.Months

library(gridGraphics)
recordedplot_to_grob <- function(recorded_plot) {
  gridGraphics::grid.echo(recorded_plot)
  grid::grid.grab()
}

HIV_neg_circos <- recordedplot_to_grob(HIV_neg_circos)
HIV_3M_circos <- recordedplot_to_grob(HIV_3M_circos)
HIV_1Y_circos <- recordedplot_to_grob(HIV_1Y_circos)


library(ggplot2)

# Example: Convert grob to ggplot
ggplot_from_grob <- function(grob) {
  ggplot() +
    annotation_custom(grob) +
    theme_void()  # Remove axes, grid lines, etc.
}

# Convert your grob
HIV_neg_circos <- ggplot_from_grob(HIV_neg_circos)
HIV_3M_circos <- ggplot_from_grob(HIV_3M_circos)
HIV_1Y_circos <- ggplot_from_grob(HIV_1Y_circos)

# Save Fig4
ggsave(HIV_neg_circos,filename="Figures/Figure 4/Immune_Neg_CelltoCellall.png",
       width = 6,height = 8,dpi = 1080)
ggsave(HIV_3M_circos,filename="Figures/Figure 4/Immune_3M_CelltoCellall.png",
       width = 6,height = 8,dpi = 1080)
ggsave(HIV_1Y_circos,filename="Figures/Figure 4/Immune_1Y_CelltoCellall.png",
       width = 6,height = 8,dpi = 1080)


# Combine Multiple ggplot Objects
library(patchwork)
combined_plot <- HIV_neg_circos + HIV_3M_circos + HIV_1Y_circos + plot_layout(ncol = 3)
print(combined_plot)

# Save Fig4
ggsave(combined_plot,filename="Figures/Figure 4/ImmuneCElltoCell.png",
       width = 30,height = 8,dpi = 1080)
ggsave(combined_plot,filename="Figures/Figure 4/ImmuneCElltoCell.pdf",
       width = 40,height = 10,dpi = 1080)

# Filtering the plots
filt_prioritized_tbl_oi <- prioritized_tbl_oi %>%
  dplyr::filter(receiver=="Neutrophils")
Neut_circos_list = make_circos_group_comparison(filt_prioritized_tbl_oi, colors_sender, colors_receiver)
Neg <- Neut_circos_list$HIV.
`3M` <- Neut_circos_list$HIV..ART.3.Months
`1Y` <- Neut_circos_list$HIV..ART.1.Year

library(gridGraphics)
recordedplot_to_grob <- function(recorded_plot) {
  gridGraphics::grid.echo(recorded_plot)
  grid::grid.grab()
}

Neg <- recordedplot_to_grob(Neg)
`3M` <- recordedplot_to_grob(`3M`)
`1Y` <- recordedplot_to_grob(`1Y`)

library(ggplot2)

# Example: Convert grob to ggplot
ggplot_from_grob <- function(grob) {
  ggplot() +
    annotation_custom(grob) +
    theme_void()  # Remove axes, grid lines, etc.
}

# Convert your grob
Neg <- ggplot_from_grob(Neg)
`3M` <- ggplot_from_grob(`3M`)
`1Y` <- ggplot_from_grob(`1Y`)

# Save Fig4
ggsave(Neg,filename="Figures/Figure 4/Immune_Neg_commun.png",
       width = 7,height = 8,dpi = 1080)
ggsave(`3M`,filename="Figures/Figure 4/Immune_3M_commun.png",
       width = 7,height = 8,dpi = 1080)
ggsave(`1Y`,filename="Figures/Figure 4/Immune_1Y_commun.png",
       width = 7,height = 8,dpi = 1080)

# Interpretable bubble plots
prioritized_tbl_oi_3M__50 <- prioritized_tbl_oi_all %>%
  dplyr::filter(group=="HIV..ART.3.Months")

prioritized_tbl_oi_1Y__50 <- prioritized_tbl_oi_all %>%
  dplyr::filter(group=="HIV..ART.1.Year")

plot_oi_3M <- multinichenetr::make_sample_lr_prod_activity_plots_Omnipath(
  Immune_multinichenet_output$prioritization_tables,
  prioritized_tbl_oi_3M__50 %>%dplyr::inner_join(lr_network)
)
````



# Further analysis of EPITHELIAL CELLS -----
all_merged_subset_labelled_new <- readRDS("Data/Single_Cell_Data/all_merged_subset_labelled_new.rds")
Epithelial_cells <- subset(all_merged_subset_labelled_new,
                       idents = c("CD3+ T cells","Mono/Mac",
                                  "B cells","Neutrophils","Dendritic cells"),
                       invert = TRUE)
num_cells_epithelial_cells <- ncol(Epithelial_cells)
num_cells_epithelial_cells
num_genes_epithelial_cells <- nrow(Epithelial_cells)
num_genes_epithelial_cells
Epithelial_cells <- Seurat::NormalizeData(Epithelial_cells)
Epithelial_cells <- Seurat::FindVariableFeatures(Epithelial_cells)
Epithelial_cells <- Seurat::ScaleData(Epithelial_cells)
Epithelial_cells <- Seurat::RunPCA(Epithelial_cells, features = VariableFeatures(Epithelial_cells))
Epithelial_cells <- harmony::RunHarmony(
  object = Epithelial_cells, 
  group.by.vars = "sample",   # Column in metadata defining the grouping variable
  assay = "RNA",
  dims.use = 1:30# The assay to use (e.g., RNA)
)
Epithelial_cells <- Seurat::RunUMAP(Epithelial_cells, reduction = "harmony", dims = 1:30)
Epithelial_cells <- Seurat::FindNeighbors(Epithelial_cells, reduction = "harmony", dims = 1:30)
Epithelial_cells <- Seurat::FindClusters(Epithelial_cells, resolution = c(.1,.2,.3,.4,.5,.6,.7))
Main_umap_Epithelial <- Seurat::DimPlot(Epithelial_cells, reduction = "umap", group.by = "RNA_snn_res.0.5",label = T,repel = F)+NoLegend()+
  labs(x="UMAP-1",y="UMAP-2")+
  NoLegend()+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank())
Main_umap_Epithelial
# Load DEG markers for each cluster
Epithelial_Markers <- FindAllMarkers(
  Epithelial_cells,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = TRUE,
  min.cells.feature = 3,
  min.cells.group = 3)

write.csv(Epithelial_Markers,'scRNAseq_Results/Epithelial_Markers.csv',row.names = F)
Epithelial_Markers <- read_csv("scRNAseq_Results/Epithelial_Markers.csv")
Epithelial_cell_marker_genes <- read_csv("scRNAseq_Results/Epithelial_Cell_marker_genes.csv")
Epithelial_Dotplot <- Seurat::DotPlot(Epithelial_cells,
                features = unique(Epithelial_cell_marker_genes$`Selected human marker genes`))+
  NoLegend()+
  theme(axis.text.x = element_text(angle=90,hjust = 1, size = 7),
        axis.title = element_blank())
Epithelial_Dotplot

# Heatmap for differential expressed gened for each cluster 
#Immune_cells <- FindVariableFeatures(Immune_cells,selection.method = 'vst',verbose = TRUE)
Epithelial_genes <- VariableFeatures(Epithelial_cells)
Epithelial_toplot <- SeuratExtend::CalcStats(Epithelial_cells,
                                         features = genes,
                                         method = 'zscore',
                                         order = 'p',
                                         n = 5)
Epithelial_cluster_heatmap <-SeuratExtend::Heatmap(Epithelial_toplot,
                                               lab_fill='zscore',
                                               plot.margin = margin(l = 30),
                                               angle = 90)    
Epithelial_cluster_heatmap <- as.ggplot(Epithelial_cluster_heatmap)
Epithelial_cluster_heatmap
# Save Figure 3c
ggsave(Epithelial_cluster_heatmap,filename="Figures/Figure 3/Epithelial_cluster_heatmap.png",
       width = 4,height = 8,dpi = 1080,units = "in")

ggsave(Epithelial_cluster_heatmap,filename="Figures/Figure 3/Epithelial_cluster_heatmap.pdf",
       width = 4,height = 8,dpi = 1080, units = "in")


library(celldex)
ref.data <- celldex::HumanPrimaryCellAtlasData()
ref.data
library(SingleR)
# Extract the assay data from seurat object
Epithelial_cells_matrix <- Seurat::GetAssayData(Epithelial_cells,
                                                assay = "RNA",layer = "data")
pred.epith.cells <- SingleR::SingleR(test = Epithelial_cells_matrix,
                                  ref = ref.data,
                                  assay.type.test = "RNA",
                                  labels = ref.data$label.fine)
pred.epith.cells
# Renaming idents for epithelial cell clusters
Epithelial_cells <- RenameIdents(object = Epithelial_cells,
                             "0" = "",
                             "1" = "",
                             "2" = "Suprabasal cells",
                             "3" = "Ciliated cells?",
                             "4" = "",
                             "5" = "Hillock cells",
                             "6" = "Mucous-ciliated cells",
                             "7" = "",
                             "8" = "Basal cells",
                             "9" = "",
                             "10" = "Deuterosomal cells",
                             "11" = "Ionocytes")

VlnPlot(Epithelial_cells,
        features = c("SCGB1A1","MUC5AC"))
# Remove the immune cells
Epithelial_cells <- subset(Epithelial_cells,
                       idents = c("Epithelial cells","Club cells"),
                       invert = T)

# Dimplot of subsetted immune cells
Main_umap_Epithelial <- Seurat::DimPlot(Epithelial_cells, reduction = "umap",label = T,repel = F)+NoLegend()+
  labs(x="UMAP-1",y="UMAP-2")+
  scale_color_manual(values = base_palette)+
  theme_minimal()+
  theme(#panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_blank())
Main_umap_Epithelial

# Load DEG markers for each cluster
Epithelial_Markers <- FindAllMarkers(
  Epithelial_cells,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = TRUE,
  min.cells.feature = 3,
  min.cells.group = 3)

write.csv(Epithelial_Markers,'scRNAseq_Results/Epithelial_Markers.csv',row.names = T)

## KEGG pathway and KEGG module analysis
Epithelial_sig_markers <- Epithelial_Markers %>%
  dplyr::filter(p_val_adj <= 0.05,abs(avg_log2FC) > 0.25)

gene_symbols <- unique(Epithelial_sig_markers$gene)

mapped_genes <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
head(mapped_genes)


Epithelial_sig_markers <- merge(
  Immune_sig_markers,
  mapped_genes,
  by.x = "gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

head(Epithelial_sig_markers)

# Split DEGs by cluster
Epithelial_DEGs_by_cluster_KEGG <- split(Epithelial_sig_markers$ENTREZID,
                                         Epithelial_sig_markers$cluster)


# Perform actual GO analysis
run_KEGG_analysis <- function(ENTREZID, cluster_name){
  enrichKEGG(
    gene = ENTREZID,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05
  )
}
# Perform KEGG analysis for all clusters
Epithelial_KEGG_results <- lapply(names(Epithelial_DEGs_by_cluster_KEGG),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_KEGG_analysis(Epithelial_DEGs_by_cluster_KEGG[[cluster]], cluster)
})

names(Epithelial_KEGG_results) <- names(Epithelial_DEGs_by_cluster_KEGG)
subsetted_KEGG_sucategories <- c("Cell growth and death","Cell motility","Cellular community - eukaryotes",
                                 "Immune disease","Immune system","Infectious disease: bacterial","Infectious disease: parasitic",
                                 "Infectious disease: viral","Signal transduction","Signalling molecules and interaction","Transport and catabolism") 
# Extract the @result slot and process data
zscore_data <- lapply(Epithelial_KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(Epithelial_KEGG_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(Epithelial_KEGG_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(Epithelial_KEGG_results)

# Subset the go_results object
KEGG_results_top20 <- lapply(Epithelial_KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 10)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
KEGG_results_top20 <- Filter(Negate(is.null), KEGG_results_top20)

# Check the first cluster's top 20 for example
head(KEGG_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(KEGG_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Epithelial_KEGG_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Epithelial_KEGG_Heatmap <- as.ggplot(Epithelial_KEGG_Heatmap)
Epithelial_KEGG_Heatmap

# Save Figure 3e
ggsave(Epithelial_KEGG_Heatmap,filename="Figures/Figure 3/Epithelial_KEGG_Heatmap.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Epithelial_KEGG_Heatmap,filename="Figures/Figure 3/Epithelial_KEGG_Heatmap.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# ggplot 2 heatmap
cluster_order <- as.data.frame(colnames(heatmap_matrix))
colnames(cluster_order) <- "Cluster"

feature_order <- as.data.frame(rownames(heatmap_matrix))
colnames(feature_order) <- "features"

heatmap_matrix_new <- as.data.frame(heatmap_matrix)
heatmap_matrix_new$feature <- rownames(heatmap_matrix_new)
heatmap_matrix_new$wrapped_feature <- str_wrap(heatmap_matrix_new$feature, width = 50)

Epithelial_gg_KEGG_heatmap <- as.data.frame(heatmap_matrix_new) %>%
  pivot_longer(cols = c(1:9),
               names_to = 'Cluster',
               values_to = 'zscore') %>%
  ggplot(aes(x=factor(Cluster,levels = cluster_order$Cluster),
             y=factor(feature,levels = feature_order$features),
             fill=zscore))+
  geom_tile(color='white')+
  scale_y_discrete(labels = heatmap_matrix_new$wrapped_feature)+
  scale_fill_gradient(low = "#FDFFFF",high = "darkred")+
  #scale_x_discrete(cluster_order)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 9),
        axis.title = element_blank(),
        #panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.7,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 5))+
  guides(fill = guide_colorbar(
    title = "zScore",
    title.theme = element_text(size = 9),
    barwidth = 0.8,
    barsize = 1))
Epithelial_gg_KEGG_heatmap

# Save Figure 3e
ggsave(Epithelial_gg_KEGG_heatmap,filename="Figures/Figure 3/Epithelial_gg_KEGG_heatmap.png",
       width = 6,height = 11,dpi = 1080,units = "in")

ggsave(Epithelial_gg_KEGG_heatmap,filename="Figures/Figure 3/Epithelial_gg_KEGG_heatmap.pdf",
       width = 6,height = 11,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Epithelial_KEGG_heatmap_SE <- SeuratExtend::Heatmap(heatmap_matrix,
                                                lab_fill='zscore',
                                                plot.margin = margin(l = 30),
                                                y_text_position = "right",
                                                angle = 90,
                                                hjust = 0.5,
                                                color_scheme = color_palette,
                                                vjust = 0.5)
Epithelial_KEGG_heatmap_SE <- as.ggplot(Epithelial_KEGG_heatmap_SE)
Epithelial_KEGG_heatmap_SE
# Save Figure 3e
ggsave(Epithelial_KEGG_heatmap_SE,filename="Figures/Figure 3/Epithelial_KEGG_heatmap_SE.png",
       width = 8,height = 16,dpi = 1080,units = "in")

ggsave(Epithelial_KEGG_heatmap_SE,filename="Figures/Figure 3/Epithelial_KEGG_heatmap_SE.pdf",
       width = 8,height = 16,dpi = 1080,units = "in")

# FURTHER ANALYSIS OF NEUTROPHILS ----
Neutrophils <- subset(all_merged_subset_labelled_new,
                      idents = c("Neutrophils"),
                      invert = FALSE)
# Setting HIV Status as idents
Idents(Neutrophils) <- Neutrophils$HIV_Status
# Finding Neutrophil markers in
Neut_Markers <- FindAllMarkers(
  Neutrophils,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = F,
  min.cells.feature = 3,
  min.cells.group = 3)
# Save Neut_Markers
write.csv(Neut_Markers,'scRNAseq_Results/Neut_Markers.csv',row.names = T)
# Filtering for significant Neut_Markers
Sig_Neut_Markers <- Neut_Markers %>%
  dplyr::filter(abs(avg_log2FC)>0.5,
                p_val_adj<0.05)
# Save Sig_Neut_Markers
write.csv(Sig_Neut_Markers,'scRNAseq_Results/Sig_Neut_Markers.csv',row.names = T)
# Top10 highly significant genes in neutrophils during HIV
Top20_Neut_Markers <- Sig_Neut_Markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=20)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("HIV-","HIV+ ART<3 Months","HIV+ ART>1 Year")
Neutrophils <- Seurat::SetIdent(Neutrophils,value = factor(Idents(Neutrophils),
                                                           levels = levels_order))

Top20_Neut_Markers_plot <- Seurat::DotPlot(Neutrophils,
                                      features = unique(Top10_Neut_Markers$gene)#,
                                      #cols = c("blue", "red")
                                      )+
  theme_bw()+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 10,face = 'bold'),
        axis.text.y = element_text(size = 7, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 1))
Top20_Neut_Markers_plot
# Save Figure 4a
ggsave(Top20_Neut_Markers_plot,filename="Figures/Figure 4/Top20_Neut_Markers_plot.png",
       width = 3.5,height = 10,dpi = 1080,units = "in")
# Save Figure 4a
ggsave(Top20_Neut_Markers_plot,filename="Figures/Figure 4/Top20_Neut_Markers_plot.pdf",
       width = 3.5,height = 10,dpi = 1080,units = "in")


# Split DEGs by cluster
DEGs_by_cluster <- split(Sig_Neut_Markers$gene,
                         Sig_Neut_Markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 5)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

GO_Terms_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
GO_Terms_Heatmap
GO_Terms_Heatmap <- as.ggplot(GO_Terms_Heatmap)

# Save Figure 3e
ggsave(GO_Terms_Heatmap,filename="Figures/Figure 4/Fig4d.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(GO_Terms_Heatmap,filename="Figures/Figure 4/Fig4d.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Fig4d2 <- SeuratExtend::Heatmap(heatmap_matrix,
                                          lab_fill='zscore',
                                          plot.margin = margin(l = 30),
                                          y_text_position = "right",
                                          angle = 90,
                                          hjust = 0.5,
                                          color_scheme = color_palette,
                                          vjust = 0.5)
Fig4d2 <- as.ggplot(Fig4d2)
Fig4d2
# Save Figure 4d
ggsave(Fig4d2,filename="Figures/Figure 4/Fig4d2.png",
       width = 5.5,height = 4,dpi = 1080,units = "in")

ggsave(Fig4d2,filename="Figures/Figure 4/Fig4d2.pdf",
       width = 5.5,height = 4,dpi = 1080, units = "in")  


## KEGG pathway and KEGG module analysis
Sig_Neut_Markers <- Neut_Markers %>%
  dplyr::filter(abs(avg_log2FC)>0.5,
                p_val_adj<0.05)

gene_symbols <- unique(Sig_Neut_Markers$gene)

mapped_genes <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
head(mapped_genes)


Sig_Neut_Markers <- merge(
  Sig_Neut_Markers,
  mapped_genes,
  by.x = "gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

head(Sig_Neut_Markers)

Sig_Neut_Markers <- Sig_Neut_Markers %>%
  dplyr::filter(ENTREZID!="NA")

# Split DEGs by cluster
Neut_DEGs_by_cluster_KEGG <- split(Sig_Neut_Markers$ENTREZID,
                                   Sig_Neut_Markers$cluster)


# Perform actual GO analysis
run_KEGG_analysis <- function(ENTREZID, cluster_name){
  enrichKEGG(
    gene = ENTREZID,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05
  )
}
# Perform KEGG analysis for all clusters
Neut_KEGG_results <- lapply(names(Neut_DEGs_by_cluster_KEGG),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_KEGG_analysis(Neut_DEGs_by_cluster_KEGG[[cluster]], cluster)
})

names(Neut_KEGG_results) <- names(Neut_DEGs_by_cluster_KEGG)
subsetted_KEGG_sucategories <- c("Cell growth and death","Cell motility","Cellular community - eukaryotes",
                                 "Immune disease","Immune system","Infectious disease: bacterial","Infectious disease: parasitic",
                                 "Infectious disease: viral","Signal transduction","Signalling molecules and interaction","Transport and catabolism") 
# Extract the @result slot and process data
zscore_data <- lapply(Neut_KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(Neut_KEGG_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(Neut_KEGG_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(Neut_KEGG_results)

# Subset the go_results object
KEGG_results_top20 <- lapply(Neut_KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 50)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
KEGG_results_top20 <- Filter(Negate(is.null), KEGG_results_top20)

# Check the first cluster's top 20 for example
head(KEGG_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(KEGG_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Neut_KEGG_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 12.83862), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = T,
  cluster_columns = F)
Neut_KEGG_Heatmap <- as.ggplot(Neut_KEGG_Heatmap)
Neut_KEGG_Heatmap

# Save Figure 3e
ggsave(Neut_KEGG_Heatmap,filename="Figures/Figure 3/Neut_KEGG_Heatmap.png",
       width = 4.5,height = 12,dpi = 1080,units = "in")

ggsave(Neut_KEGG_Heatmap,filename="Figures/Figure 3/Neut_KEGG_Heatmap.pdf",
       width = 4.5,height = 12,dpi = 1080, units = "in")

# ggplot 2 heatmap
cluster_order <- as.data.frame(colnames(heatmap_matrix))
colnames(cluster_order) <- "Cluster"

feature_order <- as.data.frame(rownames(heatmap_matrix))
colnames(feature_order) <- "features"

heatmap_matrix_new <- as.data.frame(heatmap_matrix)
heatmap_matrix_new$feature <- rownames(heatmap_matrix_new)
heatmap_matrix_new$wrapped_feature <- str_wrap(heatmap_matrix_new$feature, width = 50)

Neut_gg_KEGG_heatmap <- as.data.frame(heatmap_matrix_new) %>%
  pivot_longer(cols = c(1:2),
               names_to = 'Cluster',
               values_to = 'zscore') %>%
  ggplot(aes(x=factor(Cluster,levels = cluster_order$Cluster),
             y=factor(feature,levels = feature_order$features),
             fill=zscore))+
  geom_tile(color='white')+
  scale_y_discrete(labels = heatmap_matrix_new$wrapped_feature)+
  scale_fill_gradient(low = "#FDFFFF",high = "darkred")+
  scale_x_discrete(labels=c("HIV+ ART<3 Months"="PLHIV on ART<3m",
                            "HIV+ ART>1 Year"="PLHIV on ART>1y"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 9),
        axis.title = element_blank(),
        #panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.7,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 5))+
  guides(fill = guide_colorbar(
    title = "zScore",
    title.theme = element_text(size = 9),
    barwidth = 0.8,
    barsize = 1))
Neut_gg_KEGG_heatmap

# Save Figure 3e
ggsave(Neut_gg_KEGG_heatmap,filename="Figures/Figure 3/Neut_gg_KEGG_heatmap.png",
       width = 4.5,height = 11,dpi = 1080,units = "in")

ggsave(Neut_gg_KEGG_heatmap,filename="Figures/Figure 3/Neut_gg_KEGG_heatmap.pdf",
       width = 4.5,height = 11,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Neut_KEGG_heatmap_SE <- SeuratExtend::Heatmap(heatmap_matrix,
                                                lab_fill='zscore',
                                                plot.margin = margin(l = 30),
                                                y_text_position = "right",
                                                angle = 90,
                                                hjust = 0.5,
                                                color_scheme = color_palette,
                                                vjust = 0.5)
Neut_KEGG_heatmap_SE <- as.ggplot(Neut_KEGG_heatmap_SE)
Neut_KEGG_heatmap_SE
# Save Figure 3e
ggsave(Neut_KEGG_heatmap_SE,filename="Figures/Figure 3/Neut_KEGG_heatmap_SE.png",
       width = 5.0,height = 11,dpi = 1080, units = "in")

ggsave(Neut_KEGG_heatmap_SE,filename="Figures/Figure 3/Neut_KEGG_heatmap_SE.pdf",
       width = 5.0,height = 11,dpi = 1080, units = "in")

# Volcano plot of differentially expressed genes in Neutrophils


Sig_Neut_Markers <- Neut_Markers %>%
  dplyr::filter(abs(avg_log2FC)>0.25,p_val_adj<0.05) %>%
  dplyr::filter(pct.1>0.1) %>%
  dplyr::mutate(Diff_Expressed = case_when(
    avg_log2FC >= 1.5 & p_val_adj < 0.1 ~ "Upregulated",
    avg_log2FC < -1.5 & p_val_adj < 0.1 ~ "Downregulated",
    TRUE ~ "Not Significant" # Assign a default value if none of the conditions match
  ))


Sig_Neut_Markers_3M <- Sig_Neut_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  dplyr::select(p_val,avg_log2FC,p_val_adj,gene)
colnames(Sig_Neut_Markers_3M) <- c("pvalue","log2FoldChange","padj","gene")

EnhancedVolcano::EnhancedVolcano(
  Sig_Neut_Markers_3M,
  lab = Sig_Neut_Markers_3M$gene,
  x="log2FoldChange",
  y="pvalue",
  title = "ART<3m vs HIV-",
  pCutoff = 10e-6,
  pointSize = 0.5,
  labSize = 4,
  col=c('black', 'black', 'black', 'red3'),
  colAlpha = 1
)

Neut_Volcano_3M <- Sig_Neut_Markers %>%
  filter(cluster=="HIV+ ART<3 Months") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 15,size=4)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART<3m")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  #scale_color_manual(values = c("red",'blue',"grey"))+
  scale_y_continuous(limits = c(0,30))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Neut_Volcano_3M
# Save Figure 4j
ggsave(Neut_Volcano_3M,filename="Figures/Figure 4/Neut_Volcano_3M.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Neut_Volcano_3M,filename="Figures/Figure 4/Neut_Volcano_3M.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Neut_Volcano_1Y <- Sig_Neut_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=4)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART>1y")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,20))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Neut_Volcano_1Y

# Save Figure 4k
ggsave(Neut_Volcano_1Y,filename="Figures/Figure 4/Neut_Volcano_1Y.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Neut_Volcano_1Y,filename="Figures/Figure 4/Neut_Volcano_1Y.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Neut_Volcano_HIVNeg <- Sig_Neut_Markers %>%
  filter(cluster=="HIV-") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "HIV-")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,30))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Neut_Volcano_HIVNeg  

# Save Figure 4l
ggsave(Neut_Volcano_HIVNeg,filename="Figures/Figure 4/Fig4l.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Neut_Volcano_HIVNeg,filename="Figures/Figure 4/Fig4l.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")


# Enrichment analysis of differentially expressed genes in Neutrophils
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pathview)

Sig_Neut_Markers_3M <- Sig_Neut_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)
  
Sig_Neut_Markers_1Y <- Sig_Neut_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

sig_Neut_Markers_HIVNeg <- Sig_Neut_Markers %>%
  filter(cluster=="HIV-") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

# Perform GO analysis for each condition
# Enrich GO analysis in Neutrophils from PLHIV on ART<3m
Neut_go_results_3M <- enrichGO(
  gene = Sig_Neut_Markers_3M$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Neut_Go_bar_3M <- barplot(Neut_go_results_3M, showCategory = 20)
Neut_Go_bar_3M <- as.ggplot(Neut_Go_bar_3M)
ggsave(Neut_Go_bar_3M,filename="Figures/Figure 4/Fig4lGONeutbar3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Neut_Go_dot_3M <- clusterProfiler::dotplot(Neut_go_results_3M,showCategory=20)
Neut_Go_dot_3M <- as.ggplot(Neut_Go_dot_3M)
ggsave(Neut_Go_dot_3M,filename="Figures/Figure 4/Fig4l2GONeutdot3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from PLHIV on ART>1y
Neut_go_results_1Y <- enrichGO(
  gene = Sig_Neut_Markers_1Y$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Neut_Go_bar_1Y <- barplot(Neut_go_results_1Y, showCategory = 20)
Neut_Go_bar_1Y <- as.ggplot(Neut_Go_bar_1Y)
ggsave(Neut_Go_bar_1Y,filename="Figures/Figure 4/Fig4lGONeutbar1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Neut_Go_dot_1Y <- clusterProfiler::dotplot(Neut_go_results_1Y,showCategory=20)
Neut_Go_dot_1Y <- as.ggplot(Neut_Go_dot_1Y)
ggsave(Neut_Go_dot_1Y,filename="Figures/Figure 4/Fig4l2GONeutdot1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from HIV-
Neut_go_results_HIVneg <- enrichGO(
  gene = sig_Neut_Markers_HIVNeg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Neut_Go_bar_HIVneg <- barplot(Neut_go_results_HIVneg, showCategory = 20)
Neut_Go_bar_HIVneg <- as.ggplot(Neut_Go_bar_HIVneg)
ggsave(Neut_Go_bar_HIVneg,filename="Figures/Figure 4/Fig4lGONeutbarNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Neut_Go_dot_HIVneg <- clusterProfiler::dotplot(Neut_go_results_HIVneg,showCategory=20)
Neut_Go_dot_HIVneg <- as.ggplot(Neut_Go_dot_HIVneg)
ggsave(Neut_Go_dot_HIVneg,filename="Figures/Figure 4/Fig4l2GONeutdotNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# CEMiTool in Neutrophils
# Aggregate expression
Neutrophil_cts <- Seurat::AggregateExpression(
  object = Neutrophils,
  group.by = c("HIV_Status","sample"),
  assays = "RNA",
  slot = "counts",
  return.seurat = F)
# Convert to dataframe
Neutrophil_cts <- as.data.frame(Neutrophil_cts$RNA)
# generate sample level metadata
colData <- data.frame(samples = colnames(Neutrophil_cts))
colData <- colData %>%
  dplyr::mutate(HIV_Status = ifelse(grepl("HIV-", samples), "HIV-", 
                                    ifelse(grepl("ART<3",samples),"PLHIV on ART<3m","PLHIV on ART>1y"))) %>%
  dplyr:: mutate(Sample_Name = str_extract(samples, "(?<=_).*")) %>%
  column_to_rownames(var = "samples")
# Create a DESeq2 object
Neutrophil_dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = Neutrophil_cts,
  colData = colData,
  design = ~ HIV_Status)
# Run DESeq2
Neutrophil_dds <- DESeq2::DESeq(Neutrophil_dds)
# Extract data from the dds object
Neutrophils_out <- DESeq2::rlog(Neutrophil_dds,blind = FALSE)
# Create a matrix
Neutrophils_matrix <- assay(Neutrophil_dds)
Neutrophils_matrix <- na.omit(Neutrophils_matrix)
Neutrophils_matrix <- as.data.frame(Neutrophils_matrix)
# Create a sample annotation file
sample_annot <- colData %>%
  dplyr::mutate(SampleName = rownames(.)) %>%
  dplyr::mutate(Class = ifelse(grepl("HIV-", SampleName), "HIV-", 
                               ifelse(grepl("ART<3",SampleName),"PLHIV on ART<3m","PLHIV on ART>1y"))) %>%
  dplyr::select(SampleName,Class)
write.csv(sample_annot,"scRNAseq_Results/Sample_annotation.csv",row.names = F)
# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)
# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)
# Run CEMiTool
Neutrophil_cem <- CEMiTool::cemitool(Neutrophils_matrix,
                                     annot = sample_annot,
                                     gmt = gmt_in,
                                     interactions = int_df,
                                     filter = TRUE,
                                     filter_pval = 0.05,
                                     #apply_vst = T,
                                     cor_method = "pearson",
                                     cor_function = "cor",
                                     network_type = "signed",
                                     gsea_scale = T,
                                     force_beta = T,
                                     ora_pval = 0.05,
                                     min_ngen = 5,
                                     gsea_min_size = 3,
                                     gsea_max_size = 800,
                                     plot = T,
                                     center_func = 'median',
                                     plot_diagnostics = T)

# Plot gene set enrichment analysis modules
Neutrophil_gsea_plots <- CEMiTool::show_plot(Neutrophil_cem,"gsea")
Neutrophil_gsea_plots

Neutrophil_gsea_plot <- Neutrophil_cem@enrichment$nes %>%
  as.data.frame() %>%
  filter(pathway=="M1" | pathway=="M2" | pathway=="M3" | pathway=="M4") %>%
  pivot_longer(cols = c("HIV-",'PLHIV on ART<3m', 'PLHIV on ART>1y'), 
               names_to = "HIV_Status", values_to = "NES") %>%
  ggplot(aes(HIV_Status, pathway, fill = NES))+
  geom_point(aes(col = NES, size = NES))+
  scale_size(range = c(3, 18), guide = 'none') +
  #scale_color_gradient(low = "#0000B4", high = "#931500")+
  #scale_fill_gradient(low = "#0000B4", high = "#931500")+
  scale_color_gradient(low = "grey", high = "darkblue")+
  scale_fill_gradient(low = "grey", high = "darkblue")+
  labs(title = "Neutrophil Gene Set Enrichment Analysis")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(face = 'bold',hjust = 0.5))+
  guides(
    color = guide_colorbar(
      title = "NES",
      title.theme = element_text(size = 9),
      barwidth = 0.8))
Neutrophil_gsea_plot
ggsave(Neutrophil_gsea_plot,filename="Figures/Figure 4/Neut_gsea_plot.png",
       width = 5.5,height = 3, dpi = 1080,units = "in")

# Plot over representation analysis
Neutrophil_ora_plots <- CEMiTool::show_plot(Neutrophil_cem, "ora")
Neutrophil_ora_plots$M1
Neutrophil_ora_plots$M2
Neutrophil_ora_plots$M3
Neutrophil_ora_plots$M4
Neutrophil_ora_plots$M5

library(forcats)
Neut_M1_ora <- Neutrophil_cem@ora %>%
  dplyr::filter(Module=="M1", p.adjust<0.2) %>%
  sort("p.adjust") %>%
  slice_head(n=10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)), fill = -log10(p.adjust)))+
  geom_col()+
  #scale_fill_gradient(low = "#4A527F", high = "#931500")+
  scale_fill_gradient(low = "grey", high = "darkblue")+
  labs(x = "-log10(adjusted p value)", title = "M1 Gene Set Enrichment Analysis")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.75,0.25),
        #legend.position = "none",
        axis.title.y = element_blank(),
        axis.title = element_text(face = 'bold',size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 1))+
  guides(fill = guide_colorbar(
      title = "-log10 P-value",
      title.theme = element_text(size = 9),
      barwidth = 0.8,
      barsize = 1))
Neut_M1_ora
ggsave(Neut_M1_ora,filename="Figures/Figure 4/Neut_M1_ora.png",
       width = 6,height = 4.5,dpi = 1080,units = "in")

Neut_M2_ora <- Neutrophil_cem@ora %>%
  dplyr::filter(Module=="M2", p.adjust<0.2) %>%
  sort("p.adjust") %>%
  slice_head(n=10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)), fill = -log10(p.adjust)))+
  geom_col()+
  #scale_fill_gradient(low = "#4A527F", high = "#931500")+
  scale_fill_gradient(low = "grey", high = "darkblue")+
  labs(x = "-log10(adjusted p value)", title = "M2 Gene Set Enrichment Analysis")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.75,0.25),
        #legend.position = "none",
        axis.title.y = element_blank(),
        axis.title = element_text(face = 'bold',size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 1))+
  guides(fill = guide_colorbar(
    title = "-log10 P-value",
    title.theme = element_text(size = 9),
    barwidth = 0.8,
    barsize = 1))
Neut_M2_ora
ggsave(Neut_M2_ora,filename="Figures/Figure 4/Neut_M2_ora.png",
       width = 6,height = 4.5,dpi = 1080,units = "in")


# plot interactions
Neutrophil_interactions_plots <- CEMiTool::show_plot(Neutrophil_cem, "interaction") # view the plot for the first module
Neutrophil_interactions_plots$M1
Neutrophil_interactions_plots$M2
Neutrophil_interactions_plots$M3
Neutrophil_interactions_plots$M4

# create report as html document
generate_report(Neutrophil_cem, directory="scRNAseq_Results/Neutrophil_cemitool_Report",force=T)

# write analysis results into files
write_files(Neutrophil_cem, directory="scRNAseq_Results/Neutrophil_cemitool_Report",force=T)

# save all plots
save_plots(Neutrophil_cem, "all", directory="scRNAseq_Results/Neutrophil_cemitool_Report",force=T)

# Gene module scoring of selected pathways in Neutrophils
Neut_geneID <- Neutrophil_cem@ora %>%
  dplyr::filter(Module=="M1") %>%
  dplyr::filter(ID=="Senescence-Associated Secretory Phenotype (SASP)") %>%
  dplyr::pull("geneID")

# Split the string at "/"
Neut_geneID <- unlist(strsplit(Neut_geneID, "/"))
# Print the result
print(Neut_geneID)

# Gene module scoring in Neutrophils
Neutrophils <- AddModuleScore(
  object = Neutrophils,
  features = list(c(Neut_geneID)),
  name = "SASP",
  slot = 'data',
  ctrl = 200)

# Dotplot of SASP in Neutrophils
Seurat::DotPlot(object = Neutrophils,
                features = Neut_geneID)+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1, size = 7,face = 'bold'),
        axis.text.y = element_text(size = 9, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))

Neutrophils@meta.data %>%
  ggplot(aes(HIV_Status,SASP1, fill = HIV_Status))+
  #geom_violin()+
  geom_boxplot(width = 0.1, fill="white",outlier.shape = NA)+
  geom_pwc()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none')
# New volcano plots

```


# FURTHER ANALYSIS OF CD3+ T cells ----
Tcells <- subset(Immune_cells,
                 idents = c("CD8+ T cells","NK T cells"),
                 invert = FALSE)
# Setting HIV Status as idents
Idents(Tcells) <- Tcells$HIV_Status
# Finding Neutrophil markers in
Tcell_Markers <- FindAllMarkers(
  Tcells,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = F,
  min.cells.feature = 3,
  min.cells.group = 3)
# Save Tcells_Markers
write.csv(Tcell_Markers,'scRNAseq_Results/Tcell_Markers.csv',row.names = T)
# Filtering for significant Neut_Markers
Sig_Tcell_Markers <- Tcell_Markers %>%
  dplyr::filter(abs(avg_log2FC)>1.5,
                p_val_adj<0.05)
# Save Sig_Neut_Markers
write.csv(Sig_Tcell_Markers,'scRNAseq_Results/Sig_Tcell_Markers.csv',row.names = T)
# Top10 highly significant genes in neutrophils during HIV
Top10_Tcell_Markers <- Sig_Tcell_Markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=20)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("HIV-","HIV+ ART<3 Months","HIV+ ART>1 Year")
Tcells <- Seurat::SetIdent(Tcells,value = factor(Idents(Tcells),
                                                 levels = levels_order))
Top10_Tcell_Markers_plot <- Seurat::DotPlot(Tcells,
                                           features = unique(Top10_Tcell_Markers$gene),
                                           cols = c("blue", "red"))+
  theme_bw()+
  scale_y_discrete(labels=c("HIV-"="HIV-",
                            "HIV+ ART<3 Months"="ART<3m",
                            "HIV+ ART>1 Year"="ART>1y"))+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 7,face = 'bold',color = HIV_status_colors),
        axis.text.y = element_text(size = 7),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.05, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 0.5))
Top10_Tcell_Markers_plot
# Save Figure 4a
ggsave(Top10_Tcell_Markers_plot,filename="Figures/Figure 4/Top10_Tcell_Markers_plot.png",
       width = 3.5,height = 6,dpi = 1080,units = "in")

ggsave(Top10_Tcell_Markers_plot,filename="Figures/Figure 4/Top10_Tcell_Markers_plot.pdf",
       width = 3.5,height = 6,dpi = 1080,units = "in")

# Split DEGs by cluster
DEGs_by_cluster <- split(Sig_Tcell_Markers$gene,
                         Sig_Tcell_Markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 5)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Fig4f <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Fig4f
Fig4f <- as.ggplot(Fig4f)

# Save Figure 4f
ggsave(Fig4f,filename="Figures/Figure 4/Fig4f.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Fig4f,filename="Figures/Figure 4/Fig4f.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Fig4f2 <- SeuratExtend::Heatmap(heatmap_matrix,
                                lab_fill='zscore',
                                plot.margin = margin(l = 30),
                                y_text_position = "right",
                                angle = 90,
                                hjust = 0.5,
                                color_scheme = color_palette,
                                vjust = 0.5)
Fig4f2 <- as.ggplot(Fig4f2)
Fig4f2
# Save Figure 4d
ggsave(Fig4f2,filename="Figures/Figure 4/Fig4f2.png",
       width = 5.5,height = 5,dpi = 1080,units = "in")

ggsave(Fig4f2,filename="Figures/Figure 4/Fig4f2.pdf",
       width = 5.5,height = 5,dpi = 1080, units = "in")  

## KEGG pathway and KEGG module analysis
Sig_Tcell_Markers <- Tcell_Markers %>%
  dplyr::filter(abs(avg_log2FC)>0.25,
                p_val_adj<0.1)

gene_symbols <- unique(Sig_Tcell_Markers$gene)

mapped_genes <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
head(mapped_genes)


Sig_Tcell_Markers <- merge(
  Sig_Tcell_Markers,
  mapped_genes,
  by.x = "gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

head(Sig_Tcell_Markers)

Sig_Tcell_Markers <- Sig_Tcell_Markers %>%
  dplyr::filter(ENTREZID!="NA")

# Split DEGs by cluster
Tcell_DEGs_by_cluster_KEGG <- split(Sig_Tcell_Markers$ENTREZID,
                                    Sig_Tcell_Markers$cluster)


# Perform actual GO analysis
run_KEGG_analysis <- function(ENTREZID, cluster_name){
  enrichKEGG(
    gene = ENTREZID,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05
  )
}
# Perform KEGG analysis for all clusters
Tcell_KEGG_results <- lapply(names(Tcell_DEGs_by_cluster_KEGG),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_KEGG_analysis(Tcell_DEGs_by_cluster_KEGG[[cluster]], cluster)
})

names(Tcell_KEGG_results) <- names(Tcell_DEGs_by_cluster_KEGG)
subsetted_KEGG_sucategories <- c("Cell growth and death","Cell motility","Cellular community - eukaryotes",
                                 "Immune disease","Immune system","Infectious disease: bacterial","Infectious disease: parasitic",
                                 "Infectious disease: viral","Signal transduction","Signalling molecules and interaction","Transport and catabolism") 
# Extract the @result slot and process data
zscore_data <- lapply(Tcell_KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      #dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(Tcell_KEGG_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(Tcell_KEGG_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(Tcell_KEGG_results)

# Subset the go_results object
KEGG_results_top20 <- lapply(Tcell_KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      #dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 50)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
KEGG_results_top20 <- Filter(Negate(is.null), KEGG_results_top20)

# Check the first cluster's top 20 for example
head(KEGG_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(KEGG_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Tcell_KEGG_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 15.33777), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = T,
  cluster_columns = F)
Tcell_KEGG_Heatmap <- as.ggplot(Tcell_KEGG_Heatmap)
Tcell_KEGG_Heatmap

# Save Figure 3e
ggsave(Tcell_KEGG_Heatmap,filename="Figures/Figure 3/Tcell_KEGG_Heatmap.png",
       width = 4.5,height = 12,dpi = 1080,units = "in")

ggsave(Tcell_KEGG_Heatmap,filename="Figures/Figure 3/Tcell_KEGG_Heatmap.pdf",
       width = 4.5,height = 12,dpi = 1080, units = "in")

# ggplot 2 heatmap
cluster_order <- as.data.frame(colnames(heatmap_matrix))
colnames(cluster_order) <- "Cluster"

feature_order <- as.data.frame(rownames(heatmap_matrix))
colnames(feature_order) <- "features"

heatmap_matrix_new <- as.data.frame(heatmap_matrix)
heatmap_matrix_new$feature <- rownames(heatmap_matrix_new)
heatmap_matrix_new$wrapped_feature <- str_wrap(heatmap_matrix_new$feature, width = 50)

Tcell_gg_KEGG_heatmap <- as.data.frame(heatmap_matrix_new) %>%
  pivot_longer(cols = c(1:3),
               names_to = 'Cluster',
               values_to = 'zscore') %>%
  ggplot(aes(x=factor(Cluster,levels = cluster_order$Cluster),
             y=factor(feature,levels = feature_order$features),
             fill=zscore))+
  geom_tile(color='white')+
  scale_y_discrete(labels = heatmap_matrix_new$wrapped_feature)+
  scale_fill_gradient(low = "#FDFFFF",high = "darkred")+
  scale_x_discrete(labels=c("HIV-"="HIV",
                            "HIV+ ART<3 Months"="PLHIV on ART<3m",
                            "HIV+ ART>1 Year"="PLHIV on ART>1y",
                            ))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 9),
        axis.title = element_blank(),
        #panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.7,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 5))+
  guides(fill = guide_colorbar(
    title = "zScore",
    title.theme = element_text(size = 9),
    barwidth = 0.8,
    barsize = 1))
Tcell_gg_KEGG_heatmap

# Save Figure 3e
ggsave(Tcell_gg_KEGG_heatmap,filename="Figures/Figure 3/Tcell_gg_KEGG_heatmap.png",
       width = 4.5,height = 11,dpi = 1080,units = "in")

ggsave(Tcell_gg_KEGG_heatmap,filename="Figures/Figure 3/Tcell_gg_KEGG_heatmap.pdf",
       width = 4.5,height = 11,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Tcell_KEGG_heatmap_SE <- SeuratExtend::Heatmap(heatmap_matrix,
                                              lab_fill='zscore',
                                              plot.margin = margin(l = 30),
                                              y_text_position = "right",
                                              angle = 90,
                                              hjust = 0.5,
                                              color_scheme = color_palette,
                                              vjust = 0.5)
Tcell_KEGG_heatmap_SE <- as.ggplot(Tcell_KEGG_heatmap_SE)
Tcell_KEGG_heatmap_SE
# Save Figure 3e
ggsave(Tcell_KEGG_heatmap_SE,filename="Figures/Figure 3/Tcell_KEGG_heatmap_SE.png",
       width = 5.0,height = 11,dpi = 1080, units = "in")

ggsave(Tcell_KEGG_heatmap_SE,filename="Figures/Figure 3/Tcell_KEGG_heatmap_SE.pdf",
       width = 5.0,height = 11,dpi = 1080, units = "in")



# Volcano plot of differentially expressed genes in T cells

Sig_Tcell_Markers <- Sig_Tcell_Markers %>%
  dplyr::filter(pct.1>0.05) %>%
  dplyr::mutate(Diff_Expressed = ifelse(avg_log2FC>1.5,"Upregulated","Downregulated"))

Tcell_Volcano_3M <- Sig_Tcell_Markers %>%
  filter(cluster=="HIV+ ART<3 Months") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART<3m")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Tcell_Volcano_3M
# Save Figure 4g
ggsave(Tcell_Volcano_3M,filename="Figures/Figure 4/Fig4g.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Tcell_Volcano_3M,filename="Figures/Figure 4/Fig4g.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Tcell_Volcano_1Y <- Sig_Tcell_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART>1y")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Tcell_Volcano_1Y

# Save Figure 4h
ggsave(Tcell_Volcano_1Y,filename="Figures/Figure 4/Fig4h.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Tcell_Volcano_1Y,filename="Figures/Figure 4/Fig4h.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Tcell_Volcano_HIVNeg <- Sig_Tcell_Markers %>%
  filter(cluster=="HIV-") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "HIV-")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Tcell_Volcano_HIVNeg  

# Save Figure 4i
ggsave(Tcell_Volcano_HIVNeg,filename="Figures/Figure 4/Fig4i.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Tcell_Volcano_HIVNeg,filename="Figures/Figure 4/Fig4i.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

# Enrichment analysis of differentially expressed genes in Neutrophils
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pathview)

Sig_Tcell_Markers_3M <- Sig_Tcell_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Tcell_Markers_1Y <- Sig_Tcell_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

sig_Tcell_Markers_HIVNeg <- Sig_Tcell_Markers %>%
  filter(cluster=="HIV-") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

# Perform GO analysis for each condition
# Enrich GO analysis in Neutrophils from PLHIV on ART<3m
Tcell_go_results_3M <- enrichGO(
  gene = Sig_Tcell_Markers_3M$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Tcell_Go_bar_3M <- barplot(Tcell_go_results_3M, showCategory = 20)
Tcell_Go_bar_3M <- as.ggplot(Tcell_Go_bar_3M)
ggsave(Tcell_Go_bar_3M,filename="Figures/Figure 4/Fig4jGOTcellbar3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Tcell_Go_dot_3M <- clusterProfiler::dotplot(Tcell_go_results_3M,showCategory=20)
Tcell_Go_dot_3M <- as.ggplot(Tcell_Go_dot_3M)
ggsave(Tcell_Go_dot_3M,filename="Figures/Figure 4/Fig4j2GOTcelldot3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from PLHIV on ART>1y
Tcell_go_results_1Y <- enrichGO(
  gene = Sig_Tcell_Markers_1Y$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Tcell_Go_bar_1Y <- barplot(Tcell_go_results_1Y, showCategory = 20)
Tcell_Go_bar_1Y <- as.ggplot(Tcell_Go_bar_1Y)
ggsave(Tcell_Go_bar_1Y,filename="Figures/Figure 4/Fig4kGOTcellbar1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Tcell_Go_dot_1Y <- clusterProfiler::dotplot(Tcell_go_results_1Y,showCategory=20)
Tcell_Go_dot_1Y <- as.ggplot(Tcell_Go_dot_1Y)
ggsave(Tcell_Go_dot_1Y,filename="Figures/Figure 4/Fig4k2GOTcelldot1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from HIV-
Tcell_go_results_HIVneg <- enrichGO(
  gene = sig_Tcell_Markers_HIVNeg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Tcell_Go_bar_HIVneg <- barplot(Tcell_go_results_HIVneg, showCategory = 20)
Tcell_Go_bar_HIVneg <- as.ggplot(Tcell_Go_bar_HIVneg)
ggsave(Tcell_Go_bar_HIVneg,filename="Figures/Figure 4/Fig4lGOTcellbarNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Tcell_Go_dot_HIVneg <- clusterProfiler::dotplot(Tcell_go_results_HIVneg,showCategory=20)
Tcell_Go_dot_HIVneg <- as.ggplot(Tcell_Go_dot_HIVneg)
ggsave(Tcell_Go_dot_HIVneg,filename="Figures/Figure 4/Fig4l2GOTcelldotNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")


# CemiTool in T cells
# Aggregate expression
Tcell_cts <- Seurat::AggregateExpression(
  object = Tcells,
  group.by = c("HIV_Status","sample"),
  assays = "RNA",
  slot = "counts",
  return.seurat = F)
# Convert to dataframe
Tcell_cts <- as.data.frame(Tcell_cts$RNA)
# generate sample level metadata
colData <- data.frame(samples = colnames(Tcell_cts))
colData <- colData %>%
  dplyr::mutate(HIV_Status = ifelse(grepl("HIV-", samples), "HIV-", 
                                    ifelse(grepl("ART<3",samples),"PLHIV on ART<3m","PLHIV on ART>1y"))) %>%
  dplyr:: mutate(Sample_Name = str_extract(samples, "(?<=_).*")) %>%
  column_to_rownames(var = "samples")
# Create a DESeq2 object
Tcell_dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = Tcell_cts,
  colData = colData,
  design = ~ HIV_Status)
# Run DESeq2
Tcell_dds <- DESeq2::DESeq(Tcell_dds)
# Extract data from the dds object
Tcells_out <- DESeq2::rlog(Tcell_dds,blind = FALSE)
# Create a matrix
Tcells_matrix <- rlog(assay(Tcell_dds))
Tcells_matrix <- na.omit(Tcells_matrix)
Tcells_matrix <- as.data.frame(Tcells_matrix)
# Create a sample annotation file
sample_annot <- colData %>%
  dplyr::mutate(SampleName = rownames(.)) %>%
  dplyr::mutate(Class = ifelse(grepl("HIV-", SampleName), "HIV-", 
                               ifelse(grepl("ART<3",SampleName),"PLHIV on ART<3m","PLHIV on ART>1y"))) %>%
  dplyr::select(SampleName,Class)
write.csv(sample_annot,"scRNAseq_Results/Sample_annotation.csv",row.names = F)
# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)
# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)
# Run CEMiTool
Tcell_cem <- CEMiTool::cemitool(Tcells_matrix,
                                annot = sample_annot,
                                gmt = gmt_in,
                                interactions = int_df,
                                filter = TRUE,
                                filter_pval = 0.05,
                                apply_vst = T,
                                cor_method = "pearson",
                                cor_function = "cor",
                                network_type = "signed",
                                gsea_scale = T,
                                force_beta = T,
                                ora_pval = 0.05,
                                min_ngen = 5,
                                gsea_min_size = 3,
                                gsea_max_size = 800,
                                plot = T,
                                center_func = 'median',
                                plot_diagnostics = T)

# Plot gene set enrichment analysis modules
Tcell_gsea_plots <- CEMiTool::show_plot(Tcell_cem,"gsea")
Tcell_gsea_plots

Tcell_gsea_plot <- Tcell_cem@enrichment$nes %>%
  as.data.frame() %>%
  filter(pathway=="M1" | pathway=="M2" | pathway=="M3" | pathway=="M4") %>%
  pivot_longer(cols = c("HIV-",'PLHIV on ART<3m', 'PLHIV on ART>1y'), 
               names_to = "HIV_Status", values_to = "NES") %>%
  ggplot(aes(HIV_Status, pathway, fill = NES))+
  geom_point(aes(col = NES, size = NES))+
  scale_size(range = c(3, 20), guide = 'none') +
  scale_color_gradient(low = "#0000B4", high = "#931500")+
  scale_fill_gradient(low = "#0000B4", high = "#931500")+
  labs(title = "T cell Gene Set Enrichment Analysis")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(face = 'bold',hjust = 0.5))+
  guides(
    color = guide_colorbar(
      title = "NES",
      title.theme = element_text(size = 9),
      barwidth = 0.8))
Tcell_gsea_plot
ggsave(Tcell_gsea_plot,filename="Figures/Figure 4/Tcell_gsea_plot.png",
       width = 5.5,height = 3.0, dpi = 1080,units = "in")

  

# Plot over representation analysis
Tcell_ora_plots <- CEMiTool::show_plot(Tcell_cem, "ora")
Tcell_ora_plots$M1
Tcell_ora_plots$M2
Tcell_ora_plots$M3


library(forcats)
Tcell_M1_ora <- Tcell_cem@ora %>%
  dplyr::filter(Module=="M1", p.adjust<0.2) %>%
  sort("p.adjust") %>%
  slice_head(n=10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)), fill = -log10(p.adjust)))+
  geom_col()+
  scale_fill_gradient(low = "#4A527F", high = "#931500")+
  labs(x = "-log10(adjusted p value)", title = "M1 Gene Set Enrichment Analysis")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #legend.position = c(0.75,0.2),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title = element_text(face = 'bold',size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 1))
Tcell_M1_ora
ggsave(Tcell_M1_ora,filename="Figures/Figure 4/Tcell_M1_ora.png",
       width = 6,height = 4.5,dpi = 1080,units = "in")

Tcell_M2_ora <- Tcell_cem@ora %>%
  dplyr::filter(Module=="M2", p.adjust<0.2) %>%
  sort("p.adjust") %>%
  slice_head(n=10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)), fill = -log10(p.adjust)))+
  geom_col()+
  scale_fill_gradient(low = "#4A527F", high = "#931500")+
  labs(x = "-log10(adjusted p value)", title = "M2 Gene Set Enrichment Analysis")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #legend.position = c(0.75,0.2),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title = element_text(face = 'bold',size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 1))
Tcell_M2_ora
ggsave(Tcell_M2_ora,filename="Figures/Figure 4/Tcell_M2_ora.png",
       width = 6,height = 4.5,dpi = 1080,units = "in")

Tcell_M5_ora <- Tcell_cem@ora %>%
  dplyr::filter(Module=="M5", p.adjust<0.2) %>%
  sort("p.adjust") %>%
  slice_head(n=10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)), fill = -log10(p.adjust)))+
  geom_col()+
  scale_fill_gradient(low = "#4A527F", high = "#931500")+
  labs(x = "-log10(adjusted p value)", title = "M5 Gene Set Enrichment Analysis")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #legend.position = c(0.75,0.2),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title = element_text(face = 'bold',size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 1))
Tcell_M5_ora
ggsave(Tcell_M5_ora,filename="Figures/Figure 4/Tcell_M5_ora.png",
       width = 6,height = 4.5,dpi = 1080,units = "in")


Tcell_interactions_plots <- CEMiTool::show_plot(Tcell_cem, "interaction") # view the plot for the first module
Tcell_interactions_plots$M1
Tcell_interactions_plots$M2
Tcell_interactions_plots$M3
Tcell_interactions_plots$M4
Tcell_interactions_plots$M5

# create report as html document
generate_report(Tcell_cem, directory="scRNAseq_Results/Tcell_cemitool_Report",force=T)

# write analysis results into files
write_files(Tcell_cem, directory="scRNAseq_Results/Tcell_cemitool_Report",force=T)

# save all plots
save_plots(Tcell_cem, "all", directory="scRNAseq_Results/Tcell_cemitool_Report",force=T)

``````
# FURTHER ANALYSIS OF B cells ----
Bcells <- subset(all_merged_subset_labelled,
                 idents = c("B cells"),
                 invert = FALSE)
# Setting HIV Status as idents
Idents(Bcells) <- Bcells$HIV_Status
# Finding Neutrophil markers in
Bcell_Markers <- FindAllMarkers(
  Bcells,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = F,
  min.cells.feature = 3,
  min.cells.group = 3)
# Save Tcells_Markers
write.csv(Bcell_Markers,'scRNAseq_Results/Bcell_Markers.csv',row.names = T)
# Filtering for significant Neut_Markers
Sig_Bcell_Markers <- Bcell_Markers %>%
  dplyr::filter(abs(avg_log2FC)>1.5,
                p_val_adj<0.05)
# Save Sig_Neut_Markers
write.csv(Sig_Bcell_Markers,'scRNAseq_Results/Sig_Bcell_Markers.csv',row.names = T)
# Top10 highly significant genes in neutrophils during HIV
Top10_Bcell_Markers <- Sig_Bcell_Markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=20)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("HIV-","HIV+ ART<3 Months","HIV+ ART>1 Year")
Bcells <- Seurat::SetIdent(Bcells,value = factor(Idents(Bcells),
                                                 levels = levels_order))
Top10_Bcell_Markers_plot <- Seurat::DotPlot(Bcells,
                                            features = unique(Top10_Bcell_Markers$gene),
                                            cols = c("blue", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 7,face = 'bold'),
        axis.text.y = element_text(size = 9, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 1))
Top10_Bcell_Markers_plot
# Save Figure 4a
ggsave(Top10_Bcell_Markers_plot,filename="Figures/Figure 5/Fig5a.png",
       width = 10,height = 3.5,dpi = 1080,units = "in")

ggsave(Top10_Bcell_Markers_plot,filename="Figures/Figure 5/Fig5a.pdf",
       width = 10,height = 3.5,dpi = 1080, units = "in")

# Split DEGs by cluster
DEGs_by_cluster <- split(Sig_Bcell_Markers$gene,
                         Sig_Bcell_Markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 5)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Fig5b <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Fig5b
Fig5b <- as.ggplot(Fig5b)

# Save Figure 4f
ggsave(Fig5b,filename="Figures/Figure 5/Fig5b.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Fig5b,filename="Figures/Figure 5/Fig5b.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Fig5b2 <- SeuratExtend::Heatmap(heatmap_matrix,
                                lab_fill='zscore',
                                plot.margin = margin(l = 30),
                                y_text_position = "right",
                                angle = 90,
                                hjust = 0.5,
                                color_scheme = color_palette,
                                vjust = 0.5)
Fig5b2 <- as.ggplot(Fig5b2)
Fig5b2
# Save Figure 4d
ggsave(Fig5b2,filename="Figures/Figure 5/Fig5b2.png",
       width = 6,height = 5,dpi = 1080, units = "in")

ggsave(Fig5b2,filename="Figures/Figure 5/Fig5b2.pdf",
       width = 6,height = 5,dpi = 1080, units = "in")


# Volcano plot of differentially expressed genes in B cells

Sig_Bcell_Markers <- Sig_Bcell_Markers %>%
  dplyr::filter(pct.1>0.05) %>%
  dplyr::mutate(Diff_Expressed = ifelse(avg_log2FC>1.5,"Upregulated","Downregulated"))

Bcell_Volcano_3M <- Sig_Bcell_Markers %>%
  filter(cluster=="HIV+ ART<3 Months") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART<3m")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Bcell_Volcano_3M
# Save Figure 4g
ggsave(Bcell_Volcano_3M,filename="Figures/Figure 5/Fig5c.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Bcell_Volcano_3M,filename="Figures/Figure 5/Fig5c.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Bcell_Volcano_1Y <- Sig_Bcell_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART>1y")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Bcell_Volcano_1Y

# Save Figure 4h
ggsave(Bcell_Volcano_1Y,filename="Figures/Figure 5/Fig5d.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Bcell_Volcano_1Y,filename="Figures/Figure 5/Fig5d.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Bcell_Volcano_HIVNeg <- Sig_Bcell_Markers %>%
  filter(cluster=="HIV-") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "HIV-")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Bcell_Volcano_HIVNeg  

# Save Figure 4i
ggsave(Bcell_Volcano_HIVNeg,filename="Figures/Figure 5/Fig5e.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Bcell_Volcano_HIVNeg,filename="Figures/Figure 5/Fig5e.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

# Enrichment analysis of differentially expressed genes in B cells
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pathview)

Sig_Bcell_Markers_3M <- Sig_Bcell_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Bcell_Markers_1Y <- Sig_Bcell_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

sig_Bcell_Markers_HIVNeg <- Sig_Bcell_Markers %>%
  filter(cluster=="HIV-") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

# Perform GO analysis for each condition
# Enrich GO analysis in Neutrophils from PLHIV on ART<3m
Bcell_go_results_3M <- enrichGO(
  gene = Sig_Bcell_Markers_3M$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Bcell_Go_bar_3M <- barplot(Bcell_go_results_3M, showCategory = 20)
Bcell_Go_bar_3M <- as.ggplot(Bcell_Go_bar_1Y)
ggsave(Bcell_Go_bar_3M,filename="Figures/Figure 5/Fig5fGOBcellbar3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Bcell_Go_dot_3M <- clusterProfiler::dotplot(Bcell_go_results_3M,showCategory=20)
Bcell_Go_dot_3M <- as.ggplot(Bcell_Go_dot_3M)
ggsave(Bcell_Go_dot_3M,filename="Figures/Figure 5/Fig5f2GOBcelldot3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from PLHIV on ART>1y
Bcell_go_results_1Y <- enrichGO(
  gene = Sig_Bcell_Markers_1Y$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Bcell_Go_bar_1Y <- barplot(Bcell_go_results_1Y, showCategory = 20)
Bcell_Go_bar_1Y <- as.ggplot(Bcell_Go_bar_1Y)
ggsave(Bcell_Go_bar_1Y,filename="Figures/Figure 5/Fig5gGOBcellbar1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Bcell_Go_dot_1Y <- clusterProfiler::dotplot(Bcell_go_results_1Y,showCategory=20)
Bcell_Go_dot_1Y <- as.ggplot(Bcell_Go_dot_1Y)
ggsave(Bcell_Go_dot_1Y,filename="Figures/Figure 5/Fig5g2GOBcelldot1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from HIV-
Bcell_go_results_HIVneg <- enrichGO(
  gene = sig_Bcell_Markers_HIVNeg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Bcell_Go_bar_HIVneg <- barplot(Bcell_go_results_HIVneg, showCategory = 20)
Bcell_Go_bar_HIVneg <- as.ggplot(Bcell_Go_bar_HIVneg)
ggsave(Bcell_Go_bar_HIVneg,filename="Figures/Figure 5/Fig5hGOBcellbarNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Bcell_Go_dot_HIVneg <- clusterProfiler::dotplot(Bcell_go_results_HIVneg,showCategory=20)
Bcell_Go_dot_HIVneg <- as.ggplot(Bcell_Go_dot_HIVneg)
ggsave(Bcell_Go_dot_HIVneg,filename="Figures/Figure 5/Fig5h2GOBcelldotNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")


# CemiTool in T cells
# Aggregate expression
Bcell_cts <- Seurat::AggregateExpression(
  object = Bcells,
  group.by = c("HIV_Status","sample"),
  assays = "RNA",
  slot = "counts",
  return.seurat = F)
# Convert to dataframe
Bcell_cts <- as.data.frame(Bcell_cts$RNA)
# generate sample level metadata
colData <- data.frame(samples = colnames(Bcell_cts))
colData <- colData %>%
  dplyr::mutate(HIV_Status = ifelse(grepl("HIV-", samples), "HIV-", 
                                    ifelse(grepl("ART<3",samples),"PLHIV on ART<3m","PLHIV on ART>1y"))) %>%
  dplyr:: mutate(Sample_Name = str_extract(samples, "(?<=_).*")) %>%
  column_to_rownames(var = "samples")
# Create a DESeq2 object
Bcell_dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = Bcell_cts,
  colData = colData,
  design = ~ HIV_Status)
# Run DESeq2
Bcell_dds <- DESeq2::DESeq(Bcell_dds)
# Extract data from the dds object
Bcells_out <- DESeq2::rlog(Bcell_dds,blind = FALSE)
# Create a matrix
Bcells_matrix <- rlog(assay(Bcell_dds))
Bcells_matrix <- na.omit(Bcells_matrix)
Bcells_matrix <- as.data.frame(Bcells_matrix)
# Create a sample annotation file
sample_annot <- colData %>%
  dplyr::mutate(SampleName = rownames(.)) %>%
  dplyr::mutate(Class = ifelse(grepl("HIV-", SampleName), "HIV-", 
                               ifelse(grepl("ART<3",SampleName),"PLHIV on ART<3m","PLHIV on ART>1y"))) %>%
  dplyr::select(SampleName,Class)
write.csv(sample_annot,"scRNAseq_Results/Sample_annotation.csv",row.names = F)
# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)
# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)
# Run CEMiTool
Bcell_cem <- CEMiTool::cemitool(Bcells_matrix,
                                annot = sample_annot,
                                gmt = gmt_in,
                                interactions = int_df,
                                filter = TRUE,
                                filter_pval = 0.05,
                                #apply_vst = T,
                                cor_method = "pearson",
                                cor_function = "cor",
                                network_type = "signed",
                                gsea_scale = T,
                                force_beta = T,
                                ora_pval = 0.05,
                                min_ngen = 5,
                                gsea_min_size = 3,
                                gsea_max_size = 800,
                                plot = T,
                                center_func = 'median',
                                plot_diagnostics = T)

# Plot gene set enrichment analysis modules
Bcell_gsea_plots <- CEMiTool::show_plot(Bcell_cem,"gsea")
Bcell_gsea_plots


# Plot over representation analysis
Bcell_ora_plots <- CEMiTool::show_plot(Bcell_cem, "ora")
Bcell_ora_plots$M1
Bcell_ora_plots$M2
Bcell_ora_plots$M3


library(forcats)
Bcell_M1_ora <- Bcell_cem@ora %>%
  dplyr::filter(Module=="M1", p.adjust<0.2) %>%
  sort("p.adjust") %>%
  slice_head(n=10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)), fill = -log10(p.adjust)))+
  geom_col()+
  scale_fill_gradient(low = "#4A527F", high = "#931500")+
  labs(x = "-log10(adjusted p value)", title = "M1 Gene Set Enrichment Analysis")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #legend.position = c(0.75,0.2),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title = element_text(face = 'bold',size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 1))
Bcell_M1_ora
ggsave(Bcell_M1_ora,filename="Figures/Figure 4/Bcell_M1_ora.png",
       width = 6,height = 4.5,dpi = 1080,units = "in")

Bcell_M2_ora <- Bcell_cem@ora %>%
  dplyr::filter(Module=="M2", p.adjust<0.2) %>%
  sort("p.adjust") %>%
  slice_head(n=10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)), fill = -log10(p.adjust)))+
  geom_col()+
  scale_fill_gradient(low = "#4A527F", high = "#931500")+
  labs(x = "-log10(adjusted p value)", title = "M2 Gene Set Enrichment Analysis")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #legend.position = c(0.75,0.2),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title = element_text(face = 'bold',size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 1))
Bcell_M2_ora
ggsave(Bcell_M2_ora,filename="Figures/Figure 4/Bcell_M2_ora.png",
       width = 6,height = 4.5,dpi = 1080,units = "in")

Bcell_M3_ora <- Bcell_cem@ora %>%
  dplyr::filter(Module=="M3", p.adjust<0.2) %>%
  sort("p.adjust") %>%
  slice_head(n=10) %>%
  ggplot(aes(-log10(p.adjust),fct_reorder(ID,-log10(p.adjust)), fill = -log10(p.adjust)))+
  geom_col()+
  scale_fill_gradient(low = "#4A527F", high = "#931500")+
  labs(x = "-log10(adjusted p value)", title = "M5 Gene Set Enrichment Analysis")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        #legend.position = c(0.75,0.2),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title = element_text(face = 'bold',size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 1))
Bcell_M3_ora
ggsave(Tcell_M5_ora,filename="Figures/Figure 4/Bcell_M3_ora.png",
       width = 6,height = 4.5,dpi = 1080,units = "in")


Bcell_interactions_plots <- CEMiTool::show_plot(Bcell_cem, "interaction") # view the plot for the first module
Bcell_interactions_plots$M1
Bcell_interactions_plots$M2
Bcell_interactions_plots$M3
Bcell_interactions_plots$M4


# create report as html document
generate_report(Bcell_cem, directory="scRNAseq_Results/Bcell_cemitool_Report",force=T)

# write analysis results into files
write_files(Bcell_cem, directory="scRNAseq_Results/Bcell_cemitool_Report",force=T)

# save all plots
save_plots(Bcell_cem, "all", directory="scRNAseq_Results/Bcell_cemitool_Report",force=T)

```
# FURTHER ANALYSIS OF Monocytes ----
Mono <- subset(all_merged_subset_labelled,
                 idents = c("Mono/Mac"),
                 invert = FALSE)
# Setting HIV Status as idents
Idents(Mono) <- Mono$HIV_Status
# Finding Neutrophil markers in
Mono_Markers <- FindAllMarkers(
  Mono,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = F,
  min.cells.feature = 3,
  min.cells.group = 3)
# Save Tcells_Markers
write.csv(Mono_Markers,'scRNAseq_Results/Mono_Markers.csv',row.names = T)
# Filtering for significant Neut_Markers
Sig_Mono_Markers <- Mono_Markers %>%
  dplyr::filter(abs(avg_log2FC)>1.5,
                p_val_adj<0.05)
# Save Sig_Neut_Markers
write.csv(Sig_Mono_Markers,'scRNAseq_Results/Sig_Mono_Markers.csv',row.names = T)
# Top10 highly significant genes in neutrophils during HIV
Top10_Mono_Markers <- Sig_Mono_Markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=20)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("HIV-","HIV+ ART<3 Months","HIV+ ART>1 Year")
Mono <- Seurat::SetIdent(Mono,value = factor(Idents(Mono),
                                                 levels = levels_order))
Top10_Mono_Markers_plot <- Seurat::DotPlot(Mono,
                                            features = unique(Top10_Mono_Markers$gene),
                                            cols = c("blue", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 7,face = 'bold'),
        axis.text.y = element_text(size = 9, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 1))
Top10_Mono_Markers_plot
# Save Figure 4a
ggsave(Top10_Mono_Markers_plot,filename="Figures/Figure 6/Fig6a.png",
       width = 10,height = 3.5,dpi = 1080,units = "in")

ggsave(Top10_Mono_Markers_plot,filename="Figures/Figure 6/Fig6a.pdf",
       width = 10,height = 3.5,dpi = 1080, units = "in")

# Split DEGs by cluster
DEGs_by_cluster <- split(Sig_Mono_Markers$gene,
                         Sig_Mono_Markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 5)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Fig6b <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Fig6b
Fig6b <- as.ggplot(Fig6b)

# Save Figure 4f
ggsave(Fig6b,filename="Figures/Figure 6/Fig6b.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Fig6b,filename="Figures/Figure 5/Fig6b.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Fig6b2 <- SeuratExtend::Heatmap(heatmap_matrix,
                                lab_fill='zscore',
                                plot.margin = margin(l = 30),
                                y_text_position = "right",
                                angle = 90,
                                hjust = 0.5,
                                color_scheme = color_palette,
                                vjust = 0.5)
Fig6b2 <- as.ggplot(Fig6b2)
Fig6b2
# Save Figure 4d
ggsave(Fig6b2,filename="Figures/Figure 6/Fig6b2.png",
       width = 7,height = 5,dpi = 1080, units = "in")

ggsave(Fig6b2,filename="Figures/Figure 6/Fig6b2.pdf",
       width = 7,height = 5,dpi = 1080, units = "in")


# Volcano plot of differentially expressed genes in B cells

Sig_Mono_Markers <- Sig_Mono_Markers %>%
  dplyr::filter(pct.1>0.05) %>%
  dplyr::mutate(Diff_Expressed = ifelse(avg_log2FC>1.5,"Upregulated","Downregulated"))

Mono_Volcano_3M <- Sig_Mono_Markers %>%
  filter(cluster=="HIV+ ART<3 Months") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART<3m")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Mono_Volcano_3M
# Save Figure 4g
ggsave(Mono_Volcano_3M,filename="Figures/Figure 6/Fig6c.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Mono_Volcano_3M,filename="Figures/Figure 6/Fig6c.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Mono_Volcano_1Y <- Sig_Mono_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART>1y")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Mono_Volcano_1Y

# Save Figure 4h
ggsave(Mono_Volcano_1Y,filename="Figures/Figure 6/Fig6d.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Mono_Volcano_1Y,filename="Figures/Figure 6/Fig6d.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Mono_Volcano_HIVNeg <- Sig_Mono_Markers %>%
  filter(cluster=="HIV-") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "HIV-")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,25))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Mono_Volcano_HIVNeg  

# Save Figure 4i
ggsave(Mono_Volcano_HIVNeg,filename="Figures/Figure 6/Fig6e.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Mono_Volcano_HIVNeg,filename="Figures/Figure 6/Fig6e.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

# Enrichment analysis of differentially expressed genes in B cells
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pathview)

Sig_Mono_Markers_3M <- Sig_Mono_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Mono_Markers_1Y <- Sig_Mono_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

sig_Mono_Markers_HIVNeg <- Sig_Mono_Markers %>%
  filter(cluster=="HIV-") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

# Perform GO analysis for each condition
# Enrich GO analysis in Neutrophils from PLHIV on ART<3m
Mono_go_results_3M <- enrichGO(
  gene = Sig_Mono_Markers_3M$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Mono_Go_bar_3M <- barplot(Mono_go_results_3M, showCategory = 20)
Mono_Go_bar_3M <- as.ggplot(Mono_Go_bar_3M)
ggsave(Mono_Go_bar_3M,filename="Figures/Figure 6/Fig6fGOMonobar3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Mono_Go_dot_3M <- clusterProfiler::dotplot(Mono_go_results_3M,showCategory=20)
Mono_Go_dot_3M <- as.ggplot(Mono_Go_dot_3M)
ggsave(Mono_Go_dot_3M,filename="Figures/Figure 6/Fig6f2GOMonodot3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from PLHIV on ART>1y
Mono_go_results_1Y <- enrichGO(
  gene = Sig_Mono_Markers_1Y$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Mono_Go_bar_1Y <- barplot(Mono_go_results_1Y, showCategory = 20)
Mono_Go_bar_1Y <- as.ggplot(Mono_Go_bar_1Y)
ggsave(Mono_Go_bar_1Y,filename="Figures/Figure 6/Fig6gGOMonobar1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Mono_Go_dot_1Y <- clusterProfiler::dotplot(Mono_go_results_1Y,showCategory=20)
Mono_Go_dot_1Y <- as.ggplot(Mono_Go_dot_1Y)
ggsave(Mono_Go_dot_1Y,filename="Figures/Figure 6/Fig6g2GOMonodot1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from HIV-
Mono_go_results_HIVneg <- enrichGO(
  gene = sig_Mono_Markers_HIVNeg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Mono_Go_bar_HIVneg <- barplot(Mono_go_results_HIVneg, showCategory = 20)
Mono_Go_bar_HIVneg <- as.ggplot(Mono_Go_bar_1Y)
ggsave(Mono_Go_bar_HIVneg,filename="Figures/Figure 6/Fig6hGOMonobarNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Mono_Go_dot_HIVneg <- clusterProfiler::dotplot(Mono_go_results_HIVneg,showCategory=20)
Mono_Go_dot_HIVneg <- as.ggplot(Mono_Go_dot_HIVneg)
ggsave(Mono_Go_dot_HIVneg,filename="Figures/Figure 6/Fig6h2GOMonodotNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
```
# FURTHER ANALYSIS OF Dendritic cells ----
DC <- subset(all_merged_subset_labelled,
               idents = c("Dendritic cells"),
               invert = FALSE)
# Setting HIV Status as idents
Idents(DC) <- DC$HIV_Status
# Finding Neutrophil markers in
DC_Markers <- FindAllMarkers(
  DC,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = F,
  min.cells.feature = 3,
  min.cells.group = 3)
# Save Tcells_Markers
write.csv(DC_Markers,'scRNAseq_Results/DC_Markers.csv',row.names = T)
# Filtering for significant Neut_Markers
Sig_DC_Markers <- DC_Markers %>%
  dplyr::filter(abs(avg_log2FC)>1.5,
                p_val_adj<0.05)
# Save Sig_Neut_Markers
write.csv(Sig_DC_Markers,'scRNAseq_Results/Sig_DC_Markers.csv',row.names = T)
# Top10 highly significant genes in neutrophils during HIV
Top10_DC_Markers <- Sig_DC_Markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=20)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("HIV-","HIV+ ART<3 Months","HIV+ ART>1 Year")
DC <- Seurat::SetIdent(DC,value = factor(Idents(DC),
                                         levels = levels_order))
Top10_DC_Markers_plot <- Seurat::DotPlot(DC,
                                         features = unique(Top10_DC_Markers$gene),
                                         cols = c("blue", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 7,face = 'bold'),
        axis.text.y = element_text(size = 9, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 1))
Top10_DC_Markers_plot
# Save Figure 4a
ggsave(Top10_DC_Markers_plot,filename="Figures/Figure 7/Fig7a.png",
       width = 10,height = 3.5,dpi = 1080,units = "in")

ggsave(Top10_DC_Markers_plot,filename="Figures/Figure 7/Fig7a.pdf",
       width = 10,height = 3.5,dpi = 1080, units = "in")

# Split DEGs by cluster
DEGs_by_cluster <- split(Sig_DC_Markers$gene,
                         Sig_DC_Markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 5)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Fig7b <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Fig7b
Fig7b <- as.ggplot(Fig7b)

# Save Figure 4f
ggsave(Fig7b,filename="Figures/Figure 7/Fig7b.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Fig7b,filename="Figures/Figure 7/Fig7b.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Fig7b2 <- SeuratExtend::Heatmap(heatmap_matrix,
                                lab_fill='zscore',
                                plot.margin = margin(l = 30),
                                y_text_position = "right",
                                angle = 90,
                                hjust = 0.5,
                                color_scheme = color_palette,
                                vjust = 0.5)
Fig7b2 <- as.ggplot(Fig6b2)
Fig7b2
# Save Figure 4d
ggsave(Fig7b2,filename="Figures/Figure 7/Fig7b2.png",
       width = 7,height = 5,dpi = 1080, units = "in")

ggsave(Fig7b2,filename="Figures/Figure 7/Fig7b2.pdf",
       width = 7,height = 5,dpi = 1080, units = "in")


# Volcano plot of differentially expressed genes in B cells

Sig_DC_Markers <- Sig_DC_Markers %>%
  dplyr::filter(pct.1>0.05) %>%
  dplyr::mutate(Diff_Expressed = ifelse(avg_log2FC>1.5,"Upregulated","Downregulated"))

DC_Volcano_3M <- Sig_DC_Markers %>%
  filter(cluster=="HIV+ ART<3 Months") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART<3m")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
DC_Volcano_3M
# Save Figure 4g
ggsave(DC_Volcano_3M,filename="Figures/Figure 7/Fig7c.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(DC_Volcano_3M,filename="Figures/Figure 7/Fig7c.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

DC_Volcano_1Y <- Sig_DC_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART>1y")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
DC_Volcano_1Y

# Save Figure 4h
ggsave(DC_Volcano_1Y,filename="Figures/Figure 7/Fig7d.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(DC_Volcano_1Y,filename="Figures/Figure 7/Fig7d.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

DC_Volcano_HIVNeg <- Sig_DC_Markers %>%
  filter(cluster=="HIV-") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "HIV-")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,25))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
DC_Volcano_HIVNeg  

# Save Figure 4i
ggsave(DC_Volcano_HIVNeg,filename="Figures/Figure 7/Fig7e.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(DC_Volcano_HIVNeg,filename="Figures/Figure 7/Fig7e.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

# Enrichment analysis of differentially expressed genes in DC
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pathview)

Sig_DC_Markers_3M <- Sig_DC_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_DC_Markers_1Y <- Sig_DC_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_DC_Markers_HIVNeg <- Sig_DC_Markers %>%
  filter(cluster=="HIV-") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

# Perform GO analysis for each condition
# Enrich GO analysis in Neutrophils from PLHIV on ART<3m
DC_go_results_3M <- enrichGO(
  gene = Sig_DC_Markers_3M$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

DC_Go_bar_3M <- barplot(DC_go_results_3M, showCategory = 20)
DC_Go_bar_3M <- as.ggplot(DC_Go_bar_3M)
ggsave(DC_Go_bar_3M,filename="Figures/Figure 7/Fig7fGODCbar3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")
DC_Go_dot_3M <- clusterProfiler::dotplot(DC_go_results_3M,showCategory=20)
DC_Go_dot_3M <- as.ggplot(DC_Go_dot_3M)
ggsave(DC_Go_dot_3M,filename="Figures/Figure 7/Fig7f2GODCdot3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from PLHIV on ART>1y
DC_go_results_1Y <- enrichGO(
  gene = Sig_DC_Markers_1Y$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

DC_Go_bar_1Y <- barplot(DC_go_results_1Y, showCategory = 20)
DC_Go_bar_1Y <- as.ggplot(DC_Go_bar_1Y)
ggsave(DC_Go_bar_1Y,filename="Figures/Figure 7/Fig7gGODCbar1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")
DC_Go_dot_1Y <- clusterProfiler::dotplot(DC_go_results_1Y,showCategory=20)
DC_Go_dot_1Y <- as.ggplot(DC_Go_dot_1Y)
ggsave(DC_Go_dot_1Y,filename="Figures/Figure 7/Fig7g2GODCdot1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from HIV-
DC_go_results_HIVneg <- enrichGO(
  gene = Sig_DC_Markers_HIVNeg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

DC_Go_bar_HIVneg <- barplot(DC_go_results_HIVneg, showCategory = 20)
DC_Go_bar_HIVneg <- as.ggplot(DC_Go_bar_HIVneg)
ggsave(DC_Go_bar_HIVneg,filename="Figures/Figure 7/Fig7hGODCbarNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
DC_Go_dot_HIVneg <- clusterProfiler::dotplot(DC_go_results_HIVneg,showCategory=20)
DC_Go_dot_HIVneg <- as.ggplot(DC_Go_dot_HIVneg)
ggsave(DC_Go_dot_HIVneg,filename="Figures/Figure 7/Fig7h2GODCdotNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# FURTHER ANALYSIS OF Basal cells ----
all_merged_subset_labelled_new <- readRDS("Data/Single_Cell_Data/all_merged_subset_labelled_new.rds")
Basal <- subset(all_merged_subset_labelled_new,
             idents = c("Basal cells"),
             invert = FALSE)
# Setting HIV Status as idents
Idents(Basal) <- Basal$HIV_Status
# Finding Neutrophil markers in
Basal_Markers <- FindAllMarkers(
  Basal,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = F,
  min.cells.feature = 3,
  min.cells.group = 3)
# Save Tcells_Markers
write.csv(Basal_Markers,'scRNAseq_Results/Basal_Markers.csv',row.names = T)
# Filtering for significant Neut_Markers
Sig_Basal_Markers <- Basal_Markers %>%
  dplyr::filter(abs(avg_log2FC)>0.50,
                p_val_adj<0.05)
# Save Sig_Neut_Markers
write.csv(Sig_Basal_Markers,'scRNAseq_Results/Sig_Basal_Markers.csv',row.names = T)
# Top10 highly significant genes in neutrophils during HIV
Top10_Basal_Markers <- Sig_Basal_Markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=20)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("HIV-","HIV+ ART<3 Months","HIV+ ART>1 Year")
Basal <- Seurat::SetIdent(Basal,value = factor(Idents(Basal),
                                         levels = levels_order))
Top10_Basal_Markers_plot <- Seurat::DotPlot(Basal,
                                            features = unique(Top10_Basal_Markers$gene),
                                            cols = c("blue", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 7,face = 'bold'),
        axis.text.y = element_text(size = 9, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 1))
Top10_Basal_Markers_plot
# Save Figure 4a
ggsave(Top10_Basal_Markers_plot,filename="Figures/Figure 8/Fig8a.png",
       width = 10,height = 3.5,dpi = 1080,units = "in")

ggsave(Top10_Basal_Markers_plot,filename="Figures/Figure 8/Fig8a.pdf",
       width = 10,height = 3.5,dpi = 1080, units = "in")

# Split DEGs by cluster
DEGs_by_cluster <- split(Sig_Basal_Markers$gene,
                         Sig_Basal_Markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 10)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Fig8b <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Fig8b
Fig8b <- as.ggplot(Fig8b)

# Save Figure 4f
ggsave(Fig8b,filename="Figures/Figure 8/Fig8b.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Fig8b,filename="Figures/Figure 8/Fig8b.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Fig8b2 <- SeuratExtend::Heatmap(heatmap_matrix,
                                lab_fill='zscore',
                                plot.margin = margin(l = 30),
                                y_text_position = "right",
                                angle = 90,
                                hjust = 0.5,
                                color_scheme = color_palette,
                                vjust = 0.5)
Fig8b2 <- as.ggplot(Fig8b2)
Fig8b2
# Save Figure 4d
ggsave(Fig8b2,filename="Figures/Figure 8/Fig8b2.png",
       width = 7,height = 5,dpi = 1080, units = "in")

ggsave(Fig8b2,filename="Figures/Figure 8/Fig8b2.pdf",
       width = 7,height = 5,dpi = 1080, units = "in")


# Volcano plot of differentially expressed genes in Basal cells

Sig_Basal_Markers <- Sig_Basal_Markers %>%
  dplyr::filter(pct.1>0.05) %>%
  dplyr::mutate(Diff_Expressed = ifelse(avg_log2FC>1.5,"Upregulated","Downregulated"))

Basal_Volcano_3M <- Sig_Basal_Markers %>%
  filter(cluster=="HIV+ ART<3 Months") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART<3m")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Basal_Volcano_3M
# Save Figure 4g
ggsave(Basal_Volcano_3M,filename="Figures/Figure 8/Fig8c.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Basal_Volcano_3M,filename="Figures/Figure 8/Fig8c.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Basal_Volcano_1Y <- Sig_Basal_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART>1y")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Basal_Volcano_1Y

# Save Figure 4h
ggsave(Basal_Volcano_1Y,filename="Figures/Figure 8/Fig8d.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Basal_Volcano_1Y,filename="Figures/Figure 8/Fig8d.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Basal_Volcano_HIVNeg <- Sig_Basal_Markers %>%
  filter(cluster=="HIV-") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "HIV-")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,25))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Basal_Volcano_HIVNeg  

# Save Figure 4i
ggsave(Basal_Volcano_HIVNeg,filename="Figures/Figure 8/Fig8e.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Basal_Volcano_HIVNeg,filename="Figures/Figure 8/Fig8e.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

# Enrichment analysis of differentially expressed genes in Basal cells
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pathview)

Sig_Basal_Markers_3M <- Sig_Basal_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Basal_Markers_1Y <- Sig_Basal_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Basal_Markers_HIVNeg <- Sig_Basal_Markers %>%
  filter(cluster=="HIV-") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

# Perform GO analysis for each condition
# Enrich GO analysis in Neutrophils from PLHIV on ART<3m
Basal_go_results_3M <- enrichGO(
  gene = Sig_Basal_Markers_3M$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Basal_Go_bar_3M <- barplot(Basal_go_results_3M, showCategory = 20)
Basal_Go_bar_3M <- as.ggplot(Basal_Go_bar_3M)
print(Basal_Go_bar_3M)
ggsave(Basal_Go_bar_3M,filename="Figures/Figure 8/Fig8fGOBasalbar3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Basal_Go_dot_3M <- clusterProfiler::dotplot(Basal_go_results_3M,showCategory=20)
Basal_Go_dot_3M <- as.ggplot(Basal_Go_dot_3M)
ggsave(Basal_Go_dot_3M,filename="Figures/Figure 8/Fig8f2GOBasaldot3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")


# Enrich GO analysis in Neutrophils from PLHIV on ART>1y
Basal_go_results_1Y <- enrichGO(
  gene = Sig_Basal_Markers_1Y$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Basal_Go_bar_1Y <- barplot(Basal_go_results_1Y, showCategory = 20)
Basal_Go_bar_1Y <- as.ggplot(Basal_Go_bar_1Y)
ggsave(Basal_Go_bar_1Y,filename="Figures/Figure 8/Fig8gGOBasalbar1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Basal_Go_dot_1Y <- clusterProfiler::dotplot(Basal_go_results_1Y,showCategory=20)
Basal_Go_dot_1Y <- as.ggplot(Basal_Go_dot_1Y)
ggsave(Basal_Go_dot_1Y,filename="Figures/Figure 8/Fig8g2GOBasaldot1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from HIV-
Basal_go_results_HIVneg <- enrichGO(
  gene = Sig_Basal_Markers_HIVNeg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Basal_Go_bar_HIVneg <- barplot(Basal_go_results_HIVneg, showCategory = 20)
Basal_Go_bar_HIVneg <- as.ggplot(Basal_Go_bar_HIVneg)
ggsave(Basal_Go_bar_HIVneg,filename="Figures/Figure 8/Fig8hGOBasalbarNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Basal_Go_dot_HIVneg <- clusterProfiler::dotplot(Basal_go_results_HIVneg,showCategory=20)
Basal_Go_dot_HIVneg <- as.ggplot(Basal_Go_dot_HIVneg)
ggsave(Basal_Go_dot_HIVneg,filename="Figures/Figure 8/Fig8h2GOBasaldotNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")


# KEGG Pathways
## KEGG pathway and KEGG module analysis
Basal_sig_markers <- Basal_Markers %>%
  dplyr::filter(p_val_adj <= 0.05,abs(avg_log2FC) > 0.5)

gene_symbols <- unique(Basal_sig_markers$gene)

mapped_genes <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
head(mapped_genes)


Basal_sig_markers <- merge(
  Basal_sig_markers,
  mapped_genes,
  by.x = "gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

head(Basal_sig_markers)

# Split DEGs by cluster
Basal_DEGs_by_cluster_KEGG <- split(Basal_sig_markers$ENTREZID,
                                    Basal_sig_markers$cluster)


# Perform actual GO analysis
run_KEGG_analysis <- function(ENTREZID, cluster_name){
  enrichKEGG(
    gene = ENTREZID,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 0.05
  )
}
# Perform KEGG analysis for all clusters
Basal_KEGG_results <- lapply(names(Basal_DEGs_by_cluster_KEGG),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_KEGG_analysis(Basal_DEGs_by_cluster_KEGG[[cluster]], cluster)
})

names(Basal_KEGG_results) <- names(Basal_DEGs_by_cluster_KEGG)
subsetted_KEGG_sucategories <- c("Cell growth and death","Cell motility","Cellular community - eukaryotes",
                                 "Immune disease","Immune system","Infectious disease: bacterial","Infectious disease: parasitic",
                                 "Infectious disease: viral","Signal transduction","Signalling molecules and interaction","Transport and catabolism") 
# Extract the @result slot and process data
zscore_data <- lapply(Basal_KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.05) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(Basal_KEGG_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(Basal_KEGG_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(Basal_KEGG_results)

# Subset the Kegg_results object
KEGG_results_top20 <- lapply(Basal_KEGG_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(subcategory %in% subsetted_KEGG_sucategories) %>%
      dplyr::filter(p.adjust<0.05) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 30)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
KEGG_results_top20 <- Filter(Negate(is.null), KEGG_results_top20)

# Check the first cluster's top 20 for example
head(KEGG_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(KEGG_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Basal_KEGG_Heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Basal_KEGG_Heatmap <- as.ggplot(Basal_KEGG_Heatmap)
Basal_KEGG_Heatmap

# Save Figure 3e
ggsave(Basal_KEGG_Heatmap,filename="Figures/Figure 3/Basal_KEGG_Heatmap.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Basal_KEGG_Heatmap,filename="Figures/Figure 3/Basal_KEGG_Heatmap.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# ggplot 2 heatmap
cluster_order <- as.data.frame(colnames(heatmap_matrix))
colnames(cluster_order) <- "Cluster"

feature_order <- as.data.frame(rownames(heatmap_matrix))
colnames(feature_order) <- "features"

heatmap_matrix_new <- as.data.frame(heatmap_matrix)
heatmap_matrix_new$feature <- rownames(heatmap_matrix_new)
heatmap_matrix_new$wrapped_feature <- str_wrap(heatmap_matrix_new$feature, width = 50)

Basal_gg_KEGG_heatmap <- as.data.frame(heatmap_matrix_new) %>%
  pivot_longer(cols = c(1:3),
               names_to = 'Cluster',
               values_to = 'zscore') %>%
  ggplot(aes(x=factor(Cluster,levels = cluster_order$Cluster),
             y=factor(feature,levels = feature_order$features),
             fill=zscore))+
  geom_tile(color='white')+
  scale_y_discrete(labels = heatmap_matrix_new$wrapped_feature)+
  scale_fill_gradient(low = "#FDFFFF",high = "darkred")+
  #scale_x_discrete(cluster_order)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 9),
        axis.title = element_blank(),
        #panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.7,'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 5))+
  guides(fill = guide_colorbar(
    title = "zScore",
    title.theme = element_text(size = 9),
    barwidth = 0.8,
    barsize = 1))
Basal_gg_KEGG_heatmap

# Save Figure 3e
ggsave(Basal_gg_KEGG_heatmap,filename="Figures/Figure 3/Basal_gg_KEGG_heatmap.png",
       width = 6,height = 11,dpi = 1080,units = "in")

ggsave(Basal_gg_KEGG_heatmap,filename="Figures/Figure 3/Basal_gg_KEGG_heatmap.pdf",
       width = 6,height = 11,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Basal_KEGG_heatmap_SE <- SeuratExtend::Heatmap(heatmap_matrix,
                                                lab_fill='zscore',
                                                plot.margin = margin(l = 30),
                                                y_text_position = "right",
                                                angle = 90,
                                                hjust = 0.5,
                                                color_scheme = color_palette,
                                                vjust = 0.5)
Basal_KEGG_heatmap_SE <- as.ggplot(Basal_KEGG_heatmap_SE)
Basal_KEGG_heatmap_SE
# Save Figure 3e
ggsave(Basal_KEGG_heatmap_SE,filename="Figures/Figure 3/Basal_KEGG_heatmap_SE.png",
       width = 8,height = 16,dpi = 1080,units = "in")

ggsave(Basal_KEGG_heatmap_SE,filename="Figures/Figure 3/Basal_KEGG_heatmap_SE.pdf",
       width = 8,height = 16,dpi = 1080,units = "in")


```
# Further analysis of Goblet cells ----
Goblet <- subset(all_merged_subset_labelled,
                idents = c("Goblet cells"),
                invert = FALSE)
# Setting HIV Status as idents
Idents(Goblet) <- Goblet$HIV_Status
# Finding Neutrophil markers in
Goblet_Markers <- FindAllMarkers(
  Goblet,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = F,
  min.cells.feature = 3,
  min.cells.group = 3)
# Save Tcells_Markers
write.csv(Goblet_Markers,'scRNAseq_Results/Goblet_Markers.csv',row.names = T)
# Filtering for significant Neut_Markers
Sig_Goblet_Markers <- Goblet_Markers %>%
  dplyr::filter(abs(avg_log2FC)>1.5,
                p_val_adj<0.05)
# Save Sig_Neut_Markers
write.csv(Sig_Goblet_Markers,'scRNAseq_Results/Sig_Goblet_Markers.csv',row.names = T)
# Top10 highly significant genes in neutrophils during HIV
Top10_Goblet_Markers <- Sig_Goblet_Markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=20)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("HIV-","HIV+ ART<3 Months","HIV+ ART>1 Year")
Goblet <- Seurat::SetIdent(Goblet,value = factor(Idents(Goblet),
                                               levels = levels_order))
Top10_Goblet_Markers_plot <- Seurat::DotPlot(Goblet,
                                            features = unique(Top10_Goblet_Markers$gene),
                                            cols = c("blue", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 7,face = 'bold'),
        axis.text.y = element_text(size = 9, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 1))
Top10_Goblet_Markers_plot
# Save Figure 4a
ggsave(Top10_Goblet_Markers_plot,filename="Figures/Figure 9/Fig9a.png",
       width = 10,height = 3.5,dpi = 1080,units = "in")

ggsave(Top10_Goblet_Markers_plot,filename="Figures/Figure 9/Fig9a.pdf",
       width = 10,height = 3.5,dpi = 1080, units = "in")

# Split DEGs by cluster
DEGs_by_cluster <- split(Sig_Goblet_Markers$gene,
                         Sig_Goblet_Markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 10)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Fig9b <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Fig9b
Fig9b <- as.ggplot(Fig9b)

# Save Figure 4f
ggsave(Fig9b,filename="Figures/Figure 9/Fig9b.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Fig9b,filename="Figures/Figure 9/Fig9b.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Fig9b2 <- SeuratExtend::Heatmap(heatmap_matrix,
                                lab_fill='zscore',
                                plot.margin = margin(l = 30),
                                y_text_position = "right",
                                angle = 90,
                                hjust = 0.5,
                                color_scheme = color_palette,
                                vjust = 0.5)
Fig9b2 <- as.ggplot(Fig9b2)
Fig9b2
# Save Figure 4d
ggsave(Fig9b2,filename="Figures/Figure 9/Fig9b2.png",
       width = 6,height = 5,dpi = 1080, units = "in")

ggsave(Fig9b2,filename="Figures/Figure 9/Fig9b2.pdf",
       width = 6,height = 5,dpi = 1080, units = "in")


# Volcano plot of differentially expressed genes in B cells

Sig_Goblet_Markers <- Sig_Goblet_Markers %>%
  dplyr::filter(pct.1>0.05) %>%
  dplyr::mutate(Diff_Expressed = ifelse(avg_log2FC>1.5,"Upregulated","Downregulated"))

Goblet_Volcano_3M <- Sig_Goblet_Markers %>%
  filter(cluster=="HIV+ ART<3 Months") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART<3m")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Goblet_Volcano_3M
# Save Figure 4g
ggsave(Goblet_Volcano_3M,filename="Figures/Figure 9/Fig9c.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Goblet_Volcano_3M,filename="Figures/Figure 9/Fig9c.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Goblet_Volcano_1Y <- Sig_Goblet_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART>1y")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Goblet_Volcano_1Y

# Save Figure 4h
ggsave(Goblet_Volcano_1Y,filename="Figures/Figure 9/Fig9d.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Goblet_Volcano_1Y,filename="Figures/Figure 9/Fig9d.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Goblet_Volcano_HIVNeg <- Sig_Goblet_Markers %>%
  filter(cluster=="HIV-") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "HIV-")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,25))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Goblet_Volcano_HIVNeg  

# Save Figure 4i
ggsave(Goblet_Volcano_HIVNeg,filename="Figures/Figure 9/Fig9e.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Goblet_Volcano_HIVNeg,filename="Figures/Figure 9/Fig9e.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

# Enrichment analysis of differentially expressed genes in B cells
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pathview)

Sig_Goblet_Markers_3M <- Sig_Goblet_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Goblet_Markers_1Y <- Sig_Goblet_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Goblet_Markers_HIVNeg <- Sig_Goblet_Markers %>%
  filter(cluster=="HIV-") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

# Perform GO analysis for each condition
# Enrich GO analysis in Neutrophils from PLHIV on ART<3m
Goblet_go_results_3M <- enrichGO(
  gene = Sig_Goblet_Markers_3M$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Goblet_Go_bar_3M <- barplot(Goblet_go_results_3M, showCategory = 20)
Goblet_Go_bar_3M <- as.ggplot(Goblet_Go_bar_3M)
ggsave(Goblet_Go_bar_3M,filename="Figures/Figure 9/Fig9fGOGobletbar3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Goblet_Go_dot_3M <- clusterProfiler::dotplot(Goblet_go_results_3M, showCategory=20)
Goblet_Go_dot_3M <- as.ggplot(Goblet_Go_dot_3M)
ggsave(Goblet_Go_dot_3M,filename="Figures/Figure 9/Fig9f2GOGobletdot3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from PLHIV on ART>1y
Goblet_go_results_1Y <- enrichGO(
  gene = Sig_Goblet_Markers_1Y$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Goblet_Go_bar_1Y <- barplot(Goblet_go_results_1Y, showCategory = 20)
Goblet_Go_bar_1Y <- as.ggplot(Goblet_Go_bar_1Y)
ggsave(Goblet_Go_bar_1Y,filename="Figures/Figure 9/Fig9gGOGobletbar1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Goblet_Go_dot_1Y <- clusterProfiler::dotplot(Goblet_go_results_1Y,showCategory=20)
Goblet_Go_dot_1Y <- as.ggplot(Goblet_Go_dot_1Y)
ggsave(Goblet_Go_dot_1Y,filename="Figures/Figure 9/Fig9g2GOGobletdot1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from HIV-
Goblet_go_results_HIVneg <- enrichGO(
  gene = Sig_Goblet_Markers_HIVNeg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Goblet_Go_bar_HIVneg <- barplot(Goblet_go_results_HIVneg, showCategory = 20)
Goblet_Go_bar_HIVneg <- as.ggplot(Goblet_Go_bar_HIVneg)
ggsave(Goblet_Go_bar_HIVneg,filename="Figures/Figure 9/Fig9hGOGobletbarNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Goblet_Go_dot_HIVneg <- clusterProfiler::dotplot(Goblet_go_results_HIVneg,showCategory=20)
Goblet_Go_dot_HIVneg <- as.ggplot(Goblet_Go_dot_HIVneg)
ggsave(Goblet_Go_dot_HIVneg,filename="Figures/Figure 9/Fig9h2GOGobletdotNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
```
# Further analysis of Secretory cells ----
Secretory <- subset(all_merged_subset_labelled,
                 idents = c("Secretory cells"),
                 invert = FALSE)
# Setting HIV Status as idents
Idents(Secretory) <- Secretory$HIV_Status
# Finding Neutrophil markers in
Secretory_Markers <- FindAllMarkers(
  Secretory,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = F,
  min.cells.feature = 3,
  min.cells.group = 3)
# Save Tcells_Markers
write.csv(Secretory_Markers,'scRNAseq_Results/Secretory_Markers.csv',row.names = T)
# Filtering for significant Neut_Markers
Sig_Secretory_Markers <- Secretory_Markers %>%
  dplyr::filter(abs(avg_log2FC)>1.5,
                p_val_adj<0.05)
# Save Sig_Neut_Markers
write.csv(Sig_Secretory_Markers,'scRNAseq_Results/Sig_Secretory_Markers.csv',row.names = T)
# Top10 highly significant genes in neutrophils during HIV
Top10_Secretory_Markers <- Sig_Secretory_Markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=20)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("HIV-","HIV+ ART<3 Months","HIV+ ART>1 Year")
Secretory <- Seurat::SetIdent(Secretory,value = factor(Idents(Secretory),
                                                       levels = levels_order))
Top10_Secretory_Markers_plot <- Seurat::DotPlot(Secretory,
                                                features = unique(Top10_Secretory_Markers$gene),
                                                cols = c("blue", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 7,face = 'bold'),
        axis.text.y = element_text(size = 9, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 1))
Top10_Secretory_Markers_plot
# Save Figure 4a
ggsave(Top10_Secretory_Markers_plot,filename="Figures/Figure 10/Fig10a.png",
       width = 10,height = 3.5,dpi = 1080,units = "in")

ggsave(Top10_Secretory_Markers_plot,filename="Figures/Figure 10/Fig10a.pdf",
       width = 10,height = 3.5,dpi = 1080, units = "in")

# Split DEGs by cluster
DEGs_by_cluster <- split(Sig_Secretory_Markers$gene,
                         Sig_Secretory_Markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 10)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Fig10b <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Fig10b
Fig10b <- as.ggplot(Fig10b)

# Save Figure 4f
ggsave(Fig10b,filename="Figures/Figure 10/Fig10b.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Fig10b,filename="Figures/Figure 10/Fig10b.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Fig10b2 <- SeuratExtend::Heatmap(heatmap_matrix,
                                lab_fill='zscore',
                                plot.margin = margin(l = 30),
                                y_text_position = "right",
                                angle = 90,
                                hjust = 0.5,
                                color_scheme = color_palette,
                                vjust = 0.5)
Fig10b2 <- as.ggplot(Fig10b2)
Fig10b2
# Save Figure 4d
ggsave(Fig10b2,filename="Figures/Figure 10/Fig10b2.png",
       width = 6,height = 5,dpi = 1080, units = "in")

ggsave(Fig10b2,filename="Figures/Figure 10/Fig10b2.pdf",
       width = 6,height = 5,dpi = 1080, units = "in")


# Volcano plot of differentially expressed genes in B cells

Sig_Secretory_Markers <- Sig_Secretory_Markers %>%
  dplyr::filter(pct.1>0.05) %>%
  dplyr::mutate(Diff_Expressed = ifelse(avg_log2FC>1.5,"Upregulated","Downregulated"))

Secretory_Volcano_3M <- Sig_Secretory_Markers %>%
  filter(cluster=="HIV+ ART<3 Months") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART<3m")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Secretory_Volcano_3M
# Save Figure 4g
ggsave(Secretory_Volcano_3M,filename="Figures/Figure 10/Fig10c.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Secretory_Volcano_3M,filename="Figures/Figure 10/Fig10c.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Secretory_Volcano_1Y <- Sig_Secretory_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART>1y")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Secretory_Volcano_1Y

# Save Figure 4h
ggsave(Secretory_Volcano_1Y,filename="Figures/Figure 10/Fig10d.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Secretory_Volcano_1Y,filename="Figures/Figure 10/Fig10d.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Secretory_Volcano_HIVNeg <- Sig_Secretory_Markers %>%
  filter(cluster=="HIV-") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "HIV-")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,25))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Secretory_Volcano_HIVNeg  

# Save Figure 4i
ggsave(Secretory_Volcano_HIVNeg,filename="Figures/Figure 10/Fig10e.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Secretory_Volcano_HIVNeg,filename="Figures/Figure 10/Fig10e.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

# Enrichment analysis of differentially expressed genes in B cells
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pathview)

Sig_Secretory_Markers_3M <- Sig_Secretory_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Secretory_Markers_1Y <- Sig_Secretory_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Secretory_Markers_HIVNeg <- Sig_Secretory_Markers %>%
  filter(cluster=="HIV-") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

# Perform GO analysis for each condition
# Enrich GO analysis in Neutrophils from PLHIV on ART<3m
Secretory_go_results_3M <- enrichGO(
  gene = Sig_Secretory_Markers_3M$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Secretory_Go_bar_3M <- barplot(Secretory_go_results_3M, showCategory = 20)
Secretory_Go_bar_3M <- as.ggplot(Secretory_Go_bar_3M)
ggsave(Secretory_Go_bar_3M,filename="Figures/Figure 10/Fig10fGOSecretorybar3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Secretory_Go_dot_3M <- clusterProfiler::dotplot(Secretory_go_results_3M,showCategory=20)
Secretory_Go_dot_3M <- as.ggplot(Secretory_Go_dot_3M)
ggsave(Secretory_Go_dot_3M,filename="Figures/Figure 10/Fig10f2GOSecretorydot3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from PLHIV on ART>1y
Secretory_go_results_1Y <- enrichGO(
  gene = Sig_Secretory_Markers_1Y$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Secretory_Go_bar_1Y <- barplot(Secretory_go_results_1Y, showCategory = 20)
Secretory_Go_bar_1Y <- as.ggplot(Secretory_Go_bar_1Y)
ggsave(Secretory_Go_bar_1Y,filename="Figures/Figure 10/Fig10gGOSecretorybar1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Secretory_Go_dot_1Y <- clusterProfiler::dotplot(Secretory_go_results_1Y,showCategory=20)
Secretory_Go_dot_1Y <- as.ggplot(Secretory_Go_dot_1Y)
ggsave(Secretory_Go_dot_1Y,filename="Figures/Figure 10/Fig10g2GOSecretorydot1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from HIV-
Secretory_go_results_HIVneg <- enrichGO(
  gene = Sig_Secretory_Markers_HIVNeg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Secretory_Go_bar_HIVneg <- barplot(Secretory_go_results_HIVneg, showCategory = 20)
Secretory_Go_bar_HIVneg <- as.ggplot(Secretory_Go_bar_HIVneg)
ggsave(Secretory_Go_bar_HIVneg,filename="Figures/Figure 10/Fig10hGOSecretorybarNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Secretory_Go_dot_HIVneg <- clusterProfiler::dotplot(Secretory_go_results_HIVneg,showCategory=20)
Secretory_Go_dot_HIVneg <- as.ggplot(Secretory_Go_dot_HIVneg)
ggsave(Secretory_Go_dot_HIVneg,filename="Figures/Figure 10/Fig10h2GOSecretorydotNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
```

# Further analysis of Ciliated cells ----
Ciliated <- subset(all_merged_subset_labelled,
                    idents = c("Ciliated cells"),
                    invert = FALSE)
# Setting HIV Status as idents
Idents(Ciliated) <- Ciliated$HIV_Status
# Finding Neutrophil markers in
Ciliated_Markers <- FindAllMarkers(
  Ciliated,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = F,
  min.cells.feature = 3,
  min.cells.group = 3)
# Save Tcells_Markers
write.csv(Ciliated_Markers,'scRNAseq_Results/Ciliated_Markers.csv',row.names = T)
# Filtering for significant Neut_Markers
Sig_Ciliated_Markers <- Ciliated_Markers %>%
  dplyr::filter(abs(avg_log2FC)>1.5,
                p_val_adj<0.05)
# Save Sig_Neut_Markers
write.csv(Sig_Ciliated_Markers,'scRNAseq_Results/Sig_Ciliated_Markers.csv',row.names = T)
# Top10 highly significant genes in neutrophils during HIV
Top10_Ciliated_Markers <- Sig_Ciliated_Markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=20)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("HIV-","HIV+ ART<3 Months","HIV+ ART>1 Year")
Ciliated <- Seurat::SetIdent(Ciliated,value = factor(Idents(Ciliated),
                                                       levels = levels_order))
Top10_Ciliated_Markers_plot <- Seurat::DotPlot(Ciliated,
                                               features = unique(Top10_Ciliated_Markers$gene),
                                               cols = c("blue", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 7,face = 'bold'),
        axis.text.y = element_text(size = 9, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 1))
Top10_Ciliated_Markers_plot
# Save Figure 4a
ggsave(Top10_Ciliated_Markers_plot,filename="Figures/Figure 11/Top10_Ciliated_Markers_plot.png",
       width = 10,height = 3.5,dpi = 1080,units = "in")

ggsave(Top10_Ciliated_Markers_plot,filename="Figures/Figure 11/Fig11a.pdf",
       width = 10,height = 3.5,dpi = 1080, units = "in")

# Split DEGs by cluster
DEGs_by_cluster <- split(Sig_Ciliated_Markers$gene,
                         Sig_Ciliated_Markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 10)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Fig11b <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Fig11b
Fig11b <- as.ggplot(Fig10b)

# Save Figure 4f
ggsave(Fig11b,filename="Figures/Figure 11/Fig11b.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Fig11b,filename="Figures/Figure 11/Fig11b.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Fig11b2 <- SeuratExtend::Heatmap(heatmap_matrix,
                                 lab_fill='zscore',
                                 plot.margin = margin(l = 30),
                                 y_text_position = "right",
                                 angle = 90,
                                 hjust = 0.5,
                                 color_scheme = color_palette,
                                 vjust = 0.5)
Fig11b2 <- as.ggplot(Fig11b2)
Fig11b2
# Save Figure 4d
ggsave(Fig11b2,filename="Figures/Figure 11/Fig11b2.png",
       width = 7,height = 5,dpi = 1080, units = "in")

ggsave(Fig11b2,filename="Figures/Figure 11/Fig11b2.pdf",
       width = 6,height = 5,dpi = 1080, units = "in")


# Volcano plot of differentially expressed genes in B cells

Sig_Ciliated_Markers <- Sig_Ciliated_Markers %>%
  dplyr::filter(pct.1>0.05) %>%
  dplyr::mutate(Diff_Expressed = ifelse(avg_log2FC>1.5,"Upregulated","Downregulated"))

Ciliated_Volcano_3M <- Sig_Ciliated_Markers %>%
  filter(cluster=="HIV+ ART<3 Months") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART<3m")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Ciliated_Volcano_3M
# Save Figure 4g
ggsave(Ciliated_Volcano_3M,filename="Figures/Figure 11/Fig11c.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Ciliated_Volcano_3M,filename="Figures/Figure 11/Fig11c.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Ciliated_Volcano_1Y <- Sig_Ciliated_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART>1y")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Ciliated_Volcano_1Y

# Save Figure 4h
ggsave(Ciliated_Volcano_1Y,filename="Figures/Figure 11/Fig11d.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Ciliated_Volcano_1Y,filename="Figures/Figure 11/Fig11d.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Ciliated_Volcano_HIVNeg <- Sig_Ciliated_Markers %>%
  filter(cluster=="HIV-") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "HIV-")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,25))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Ciliated_Volcano_HIVNeg  

# Save Figure 4i
ggsave(Ciliated_Volcano_HIVNeg,filename="Figures/Figure 11/Fig11e.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Ciliated_Volcano_HIVNeg,filename="Figures/Figure 11/Fig11e.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

# Enrichment analysis of differentially expressed genes in B cells
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pathview)

Sig_Ciliated_Markers_3M <- Sig_Ciliated_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Ciliated_Markers_1Y <- Sig_Ciliated_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Ciliated_Markers_HIVNeg <- Sig_Ciliated_Markers %>%
  filter(cluster=="HIV-") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

# Perform GO analysis for each condition
# Enrich GO analysis in Neutrophils from PLHIV on ART<3m
Ciliated_go_results_3M <- enrichGO(
  gene = Sig_Ciliated_Markers_3M$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Ciliated_Go_bar_3M <- barplot(Ciliated_go_results_3M, showCategory = 20)
Ciliated_Go_bar_3M <- as.ggplot(Ciliated_Go_bar_3M)
ggsave(Ciliated_Go_bar_3M,filename="Figures/Figure 11/Fig11fGOCiliatedbar3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Ciliated_Go_dot_3M <- clusterProfiler::dotplot(Ciliated_go_results_3M,showCategory=20)
Ciliated_Go_dot_3M <- as.ggplot(Ciliated_Go_dot_3M)
ggsave(Ciliated_Go_dot_3M,filename="Figures/Figure 11/Fig11f2GOCiliateddot3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from PLHIV on ART>1y
Ciliated_go_results_1Y <- enrichGO(
  gene = Sig_Ciliated_Markers_1Y$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Ciliated_Go_bar_1Y <- barplot(Ciliated_go_results_1Y, showCategory = 20)
Ciliated_Go_bar_1Y <- as.ggplot(Ciliated_Go_bar_1Y)
ggsave(Ciliated_Go_bar_1Y,filename="Figures/Figure 11/Fig11gGOCiliatedbar1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Ciliated_Go_dot_1Y <- clusterProfiler::dotplot(Ciliated_go_results_3M,showCategory=20)
Ciliated_Go_dot_1Y <- as.ggplot(Ciliated_Go_dot_1Y)
ggsave(Ciliated_Go_dot_3M,filename="Figures/Figure 11/Fig11g2GOCiliateddot1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from HIV-
Ciliated_go_results_HIVneg <- enrichGO(
  gene = Sig_Ciliated_Markers_HIVNeg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Ciliated_Go_bar_HIVneg <- barplot(Ciliated_go_results_HIVneg, showCategory = 20)
Ciliated_Go_bar_HIVneg <- as.ggplot(Ciliated_Go_bar_HIVneg)
ggsave(Ciliated_Go_bar_HIVneg,filename="Figures/Figure 11/Fig11hGOCiliatedbarNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Ciliated_Go_dot_HIVneg <- clusterProfiler::dotplot(Ciliated_go_results_HIVneg,showCategory=20)
Ciliated_Go_dot_HIVneg <- as.ggplot(Ciliated_Go_dot_HIVneg)
ggsave(Ciliated_Go_dot_HIVneg,filename="Figures/Figure 11/Fig11h2GOCiliateddotNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
```
# Further analysis of Squamous cells ----
Squamous <- subset(all_merged_subset_labelled,
                   idents = c("Squamous cells"),
                   invert = FALSE)
# Setting HIV Status as idents
Idents(Squamous) <- Squamous$HIV_Status
# Finding Neutrophil markers in
Squamous_Markers <- FindAllMarkers(
  Squamous,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = F,
  min.cells.feature = 3,
  min.cells.group = 3)
# Save Tcells_Markers
write.csv(Squamous_Markers,'scRNAseq_Results/Squamous_Markers.csv',row.names = T)
# Filtering for significant Neut_Markers
Sig_Squamous_Markers <- Squamous_Markers %>%
  dplyr::filter(abs(avg_log2FC)>1.5,
                p_val_adj<0.05)
# Save Sig_Neut_Markers
write.csv(Sig_Squamous_Markers,'scRNAseq_Results/Sig_Squamous_Markers.csv',row.names = T)
# Top10 highly significant genes in neutrophils during HIV
Top10_Squamous_Markers <- Sig_Squamous_Markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=20)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("HIV-","HIV+ ART<3 Months","HIV+ ART>1 Year")
Squamous <- Seurat::SetIdent(Squamous,value = factor(Idents(Squamous),
                                                     levels = levels_order))
Top10_Squamous_Markers_plot <- Seurat::DotPlot(Squamous,
                                               features = unique(Top10_Squamous_Markers$gene),
                                               cols = c("blue", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 7,face = 'bold'),
        axis.text.y = element_text(size = 9, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 1))
Top10_Squamous_Markers_plot
# Save Figure 4a
ggsave(Top10_Squamous_Markers_plot,filename="Figures/Figure 12/Fig12a.png",
       width = 10,height = 3.5,dpi = 1080,units = "in")

ggsave(Top10_Squamous_Markers_plot,filename="Figures/Figure 12/Fig12a.pdf",
       width = 10,height = 3.5,dpi = 1080, units = "in")

# Split DEGs by cluster
DEGs_by_cluster <- split(Sig_Squamous_Markers$gene,
                         Sig_Squamous_Markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 10)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Fig12b <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Fig12b
Fig12b <- as.ggplot(Fig12b)

# Save Figure 4f
ggsave(Fig12b,filename="Figures/Figure 12/Fig12b.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Fig12b,filename="Figures/Figure 12/Fig12b.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Fig12b2 <- SeuratExtend::Heatmap(heatmap_matrix,
                                 lab_fill='zscore',
                                 plot.margin = margin(l = 30),
                                 y_text_position = "right",
                                 angle = 90,
                                 hjust = 0.5,
                                 color_scheme = color_palette,
                                 vjust = 0.5)
Fig12b2 <- as.ggplot(Fig12b2)
Fig12b2
# Save Figure 4d
ggsave(Fig12b2,filename="Figures/Figure 12/Fig12b2.png",
       width = 7,height = 5,dpi = 1080, units = "in")

ggsave(Fig12b2,filename="Figures/Figure 12/Fig12b2.pdf",
       width = 6,height = 5,dpi = 1080, units = "in")


# Volcano plot of differentially expressed genes in B cells

Sig_Squamous_Markers <- Sig_Squamous_Markers %>%
  dplyr::filter(pct.1>0.05) %>%
  dplyr::mutate(Diff_Expressed = ifelse(avg_log2FC>1.5,"Upregulated","Downregulated"))

Squamous_Volcano_3M <- Sig_Squamous_Markers %>%
  filter(cluster=="HIV+ ART<3 Months") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART<3m")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Squamous_Volcano_3M
# Save Figure 4g
ggsave(Squamous_Volcano_3M,filename="Figures/Figure 12/Fig12c.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Squamous_Volcano_3M,filename="Figures/Figure 12/Fig12c.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Squamous_Volcano_1Y <- Sig_Squamous_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART>1y")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Squamous_Volcano_1Y

# Save Figure 4h
ggsave(Squamous_Volcano_1Y,filename="Figures/Figure 12/Fig12d.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Squamous_Volcano_1Y,filename="Figures/Figure 12/Fig12d.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Squamous_Volcano_HIVNeg <- Sig_Squamous_Markers %>%
  filter(cluster=="HIV-") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "HIV-")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,25))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Squamous_Volcano_HIVNeg  

# Save Figure 4i
ggsave(Squamous_Volcano_HIVNeg,filename="Figures/Figure 12/Fig12e.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Squamous_Volcano_HIVNeg,filename="Figures/Figure 12/Fig12e.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

# Enrichment analysis of differentially expressed genes in B cells
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pathview)

Sig_Squamous_Markers_3M <- Sig_Squamous_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Squamous_Markers_1Y <- Sig_Squamous_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Squamous_Markers_HIVNeg <- Sig_Squamous_Markers %>%
  filter(cluster=="HIV-") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

# Perform GO analysis for each condition
# Enrich GO analysis in Neutrophils from PLHIV on ART<3m
Squamous_go_results_3M <- enrichGO(
  gene = Sig_Squamous_Markers_3M$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Squamous_Go_bar_3M <- barplot(Squamous_go_results_3M, showCategory = 20)
Squamous_Go_bar_3M <- as.ggplot(Squamous_Go_bar_3M)
ggsave(Squamous_Go_bar_3M,filename="Figures/Figure 12/Fig12fGOSquamousbar3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Squamous_Go_dot_3M <- clusterProfiler::dotplot(Squamous_go_results_3M,showCategory=20)
Squamous_Go_dot_3M <- as.ggplot(Squamous_Go_dot_3M)
ggsave(Squamous_Go_dot_3M,filename="Figures/Figure 12/Fig12f2GOSquamousdot3M.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from PLHIV on ART>1y
Squamous_go_results_1Y <- enrichGO(
  gene = Sig_Squamous_Markers_1Y$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Squamous_Go_bar_1Y <- barplot(Squamous_go_results_1Y, showCategory = 20)
Squamous_Go_bar_1Y <- as.ggplot(Squamous_Go_bar_1Y)
ggsave(Squamous_Go_bar_1Y,filename="Figures/Figure 12/Fig12gGOSquamousbar1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Squamous_Go_dot_1Y <- clusterProfiler::dotplot(Squamous_go_results_1Y,showCategory=20)
Squamous_Go_dot_1Y <- as.ggplot(Squamous_Go_dot_1Y)
ggsave(Squamous_Go_dot_1Y,filename="Figures/Figure 12/Fig12g2GOSquamousdot1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from HIV-
Squamous_go_results_HIVneg <- enrichGO(
  gene = Sig_Squamous_Markers_HIVNeg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Squamous_Go_bar_HIVneg <- barplot(Squamous_go_results_HIVneg, showCategory = 20)
Squamous_Go_bar_HIVneg <- as.ggplot(Squamous_Go_bar_HIVneg)
ggsave(Squamous_Go_bar_HIVneg,filename="Figures/Figure 12/Fig12gGOSquamousbarNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Squamous_Go_dot_HIVneg <- clusterProfiler::dotplot(Squamous_go_results_HIVneg,showCategory=20)
Squamous_Go_dot_HIVneg<- as.ggplot(Squamous_Go_dot_HIVneg)
ggsave(Squamous_Go_dot_HIVneg,filename="Figures/Figure 12/Fig12g2GOSquamousdotNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
```
# Further analysis of Ionocytes cells ----
Ionocytes <- subset(all_merged_subset_labelled,
                   idents = c("Ionocytes"),
                   invert = FALSE)
# Setting HIV Status as idents
Idents(Ionocytes) <- Ionocytes$HIV_Status
# Finding Neutrophil markers in
Ionocytes_Markers <- FindAllMarkers(
  Ionocytes,
  assay = "RNA",
  logfc.threshold = 0.1,
  test.use = "wilcox",
  slot = "data",
  min.pct = 0.01,
  verbose = TRUE,
  only.pos = F,
  min.cells.feature = 3,
  min.cells.group = 3)
# Save Tcells_Markers
write.csv(Ionocytes_Markers,'scRNAseq_Results/Ionocytes_Markers.csv',row.names = T)
# Filtering for significant Neut_Markers
Sig_Ionocytes_Markers <- Ionocytes_Markers %>%
  dplyr::filter(abs(avg_log2FC)>1.5,
                p_val_adj<0.05)
# Save Sig_Neut_Markers
write.csv(Sig_Ionocytes_Markers,'scRNAseq_Results/Sig_Ionocytes_Markers.csv',row.names = T)
# Top10 highly significant genes in neutrophils during HIV
Top10_Ionocytes_Markers <- Sig_Ionocytes_Markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=20)
# Dotplot of highly significant genes for each condition in neutrophils
levels_order <- c("HIV-","HIV+ ART<3 Months","HIV+ ART>1 Year")
Ionocytes <- Seurat::SetIdent(Ionocytes,value = factor(Idents(Ionocytes),
                                                     levels = levels_order))
Top10_Ionocytes_Markers_plot <- Seurat::DotPlot(Ionocytes,
                                               features = unique(Top10_Ionocytes_Markers$gene),
                                               cols = c("blue", "red"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size = 7,face = 'bold'),
        axis.text.y = element_text(size = 9, face = 'bold'),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 9))+
  guides(
    size = guide_legend(
      title = "Percent Expressed",
      title.theme = element_text(size = 10),
      keyheight = unit(0.1, "cm")),
    color = guide_colorbar(
      title = "Average Expression",
      title.theme = element_text(size = 10),
      barwidth = 1))
Top10_Ionocytes_Markers_plot
# Save Figure 4a
ggsave(Top10_Ionocytes_Markers_plot,filename="Figures/Figure 13/Fig13a.png",
       width = 10,height = 3.5,dpi = 1080,units = "in")

ggsave(Top10_Ionocytes_Markers_plot,filename="Figures/Figure 13/Fig13a.pdf",
       width = 10,height = 3.5,dpi = 1080, units = "in")

# Split DEGs by cluster
DEGs_by_cluster <- split(Sig_Ionocytes_Markers$gene,
                         Sig_Ionocytes_Markers$cluster)

# Perform actual GO analysis
run_go_analysis <- function(genes, cluster_name){
  enrichGO(
    gene = genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.1
  )
}
# Perform GO analysis for all clusters
go_results <- lapply(names(DEGs_by_cluster),function(cluster){
  cat("Running GO analysis for cluster:", cluster, "\n")
  run_go_analysis(DEGs_by_cluster[[cluster]], cluster)
})

names(go_results) <- names(DEGs_by_cluster)

# Extract the @result slot and process data
zscore_data <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    data <- as.data.frame(cluster@result)
    data <- data %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::select(Description, zScore)
    rownames(data) <- data$Description
    data <- data %>%
      dplyr::select(-Description)
    return(data)
  } else {
    return(NULL)
  }
})

# Start with the first cluster's z-scores
zscore_matrix <- zscore_data[[1]]
colnames(zscore_matrix) <- names(go_results)[1]  # Name the first column by the cluster name

# Loop through remaining clusters and merge
for (i in 2:length(zscore_data)) {
  current_cluster <- zscore_data[[i]]
  colname <- names(go_results)[i]  # Get cluster name for column
  colnames(current_cluster) <- colname
  zscore_matrix <- merge(zscore_matrix, current_cluster, 
                         by = "row.names", all = TRUE)
  rownames(zscore_matrix) <- zscore_matrix$Row.names
  zscore_matrix <- zscore_matrix[, -1]  # Drop redundant Row.names column
}

# Replace NA values with 0
zscore_matrix[is.na(zscore_matrix)] <- 0

# Optional: Rename columns to ensure clarity
colnames(zscore_matrix) <- names(go_results)

# Subset the go_results object
go_results_top20 <- lapply(go_results, function(cluster) {
  if ("enrichResult" %in% class(cluster)) {
    cluster@result %>%
      dplyr::filter(p.adjust<0.1) %>%
      dplyr::arrange(p.adjust,-zScore) %>% # Sort by adjusted p-value
      dplyr::slice_head(n = 10)   # Select the top 20 rows
  } else {
    NULL
  }
})

# Filter out NULL values if any cluster doesn't have results
go_results_top20 <- Filter(Negate(is.null), go_results_top20)

# Check the first cluster's top 20 for example
head(go_results_top20[[1]])

# Extract the @result slot and process data
merged_top20_results <- dplyr::bind_rows(go_results_top20, .id="Cluster") %>%
  dplyr::select(Cluster,Description,zScore)
rownames(merged_top20_results)<-NULL

merged_top20_results <- merged_top20_results %>%
  as.data.frame() %>%
  pivot_wider(names_from = Cluster, values_from = zScore)

#rownames(merged_top20_results)<-merged_top20_results$Description
merged_top20_results[is.na(merged_top20_results)] <- 0


# Convert to matrix for ComplexHeatmap
heatmap_matrix <- as.matrix(merged_top20_results)
heatmap_matrix <- heatmap_matrix[,-1]
all(is.numeric(heatmap_matrix))
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
sum(is.na(heatmap_matrix))
heatmap_matrix[is.na(heatmap_matrix)] <- 0
range(heatmap_matrix, na.rm = TRUE)

rownames(heatmap_matrix) <- merged_top20_results$Description


# Plot the heatmap
library(ComplexHeatmap)
library(circlize)

Fig13b <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Z-Score",
  clustering_distance_rows = "spearman",
  #col = colorRamp2(
  #c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
  #c("lightblue", "white", "red")),
  #col = rev(rainbow(10)),
  col = circlize::colorRamp2(c(0.00000, 24.95776), c("#EBEBEB", "#C3642E")),
  #column_title = "cell clusters",
  #row_title = "GO Terms",
  cluster_rows = F,
  cluster_columns = F)
Fig13b
Fig13b <- as.ggplot(Fig13b)

# Save Figure 4f
ggsave(Fig13b,filename="Figures/Figure 13/Fig13b.png",
       width = 6,height = 12,dpi = 1080,units = "in")

ggsave(Fig13b,filename="Figures/Figure 13/Fig13b.pdf",
       width = 6,height = 12,dpi = 1080, units = "in")

# Heatmap using SeuratExtend
color_palette <- RColorBrewer::brewer.pal(9, "Reds")
Fig13b2 <- SeuratExtend::Heatmap(heatmap_matrix,
                                 lab_fill='zscore',
                                 plot.margin = margin(l = 30),
                                 y_text_position = "right",
                                 angle = 90,
                                 hjust = 0.5,
                                 color_scheme = color_palette,
                                 vjust = 0.5)
Fig13b2 <- as.ggplot(Fig13b2)
Fig13b2
# Save Figure 4d
ggsave(Fig13b2,filename="Figures/Figure 13/Fig13b2.png",
       width = 5.0,height = 5,dpi = 1080, units = "in")

ggsave(Fig13b2,filename="Figures/Figure 13/Fig13b2.pdf",
       width = 5,height = 5,dpi = 1080, units = "in")


# Volcano plot of differentially expressed genes in B cells

Sig_Ionocytes_Markers <- Sig_Ionocytes_Markers %>%
  dplyr::filter(pct.1>0.05) %>%
  dplyr::mutate(Diff_Expressed = ifelse(avg_log2FC>1.5,"Upregulated","Downregulated"))

Ionocytes_Volcano_3M <- Sig_Ionocytes_Markers %>%
  filter(cluster=="HIV+ ART<3 Months") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART<3m")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Ionocytes_Volcano_3M
# Save Figure 4g
ggsave(Ionocytes_Volcano_3M,filename="Figures/Figure 13/Fig13c.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Ionocytes_Volcano_3M,filename="Figures/Figure 13/Fig13c.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Ionocytes_Volcano_1Y <- Sig_Ionocytes_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "PLHIV on ART>1y")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,50))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Ionocytes_Volcano_1Y

# Save Figure 4h
ggsave(Ionocytes_Volcano_1Y,filename="Figures/Figure 13/Fig13d.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Ionocytes_Volcano_1Y,filename="Figures/Figure 13/Fig13d.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

Ionocytes_Volcano_HIVNeg <- Sig_Ionocytes_Markers %>%
  filter(cluster=="HIV-") %>%
  ggplot(aes(avg_log2FC,-log10(p_val_adj), 
             col=factor(Diff_Expressed,levels = c("Upregulated","Downregulated")),
             label=gene))+
  geom_point(size=1)+
  geom_text_repel(max.overlaps = 7,size=2)+
  #geom_vline(xintercept = c(-1.5,1.5), col="red")+
  #geom_hline(yintercept = -log10(0.1), col="blue")+
  labs(x="Avg Log2FC",y="-Log10(Adjusted P value)",title = "HIV-")+
  #scale_x_continuous(limits = c(-7.5,7.5))+
  scale_y_continuous(limits = c(0,25))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.title = element_text(face = 'bold',hjust = 0.5,size = 8))
Ionocytes_Volcano_HIVNeg  

# Save Figure 4i
ggsave(Ionocytes_Volcano_HIVNeg,filename="Figures/Figure 13/Fig13e.png",
       width = 3,height = 5,dpi = 1080,units = "in")

ggsave(Ionocytes_Volcano_HIVNeg,filename="Figures/Figure 13/Fig13e.pdf",
       width = 3,height = 5,dpi = 1080,units = "in")

# Enrichment analysis of differentially expressed genes in B cells
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(pathview)

Sig_Ionocytes_Markers_3M <- Sig_Ionocytes_Markers %>%
  dplyr::filter(cluster=="HIV+ ART<3 Months") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Ionocytes_Markers_1Y <- Sig_Ionocytes_Markers %>%
  filter(cluster=="HIV+ ART>1 Year") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

Sig_Ionocytes_Markers_HIVNeg <- Sig_Ionocytes_Markers %>%
  filter(cluster=="HIV-") %>%
  sort('avg_log2FC') %>%
  dplyr::select(gene)

# Perform GO analysis for each condition
# Enrich GO analysis in Neutrophils from PLHIV on ART<3m
Ionocytes_go_results_3M <- enrichGO(
  gene = Sig_Ionocytes_Markers_3M$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Ionocytes_Go_bar_3M <- barplot(Ionocytes_go_results_3M, showCategory = 20)
Ionocytes_Go_bar_3M <- as.ggplot(Ionocytes_Go_bar_3M)
ggsave(Ionocytes_Go_bar_3M,filename="Figures/Figure 13/Fig13fGOIonocytesbar3M.png",
       width = 6,height = 10,dpi = 1080,units = "in")

Ionocytes_Go_dot_3M <- clusterProfiler::dotplot(Ionocytes_go_results_3M, showCategory = 20)
Ionocytes_Go_dot_3M <- as.ggplot(Ionocytes_Go_dot_3M)
ggsave(Ionocytes_Go_dot_3M,filename="Figures/Figure 13/Fig13f2GOIonocytesdot3M.png",
       width = 6,height = 10,dpi = 1080,units = "in")

# Enrich GO analysis in Neutrophils from PLHIV on ART>1y
Ionocytes_go_results_1Y <- enrichGO(
  gene = Sig_Ionocytes_Markers_1Y$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Ionocytes_Go_bar_1Y <- barplot(Ionocytes_go_results_1Y, showCategory = 20)
Ionocytes_Go_bar_1Y <- as.ggplot(Ionocytes_Go_bar_1Y)
ggsave(Ionocytes_Go_bar_1Y,filename="Figures/Figure 13/Fig13gGOIonocytesbar1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Ionocytes_Go_dot_1Y <- clusterProfiler::dotplot(Ionocytes_go_results_1Y,showCategory=20)
Ionocytes_Go_dot_1Y <- as.ggplot(Ionocytes_Go_dot_1Y)
ggsave(Ionocytes_Go_dot_1Y,filename="Figures/Figure 13/Fig13g2GOIonocytesdot1Y.png",
       width = 6,height = 9,dpi = 1080,units = "in")


# Enrich GO analysis in Neutrophils from HIV-
Ionocytes_go_results_HIVneg <- enrichGO(
  gene = Sig_Ionocytes_Markers_HIVNeg$gene,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1
)

Ionocytes_Go_bar_HIVneg <- barplot(Ionocytes_go_results_HIVneg, showCategory = 20)
Ionocytes_Go_bar_HIVneg <- as.ggplot(Ionocytes_Go_bar_HIVneg)
ggsave(Ionocytes_Go_bar_HIVneg,filename="Figures/Figure 13/Fig13hGOIonocytesbarNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")
Ionocytes_Go_dot_HIVneg <- clusterProfiler::dotplot(Ionocytes_go_results_HIVneg,showCategory=20)
Ionocytes_Go_dot_HIVneg <- as.ggplot(Ionocytes_Go_dot_HIVneg)
ggsave(Ionocytes_Go_dot_HIVneg,filename="Figures/Figure 13/Fig13h2GOIonocytesdotNeg.png",
       width = 6,height = 9,dpi = 1080,units = "in")

# Regulatory Networks in Neutrophils using SCENIC-----
library(SCENIC)
library(GENIE3)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)

# Load seurat files
load("Data/Single_Cell_Data/Immune_cells.RData")
Immune_cells_sce <-as.SingleCellExperiment(Immune_cells)
# Get cluster information from sce object
cellInfo <- colData(Immune_cells_sce)
# Getting the expression matrix from the sce object
exprMat <- counts(Immune_cells_sce)
# Saving into a loom file
loom <- SCopeLoomR::build_loom("Data/Single_Cell_Data/Immune_cells.loom",
                               dgem = exprMat)
# Add cell annotations
loom <- SCENIC::add_cell_annotation(loom,cellInfo)
SCopeLoomR::close_loom(loom)
# Load .loom file
loomPath <- system.file(package = "SCENIC","Data/Single_Cell_Data/Immune_cells.loom")
loom <- SCopeLoomR::open_loom("Data/Single_Cell_Data/Immune_cells.loom")
exprMat <- SCopeLoomR::get_dgem(loom)
#cellInfo <- SCopeLoomR::get_cell_annotation(loom)
SCopeLoomR::close_loom(loom)
# Initialize settings
library(SCENIC)
library(RcisTarget)
data("motifAnnotations_hgnc_v9",package = "RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
ScenicOptions <- SCENIC::initializeScenic(org = 'hgnc',dbDir = "Data/Single_Cell_Data/cisTarget_databases",nCores = 1)
saveRDS(ScenicOptions,file = "Data/Single_Cell_Data/scenicOptions.rds")

# Coexpression network
genesKept <- SCENIC::geneFiltering(exprMat,ScenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
SCENIC::runCorrelation(exprMat_filtered, ScenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, ScenicOptions)
exprMat <- getNormalizedMatrix(ScenicOptions)


# Build and score the GRN
exprMat_log <- log2(exprMat+1)
ScenicOptions@settings$dbs <- ScenicOptions@settings$dbs["10kb"]
#ScenicOptions <- SCENIC::runSCE (ScenicOptions)
ScenicOptions <- SCENIC::runSCENIC_1_coexNetwork2modules(ScenicOptions)
ScenicOptions <- SCENIC::runSCENIC_2createRegulons(ScenicOptions, coexMethod=c("top5perTarget"))
ScenicOptions <- SCENIC::runSCENIC_3_scoreCells(ScenicOptions, exprMat_log)

# Optional: Binarize activity
aucellApp <- plotTsne_AUCellApp(ScenicOptions, exprMat_log)
savedSelection <- shinny::runAPP(aucellApp)
newThresholds <- savedSelection$thresholds
ScenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(ScenicOptions, "aucell_thresholds"))
ScenicOption <- SCENIC::runSCENIC_4_aucell_binarize(ScenicOptions)
tsneAUC(ScenicOptions,aucType="AUC")

# Export:
SaveRDS(cellInfo, file=getDatasetInfo(ScenicOptions, "cellInfo"))
export2loom(ScenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object agains
saveRDS(ScenicOptions,file = "HIV PAPER/Data/Single_Cell_Data/ScenicOptions.Rds")

# exploring output
motifEnrichment_selfMotifs_wGenes <- loadInt(ScenicOptions,"motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedRFs=="Sox8"]
viewMotifs(tableSubset)

# projecting cells onto a dataset using SINGLER----
library(data.table)
library(SingleR)
# Load Moxon data
Moxon_data <- readRDS("Data/Single_Cell_Data/integrated_annotated_nasal.rds")
DimPlot(Moxon_data,reduction = "umap")
# Load Alex Shalek Data
counts <- fread("Data/Single_Cell_Data/Alex_Shalek_Data/expression/20210220_NasalSwab_RawCounts.txt",
                sep = "\t")
names(counts)[1] <- ""
metadata <- fread("Data/Single_Cell_Data/Alex_Shalek_Data/metadata/NasalSwab_MetaData.txt",
                  sep = "\t")
names(metadata)[1] <- ""
metadata <- metadata[-1]
normalised <- fread("Data/Single_Cell_Data/Alex_Shalek_Data/expression/20210220_NasalSwab_NormCounts.txt",
                    sep = "\t")

Alex_Shalek <- Seurat::CreateSeuratObject(counts = counts,
                                          meta.data = metadata)  


anchors <- Seurat::FindTransferAnchors(
  reference = Moxon_data,
  query = Epithelial_cells,
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  dims = 1:30
)

Epithelial_cells <- Seurat::MapQuery(
  anchorset = anchors,
  reference = Moxon_data,
  query = Epithelial_cells,
  refdata = list(celltype = "manual_annotation"),
  reference.reduction = "pca"
)

DimPlot(Epithelial_cells,
        reduction = 'umap.harmony',
        group.by = "predicted.celltype")

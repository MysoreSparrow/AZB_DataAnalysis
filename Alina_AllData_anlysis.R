library(tidyverse)

# have created a categorical variable for number of replicates to be used later!
# have created a categorical variable for strain number to be used later!

# STRAIN 476

t476 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_476_vs_C57BL6_ctl.csv")

# Adding the table number as a column
t476$Strain_name <- "S476"
t476$N_infected_Replicates <- "Double"

# Removing Unwanted Columns
t476 <- select(t476, -lfcSE, -stat, -symbol)

# Renaming Columns
t476 <- t476 %>% rename(CPM_infected_R1 = CPM_.C57BL6_476_R1.) %>% rename(CPM_infected_R2 = CPM_.C57BL6_476_R2.) %>% rename(CPM_ctl_R1 = CPM_.C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM_.C57BL6_ctl_R2.)
head(t476)

dim(t476)

# STRAIN 754
t754 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_754_vs_C57BL6_ctl.csv")

# Adding the table number as a column
t754$Strain_name <- "S754"
t754$CPM_infected_R2 <- ""
t754$N_infected_Replicates <- "Single"
head(t754)

# Removing Unwanted Columns
t754 <- select(t754, -lfcSE, -stat, -symbol)
head(t754)

str(t754)

# Renaming Columns
t754 <- t754 %>% rename(CPM_infected_R1 = CPM..C57BL6_754_R1.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)

# Reordering Columns to be or correct order
t754 <- t754[, c(1, 2, 3, 4, 5, 6, 10, 7, 8, 9,11)]
head(t754)
dim(t754)


# STRAIN 755

t755 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_755_vs_C57BL6_ctl.csv")

# Adding the table number as a column
t755$Strain_name <- "S755"
t755$CPM_infected_R2 <- ""
t755$N_infected_Replicates <- "Single"
head(t755)

# Removing Unwanted Columns
t755 <- select(t755, -lfcSE, -stat, -symbol)
head(t755)

str(t755)

# Renaming Columns
t755 <- t755 %>% rename(CPM_infected_R1 = CPM..C57BL6_755_R1.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)

# Reordering Columns to be or correct order
t755 <- t755[, c(1, 2, 3, 4, 5, 6, 10, 7, 8, 9,11)]
head(t755)
dim(t755)



# STRAIN 757

t757 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_757_vs_C57BL6_ctl.csv")

# Removing Unwanted Columns
t757 <- select(t757, -lfcSE, -stat, -symbol)

# Renaming Columns
t757 <- t757 %>% rename(CPM_infected_R1 = CPM..C57BL6_757_R1.) %>% rename(CPM_infected_R2 = CPM..C57BL6_757_R2.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)

# Adding the table number as a column
t757$Strain_name <- "S757"
t757$N_infected_Replicates <- "Double"
head(t757)
dim(t757)


# Strain 758

t758 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_758_vs_C57BL6_ctl.csv")

# Adding the table number as a column
t758$Strain_name <- "S758"
t758$CPM_infected_R2 <- ""
t758$N_infected_Replicates <- "Single"
head(t758)
str(t758)
# Removing Unwanted Columns
t758 <- select(t758, -lfcSE, -stat, -symbol)

# Renaming Columns
t758 <- t758 %>% rename(CPM_infected_R1 = CPM..C57BL6_758_R1.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)

# Reordering Columns to be or correct order
t758 <- t758[, c(1, 2, 3, 4, 5, 6, 10, 7, 8, 9,11)]
head(t758)
dim(t758)

# Strain 760

t760 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_760_vs_C57BL6_ctl.csv")

# Adding the table number as a column
t760$Strain_name <- "S760"
t760$CPM_infected_R2 <- ""
t760$N_infected_Replicates <- "Single"
head(t760)
str(t760)
# Removing Unwanted Columns
t760 <- select(t760, -lfcSE, -stat, -symbol)

# Renaming Columns
t760 <- t760 %>% rename(CPM_infected_R1 = CPM..C57BL6_760_R1.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)

# Reordering Columns to be or correct order
t760 <- t760[, c(1, 2, 3, 4, 5, 6, 10, 7, 8, 9,11)]
head(t760)
dim(t760)

# Strain 761

t761 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_761_vs_C57BL6_ctl.csv")

# Adding the table number as a column
t761$Strain_name <- "S761"
t761$CPM_infected_R2 <- ""
t761$N_infected_Replicates <- "Single"
head(t761)
str(t761)
# Removing Unwanted Columns
t761 <- select(t761, -lfcSE, -stat, -symbol)

# Renaming Columns
t761 <- t761 %>% rename(CPM_infected_R1 = CPM..C57BL6_761_R1.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)

# Reordering Columns to be or correct order
t761 <- t761[, c(1, 2, 3, 4, 5, 6, 10, 7, 8, 9,11)]
head(t761)
dim(t761)


# Strain 762

t762 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_762_vs_C57BL6_ctl.csv")

# Adding the table number as a column
t762$Strain_name <- "S762"
t762$CPM_infected_R2 <- ""
t762$N_infected_Replicates <- "Single"
str(t762)
# Removing Unwanted Columns
t762 <- select(t762, -lfcSE, -stat, -symbol)
# Renaming Columns
t762 <- t762 %>% rename(CPM_infected_R1 = CPM..C57BL6_762_R1.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)
# Reordering Columns to be or correct order
t762 <- t762[, c(1, 2, 3, 4, 5, 6, 10, 7, 8, 9,11)]
head(t762)


# Strain 763
t763 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_763_vs_C57BL6_ctl.csv")

# Adding the table number as a column
t763$Strain_name <- "S763"
t763$CPM_infected_R2 <- ""
t763$N_infected_Replicates <- "Single"
str(t763)
# Removing Unwanted Columns
t763 <- select(t763, -lfcSE, -stat, -symbol)
# Renaming Columns
t763 <- t763 %>% rename(CPM_infected_R1 = CPM..C57BL6_763_R1.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)
# Reordering Columns to be or correct order
t763 <- t763[, c(1, 2, 3, 4, 5, 6, 10, 7, 8, 9,11)]
head(t763)


# Strain 764
t764 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_764_vs_C57BL6_ctl.csv")

# Adding the table number as a column
t764$Strain_name <- "S764"
t764$CPM_infected_R2 <- ""
t764$N_infected_Replicates <- "Single"
str(t764)
# Removing Unwanted Columns
t764 <- select(t764, -lfcSE, -stat, -symbol)
# Renaming Columns
t764 <- t764 %>% rename(CPM_infected_R1 = CPM..C57BL6_764_R1.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)
# Reordering Columns to be or correct order
t764 <- t764[, c(1, 2, 3, 4, 5, 6, 10, 7, 8, 9,11)]
head(t764)


# Strain 765
t765 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_765_vs_C57BL6_ctl.csv")

# Adding the table number as a column
t765$Strain_name <- "S765"
t765$CPM_infected_R2 <- ""
t765$N_infected_Replicates <- "Single"
str(t765)
# Removing Unwanted Columns
t765 <- select(t765, -lfcSE, -stat, -symbol)
# Renaming Columns
t765 <- t765 %>% rename(CPM_infected_R1 = CPM..C57BL6_765_R1.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)
# Reordering Columns to be or correct order
t765 <- t765[, c(1, 2, 3, 4, 5, 6, 10, 7, 8, 9,11)]
head(t765)

# Strain 766
t766 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_766_vs_C57BL6_ctl.csv")

# Removing Unwanted Columns
t766 <- select(t766, -lfcSE, -stat, -symbol)

# Renaming Columns
t766 <- t766 %>% rename(CPM_infected_R1 = CPM..C57BL6_766_R1.) %>% rename(CPM_infected_R2 = CPM..C57BL6_766_R2.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)

# Adding the table number as a column
t766$Strain_name <- "S766"
t766$N_infected_Replicates <- "Double"
head(t766)
dim(t766)


# Strain 768
t768 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_768_vs_C57BL6_ctl.csv")

# Removing Unwanted Columns
t768 <- select(t768, -lfcSE, -stat, -symbol)

# Renaming Columns
t768 <- t768 %>% rename(CPM_infected_R1 = CPM..C57BL6_768_R1.) %>% rename(CPM_infected_R2 = CPM..C57BL6_768_R2.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)

# Adding the table number as a column
t768$Strain_name <- "S768"
t768$N_infected_Replicates <- "Double"
head(t768)
dim(t768)

# Strain 769
t769 <- read.csv("~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/C57BL6_769_vs_C57BL6_ctl.csv")

# Adding the table number as a column
t769$Strain_name <- "S769"
t769$CPM_infected_R2 <- ""
t769$N_infected_Replicates <- "Single"
str(t769)
# Removing Unwanted Columns
t769 <- select(t769, -lfcSE, -stat, -symbol)
# Renaming Columns
t769 <- t769 %>% rename(CPM_infected_R1 = CPM..C57BL6_769_R1.) %>% rename(CPM_ctl_R1 = CPM..C57BL6_ctl_R1.)  %>% rename(CPM_ctl_R2 = CPM..C57BL6_ctl_R2.)
# Reordering Columns to be or correct order
t769 <- t769[, c(1, 2, 3, 4, 5, 6, 10, 7, 8, 9,11)]
head(t769)

#rbind two data frames into one data frame
binded_frame <- rbind(t476,t754,t755,t757,t758,t760,t761,t762,t763,t764,t765,t766, t768, t769)
write.csv(binded_frame,"~/R/Rtuts/Data/WG__Discussion_Alina's_EPEC_project/csv/Master_Data_Alina.csv")

dim(binded_frame)

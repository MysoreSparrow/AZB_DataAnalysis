library(tidyverse)
C57BL6_766_vs_C57BL6_ctl <- read_csv("D:/C57BL6_766_vs_C57BL6_ctl.csv")

C57BL6_766_vs_C57BL6_ctl <- as_tibble(C57BL6_766_vs_C57BL6_ctl)

head(C57BL6_766_vs_C57BL6_ctl)

str(C57BL6_766_vs_C57BL6_ctl)

#plot1 <- ggplot(data = C57BL6_766_vs_C57BL6_ctl, mapping = aes(x = log2FoldChange, y = pvalue  ))
#plot1 + geom_point(alpha=0.1)

#plot2 <- ggplot(data = C57BL6_766_vs_C57BL6_ctl, mapping = aes(x = log2FoldChange, y = baseMean  ))
#plot2 + geom_point()

#plot3 <- ggplot(data = C57BL6_766_vs_C57BL6_ctl, mapping = aes(x = CPM_C57BL6_766_R1, y= CPM_C57BL6_ctl_R1 ))
#plot3 + geom_point()+ geom_smooth()

#plot4 <- ggplot(data = C57BL6_766_vs_C57BL6_ctl, mapping = aes(x = log2FoldChange, y= CPM_C57BL6_766_R2 ))
#plot4 + geom_point() + geom_smooth()


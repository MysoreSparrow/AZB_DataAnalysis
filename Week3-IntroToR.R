library(tidyverse)

metabric <- read.csv("D:/METABRIC/METABRIC_RNA_Mutation.csv")
#str(metabric)

metabric_mini <- metabric[1:1000,20:50]
#str(metabric_mini)
#view(metabric_mini)


#Basic Scatter Plot


#ggplot(data = metabric_mini) + geom_point(mapping = aes(x = lymph_nodes_examined_positive, y= tumor_size, colour = pr_status))

##contouring points based on a continous variable
#ggplot(data = metabric_mini) + geom_point(mapping = aes(x = overall_survival_months, y= tumor_size, colour = tumor_stage))
ggplot(data = metabric_mini) + geom_point(mapping = aes(x = overall_survival_months, y= tumor_size, colour = tumor_stage, shape = pr_status), size=0.9, alpha = 0.65)

#?geom_bar

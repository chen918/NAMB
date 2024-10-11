library(reshape)
library(ggplot2)

dat <- openxlsx::read.xlsx("robust_dat.xlsx")

plotdata <- melt(dat,
                 id.vars=c("rho","group"))

plotdata$variable <- gsub("_"," ",plotdata$variable)
plotdata$group <- gsub("_"," ",plotdata$group)

#plotdata$rho <- as.character(plotdata$rho)

colnames(plotdata)[3] <- "Model for estimation"
colnames(plotdata)[1] <- "Rho"

plotdata$cut <- 0

plotdata$group <- factor(plotdata$group,levels=c("Probit Latent","Logistic Latent","Probit Mean","Logistic Mean"))
plotdata$`Model for estimation` <- factor(plotdata$`Model for estimation`,levels=c("Probit Latent","Logistic Latent","Probit Mean","Logistic Mean"))

plotdata$l <- ifelse(plotdata$group==plotdata$`Model for estimation`,plotdata$value,NA)
plotdata$m <- ifelse(plotdata$group!=plotdata$`Model for estimation`,plotdata$value,NA)



levels(plotdata$group)=c("Figure 1.a. Generation Probit Latent","Figure 1.b. Generation Logistic Latent","Figure 1.c. Generation Probit Mean","Figure 1.d. Generation Logistic Mean")


ggplot(plotdata) +
  geom_point(aes(x=Rho,y=value,col=`Model for estimation`,shape=`Model for estimation`),size=3) +
  geom_line(aes(x=Rho,y=m,col=`Model for estimation`),linetype="dashed", alpha=0.3) +
  geom_line(aes(x=Rho,y=l,col=`Model for estimation`)) +
  geom_hline(aes(yintercept = cut),show.legend = FALSE, linetype="dashed") +
  facet_wrap(~group,nrow=2,scales = "free") +
  ylab("Bias") +
  theme(text=element_text(size=15)) +
  scale_x_continuous(limits = c(-0.2,0.2), breaks =c(-0.2,0,0.2))

ggsave("/Users/CongMeow/Desktop/Robust_plot.png",dpi = 600,width = 12,height = 8) 

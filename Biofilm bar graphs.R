##Bar graphs for biofilm PC
my.data<-read_excel("D5F3.xlsx", sheet=2, col_names = TRUE)
# Calculates mean, sd, se and IC
my_sum <- my.data %>%
  group_by(Treatment) %>%
  summarise(
    n=n(),
    mean=mean(Concentration),
    sd=sd(Concentration),
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
# insert mean and standard error APC 
#add PC
p115<-c(8.427,8.789,8.882,8.903)
s115<-c(0.039,0.078,0.134,0.103)
p315<-c(8.452,8.630,8.980,8.721)
s315<-c(0.123,0.093,0.159,0.174)
p215<-c(8.700,8.658,8.961,9.019)
s215<-c(0.095,0.057,0.036,0.001)
p13<-c(6.107,7.163,8.292,8.308)
s13<-c(0.336,0.365,0.184,0.184)
p23<-c(5.601,6.725,7.816,7.946)
s23<-c(0.092,0.280,0.077,0.119)
p33<-c(4.989,7.026,7.944,8.047)
s33<-c(0.273,0.348,0.676,0.588)
p15<-c(5.837,7.978,7.471,7.256)
s15<-c(0.311,0.093,0.247,0.207)
p25<-c(5.948,7.050,7.153,6.980)
s25<-c(0.335,0.212,0.192,0.270)
p35<-c(6.501,7.817,7.651,7.559)
s35<-c(0.118,0.127,0.278,0.313)
my_sum<-add_column(my_sum, d=p315)
my_sum<-add_column(my_sum, e=s315)



# Standard Error
x<-ggplot(data=my_sum, aes(x=Treatment, y=10)) +
  geom_bar( aes(x=Treatment, y=mean), stat="identity", fill="blue") +
  geom_errorbar( aes(x=Treatment, ymin=mean-se, ymax=mean+se), width=0.4, colour="orange", alpha=0.9, size=1.5) +
  geom_line(aes(x=Treatment, y=d), color="black", group=1) +
  geom_errorbar( aes(x=Treatment, ymin=d-e, ymax=d+e), width=0.3, colour="brown", alpha=0.9, size=1) +
  geom_text(aes(label=mean), vjust=25, color="black", size=3.5) +
  geom_text(aes(label=d), vjust=7, color="black", size=3.5) +
  geom_hline(yintercept = 1.52, color='grey30', linetype="dashed", size=1) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))
x

#spot innoc analysis

#ed 107
my.data<-read_excel("TempTrialresultswithFDAstrains.xlsx", sheet=4, col_names = TRUE)
my.data$Temp <- factor(my.data$Temp, levels = rev(levels(factor(my.data$Temp))))
one.way<-aov(Concentration~Temp, data=my.data)
summary(one.way)
tukey.two.way<-TukeyHSD(one.way)

#ed 108

my.data<-read_excel("TempTrialresultswithFDAstrains.xlsx", sheet=3, col_names = TRUE)
my.data$Temp <- factor(my.data$Temp, levels = rev(levels(factor(my.data$Temp))))
one.way<-aov(Concentration~Temp, data=my.data)
summary(one.way)
tukey.two.way<-TukeyHSD(one.way)
tukey.two.way

#both at 107

my.data<-read_excel("TempTrialresultswithFDAstrains.xlsx", sheet=6, col_names = TRUE)

ed.test<-t.test(Conc~ED, data=my.data, paired = TRUE, alternative = "two.sided")
summary(ed.test)

#both at 108

my.data<-read_excel("TempTrialresultswithFDAstrains.xlsx", sheet=7, col_names = TRUE)
ll.test<-t.test(Conc~LL, data=my.data, paired = TRUE, alternative = "two.sided")
summary(ll.test)

# spot innoc graph
  fdata<- read.csv("LLspot.csv", check.names = FALSE)

  LL8my_sum <- fdata %>%
  group_by(Temp) %>%
  summarise(
    n=n(),
    mean=mean(LL),
    sd=sd(LL),
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

  a=as.numeric(fdata$mean)-as.numeric(fdata$se)
  b=as.numeric(fdata$mean)+as.numeric(fdata$se)



#PS01155 barplot

  my.data<-read_excel("TempTrialresultswithFDAstrains.xlsx", sheet=7, col_names = TRUE)
  APC_stat<-describeBy(my.data$LL, list(my.data$Conc, my.data$Temp), mat = TRUE)

  LL<-ggplot(APC_stat, aes(x=group2, y=mean, fill=group1))+
    geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(-.9)) +
    scale_color_brewer(palette="Dark2") +
    theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
    theme(axis.title = element_text(size=15,color='black')) +
    theme(panel.background = element_rect(fill='transparent', color = NA),
          plot.background = element_rect(fill = 'transparent',color = NA),
          panel.border = element_rect(color='black',fill = NA,size=1))

#PS01156 barplot

  my.data<-read_excel("TempTrialresultswithFDAstrains.xlsx", sheet=6, col_names = TRUE)
  edd<-describeBy(my.data$ED, list(my.data$Conc, my.data$Temp), mat = TRUE)

  ED<-ggplot(edd, aes(x=group2, y=mean, fill=group1))+
    geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(-.9)) +
    theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
    theme(axis.title = element_text(size=15,color='black')) +
    theme(panel.background = element_rect(fill='transparent', color = NA),
          plot.background = element_rect(fill = 'transparent',color = NA),
          panel.border = element_rect(color='black',fill = NA,size=1))



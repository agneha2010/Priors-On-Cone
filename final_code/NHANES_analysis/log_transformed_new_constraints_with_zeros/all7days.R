result1 = read.csv("ydata1.log.csv")
result2 = read.csv("ydata2.log.csv")
result3 = read.csv("ydata3.log.csv")
result4 = read.csv("ydata4.log.csv")
result5 = read.csv("ydata5.log.csv")
result6 = read.csv("ydata6.log.csv")
result7 = read.csv("ydata7.log.csv")

result = bind_rows(result1[,-1],result2[,-1],result3[,-1],result4[,-1],result5[,-1],
                   result6[,-1],result7[,-1],.id="id")
result.exp = exp(result[,3:6])
result$id = as.factor(result$id)
#levels(result$id) = c("Day 1","Day 2","Day 3","Day 4","Day 5","Day 6","Day 7")
levels(result$id) = c("Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday")


df.plot = melt(result,id.vars = c("x","id"))
df.plot$Method = as.factor(df.plot$variable)
levels(df.plot$Method) = c("Restricted","Unrestricted spline",
                             "Unrestricted","Restricted spline")


day1 = data.frame(x=1:Ti,y=exp(colMeans(ydata1.log)))
day2 = data.frame(x=1:Ti,y=exp(colMeans(ydata2.log)))
day3 = data.frame(x=1:Ti,y=exp(colMeans(ydata3.log)))
day4 = data.frame(x=1:Ti,y=exp(colMeans(ydata4.log)))
day5 = data.frame(x=1:Ti,y=exp(colMeans(ydata5.log)))
day6 = data.frame(x=1:Ti,y=exp(colMeans(ydata6.log)))
day7 = data.frame(x=1:Ti,y=exp(colMeans(ydata7.log)))

df.true = bind_rows(day1,day2,day3,day4,day5,day6,day7,.id="id")
df.true$id = as.factor(df.true$id)
#levels(df.true$id) = c("Day 1","Day 2","Day 3","Day 4","Day 5","Day 6","Day 7")
levels(df.true$id) = c("Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday")

p1 = ggplot() +
     geom_line(data = df.plot,aes(x=x,y=value,color=Method)) +
     #  scale_linetype_manual(values=c("dashed","solid","solid","solid"))+
     geom_point(data = df.true,aes(x=x,y=y)) + 
     ylab("Activity count") + 
     theme_bw() +
     xlab("Time (in hours)") + facet_wrap(~id, dir="h") +
     theme(text = element_text(size = 17),
           legend.position=c(.7,.1),
           plot.title = element_text(size = 15, face = "bold"),
           legend.title=element_text(size=15),
           legend.text=element_text(size=15))
p1  

ggsave('all7days.png', p1, width = 17,height=11, dpi = 300)

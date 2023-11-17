result1 = read.csv("result_day1.csv")
result2 = read.csv("result_day2.csv")
result3 = read.csv("result_day3.csv")
result4 = read.csv("result_day4.csv")
result5 = read.csv("result_day5.csv")
result6 = read.csv("result_day6.csv")
result7 = read.csv("result_day7.csv")

result = bind_rows(result1[,-1],result2[,-1],result3[,-1],result4[,-1],result5[,-1],
                   result6[,-1],result7[,-1],.id="id")
result$id = as.factor(result$id)
levels(result$id) = c("Day 1","Day 2","Day 3","Day 4","Day 5","Day 6","Day 7")

df.plot = melt(result,id.vars = c("x","id"))
df.plot$Method = as.factor(df.plot$variable)
levels(df.plot$variable) = c("Restricted","Unrestricted spline",
                             "Unrestricted","Restricted spline")


day1 = data.frame(x=1:Ti,y=colMeans(ydata1))
day2 = data.frame(x=1:Ti,y=colMeans(ydata2))
day3 = data.frame(x=1:Ti,y=colMeans(ydata3))
day4 = data.frame(x=1:Ti,y=colMeans(ydata4))
day5 = data.frame(x=1:Ti,y=colMeans(ydata5))
day6 = data.frame(x=1:Ti,y=colMeans(ydata6))
day7 = data.frame(x=1:Ti,y=colMeans(ydata7))

df.true = bind_rows(day1,day2,day3,day4,day5,day6,day7,.id="id")
df.true$id = as.factor(df.true$id)
levels(df.true$id) = c("Day 1","Day 2","Day 3","Day 4","Day 5","Day 6","Day 7")

p1 = ggplot() +
     geom_line(data = df.plot,aes(x=x,y=value,color=Method)) +
     #  scale_linetype_manual(values=c("dashed","solid","solid","solid"))+
     geom_point(data = df.true,aes(x=x,y=y)) + 
     ylab("Activity count") + 
     theme_bw() +
     xlab("Time (in hours)") + facet_wrap(~id) +
     theme(text = element_text(size = 17),
           legend.position=c(.8,.1),
           plot.title = element_text(size = 15, face = "bold"),
           legend.title=element_text(size=15),
           legend.text=element_text(size=15))
p1  

ggsave('all7days.png', p1, width = 17,height=11, dpi = 300)

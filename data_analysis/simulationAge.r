
dataLine <- filter(groupedData,survivalRate <= 0.90,survivalRate >= 0.10)

# time analysis
x = (0:(1500/25))*25
y= data.frame(w=0)
for (a in x) {
   y[nrow(y)+1,] = nrow(data[data$Time>a,])/nrow(data)
}
plot(x,y[2:62,])
# print(filter(data, Time>5000, Result == "EXTINCT"))
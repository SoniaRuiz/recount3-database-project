library(pwr)


pwr::pwr.2p.test(h = ES.h(p1 = 0.60, p2 = 0.5), 
                 sig.level = 0.05,
                 power = 0.8, 
                 alternative = "greater")


pwr::pwr.2p.test(h = ES.h(p1 = 0.95, p2 = 0.50), 
                 sig.level = 0.05,
                 power = 0.8, 
                 alternative = "greater")



# powerchange <- data.frame(p1, power = power1$power * 100)
# plot(powerchange$p1, 
#      powerchange$power, 
#      type = "b", 
#      xlab = "Proportion of Responders in Treatment A", 
#      ylab = "Power (%)")
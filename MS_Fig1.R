#################################################################################
#### Coexistence theory and the frequency-dependence of priority effects 
#### Ke & Letten (2018) Nature Ecology & Evolution
#### This R script creates Figure 1 in the main text
#### Effect of changing species' consumption vector and the supply ratio of two resources in a consumer-resource model on the fitness ratio and stabilization potential (niche difference) of coexistence theory.
#################################################################################




######################################
#### Load packages
######################################
library(tidyr)
library(ggplot2)
library(cowplot)
library(Cairo)




######################################
#### Parameter and saving space setup
######################################

#### Setup Tilman's consumer resource model parameters 
D = 0.7
c11 = 2; c12 = 4; c21 = 4; c22 = 2
w11 = 2; w12 = 4; w21 = 4; w22 = 2
k1 = 0.4; k2 = 0.4
r1 = 1; r2 = 1
T1 = 0.1; T2 = 0.1


#### Calculate slope and intercept for the ZNGI's
slope.1 = - (w11 / w12)
slope.2 = - (w21 / w22)
y.inter.1 = (D * (k1 - T1) + r1 * T1) / (w12 * (r1 - D))
y.inter.2 = (w21 / w22) * (D * (k2 - T2) + r2 * T2) / (w21 * (r2 - D))


#### Calculate equilibrium resource value
B1 = y.inter.1
B2 = y.inter.2
Lamda.1 = (w11 / w12)
Lamda.2 = (w21 / w22)
R1.star = (B1 - B2) / (Lamda.1 - Lamda.2)
R2.star = (B2 * Lamda.1 - B1 * Lamda.2) / (Lamda.1 - Lamda.2)  


#### Setup simulation angle (here, flip theta.1 and theta.2 to create priority effects)
theta.1 = acos(c12 / sqrt(c11^2 + c12^2)) * (180 / pi)
theta.2 = acos(c22 / sqrt(c21^2 + c22^2)) * (180 / pi)
theta.diff.1 = seq(theta.1, theta.2, length=99*2)
theta.diff.2 = seq(theta.2, theta.1, length=99*2)


#### Create saving space
#### Set 1 (square, a supply point with S1 > S2)
c11.save = c(); c12.save = c(); c21.save = c(); c22.save = c()
K1.save = c(); K2.save = c(); alpha.save = c(); beta.save = c()
impact.overlap.save = c(); status.save = c()
LV.slope.1.save = c(); LV.slope.2.save = c(); LV.inter.1.save = c(); LV.inter.2.save = c()

#### Set 2 (circle, a supply point with S1 < S2)
c11.save.2 = c(); c12.save.2 = c(); c21.save.2 = c(); c22.save.2 = c()
K1.save.2 = c(); K2.save.2 = c(); alpha.save.2 = c(); beta.save.2 = c()
impact.overlap.save.2 = c(); status.save.2 = c()
LV.slope.1.save.2 = c(); LV.slope.2.save.2 = c(); LV.inter.1.save.2 = c(); LV.inter.2.save.2 = c()


#### Setup resource supply rate (this affects fitness difference between species)
#### Set 1 (square, a supply point with S1 > S2)
S1 = 0.248; S2 = 0.245

#### Set 2 (circle, a supply point with S1 < S2)
S1.2 = 0.3; S2.2 = 0.38




######################################
#### Simulation for impact vector
######################################

#### Run calculation on different impact vector setup
for(i in 1:(99*2))
{
  #### Set 1 (square, a supply point with S1 > S2)
  # 1st -- get new impact vector for the new angle setting
  c11.save[i]=sin(theta.diff.1[i] * (pi / 180)) * (sqrt(c11^2 + c12^2))
  c12.save[i]=cos(theta.diff.1[i] * (pi / 180)) * (sqrt(c11^2 + c12^2))
  c21.save[i]=sin(theta.diff.2[i] * (pi / 180)) * (sqrt(c21^2 + c22^2))
  c22.save[i]=cos(theta.diff.2[i] * (pi / 180)) * (sqrt(c21^2 + c22^2))
  
  # 2nd -- get new LV-components for the new angle setting
  K1.save[i] = (D * (S2 + S1 * Lamda.1 - B1)) / (c12.save[i] + c11.save[i] * Lamda.1)
  K2.save[i] = (D * (S2 + S1 * Lamda.2 - B2)) / (c22.save[i] + c21.save[i] * Lamda.2)
  alpha.save[i] = (c22.save[i] + c21.save[i] * Lamda.1) / (c12.save[i] + c11.save[i] * Lamda.1)
  beta.save[i]  = (c12.save[i] + c11.save[i] * Lamda.2) / (c22.save[i] + c21.save[i] * Lamda.2)
  
  # 3rd -- get new LV-ZNGI for the new angle setting
  LV.slope.1.save[i] = -(1 / alpha.save[i])
  LV.slope.2.save[i] = -(beta.save[i])
  LV.inter.1.save[i] = K1.save[i] / alpha.save[i]
  LV.inter.2.save[i] = K2.save[i]
  
  # 4th -- get impact.niche.overlap for the new angle setting
  status.save[i] = ifelse(theta.diff.2[i] > theta.diff.1[i], 0, 1)  
  impact.overlap.save[i] = (c11.save[i] * c21.save[i] + c12.save[i] * c22.save[i]) / (sqrt(c11.save[i]^2 + c12.save[i]^2) * sqrt(c21.save[i]^2 + c22.save[i]^2))
  
  
  #### Set 2 (circle, a supply point with S1 < S2)
  # 1st -- get new impact vector for the new angle setting
  c11.save.2[i]=sin(theta.diff.1[i] * (pi / 180)) * (sqrt(c11^2 + c12^2))
  c12.save.2[i]=cos(theta.diff.1[i] * (pi / 180)) * (sqrt(c11^2 + c12^2))
  c21.save.2[i]=sin(theta.diff.2[i] * (pi / 180)) * (sqrt(c21^2 + c22^2))
  c22.save.2[i]=cos(theta.diff.2[i] * (pi / 180)) * (sqrt(c21^2 + c22^2))
  
  # 2nd -- get new LV-components for the new angle setting
  K1.save.2[i] = (D * (S2.2 + S1.2 * Lamda.1 - B1)) / (c12.save.2[i] + c11.save.2[i] * Lamda.1)
  K2.save.2[i] = (D * (S2.2 + S1.2 * Lamda.2 - B2)) / (c22.save.2[i] + c21.save.2[i] * Lamda.2)
  alpha.save.2[i] = (c22.save.2[i] + c21.save.2[i] * Lamda.1) / (c12.save.2[i] + c11.save.2[i] * Lamda.1)
  beta.save.2[i]  = (c12.save.2[i] + c11.save.2[i] * Lamda.2) / (c22.save.2[i] + c21.save.2[i] * Lamda.2)
  
  # 3rd -- get new LV-ZNGI for the new angle setting
  LV.slope.1.save.2[i] = -(1 / alpha.save.2[i])
  LV.slope.2.save.2[i] = -(beta.save.2[i])
  LV.inter.1.save.2[i] = K1.save.2[i] / alpha.save.2[i]
  LV.inter.2.save.2[i] = K2.save.2[i]
  
  # 4th -- get impact.niche.overlap for the new angle setting
  status.save.2[i] = ifelse(theta.diff.2[i] > theta.diff.1[i], 0, 1)  
  impact.overlap.save.2[i] = (c11.save.2[i] * c21.save.2[i] + c12.save.2[i] * c22.save.2[i]) / (sqrt(c11.save.2[i]^2 + c12.save.2[i]^2) * sqrt(c21.save.2[i]^2 + c22.save.2[i]^2))
}




######################################
#### Plot ZNGI state-space diagram
######################################

#### Gather data for ploting 
moving.impact = data.frame(c11.save, 
                           c21.save, 
                           c22.save, 
                           c12.save, 
                           R1.star = rep(R1.star, length(c11.save)), 
                           R2.star = rep(R2.star, length(c11.save)))
ZNGI.df = data.frame(bl = c(R1.star), ob = c(R2.star))
supply.seg = data.frame(x1 = S1.2-0.125*(0.3-0.248), y1 = S2.2-0.125*(0.38-0.245), 
                        x2 = S1+0.11*(0.3-0.248), y2 = S2+0.11*(0.38-0.245))

#### Plot setting 
impact.subst.ZNGI = list()
impact.params = c(1, 70, 140, 140)
alpha.vec = c(1.0, 0.8, 0.5, 0.5)
curve.vec = c(0.5, 0.5, -0.5, -0.5)
pos.vec = c(0.008, 0.02, 0.048, 0.048)
label.vec = c("theta[1]", "theta[2]", "theta[3]", "theta[3]")
angle.sep = c(29, 29, 23, 23)
mylims.x = 0.55
mylims.y = 0.55


#### Plot !
for(i in 1:4){
  
  temp = eval(substitute(
      
      # Plot ZNGIs
      ggplot(ZNGI.df,
             aes(x=bl, y=ob)) +   
      geom_point(x = S1.2, y = S2.2, size = 3) +
      geom_abline(intercept = y.inter.1, 
                  slope = slope.1, 
                  col ='#d1495b', 
                  size = 0.5) +
      geom_abline(intercept = y.inter.2, 
                  slope = slope.2, 
                  col ='#30638e', 
                  size = 0.5) + 
      coord_cartesian(expand = 0, 
                      ylim=c(0, mylims.y), 
                      xlim = c(0, mylims.x)) +
      
      # Plot impact vectors and its extensions for selected angles 
      geom_segment(data = moving.impact[impact.params[i], ], 
                   aes(x = R1.star, y = R2.star, xend = R1.star-c11.save/100, yend = R2.star-c12.save/100), 
                   size = 0.5, 
                   col='#d1495b', 
                   arrow = arrow(type = "closed", length = unit(0.05, "inches"))) + 
      geom_segment(data = moving.impact[impact.params[i], ], 
                   aes(x = R1.star, y = R2.star, xend = R1.star-c21.save/100, yend = R2.star-c22.save/100), 
                   size = 0.5, col='#30638e', 
                   arrow = arrow(type = "closed", length = unit(0.05, "inches"))) + 
      geom_segment(data = moving.impact[impact.params[i], ], 
                   aes(x = R1.star, y = R2.star, xend = R1.star+c11.save, yend = R2.star+c12.save), 
                   size = 0.5, 
                   col='#d1495b', 
                   linetype = 2) + 
      geom_segment(data = moving.impact[impact.params[i],], 
                   aes(x = R1.star, y = R2.star, xend = R1.star+c21.save, yend = R2.star+c22.save), 
                   size = 0.5, 
                   col='#30638e', 
                   linetype = 2) + 
      
      # Plot selected angles with curves
      geom_curve(x = R1.star+c21.save[impact.params[i]]/angle.sep[i],
                 xend = R1.star+c11.save[impact.params[i]]/angle.sep[i],
                 y = R2.star+c22.save[impact.params[i]]/angle.sep[i],
                 yend = R2.star+c12.save[impact.params[i]]/angle.sep[i],
                 curvature = curve.vec[i],
                 size = 0.5,
                 col='black',
                 linetype = 1) +
      
      # Plot text indicating angles
      geom_text(x = R1.star+c21.save[impact.params[i]]/angle.sep[i] + pos.vec[i], 
                y = R2.star+c12.save[impact.params[i]]/angle.sep[i] + pos.vec[i], 
                label = label.vec[i], 
                parse = TRUE,  
                size = 5.5) +
      
      # Other settings 
      xlab(expression(R[1])) + 
      ylab(expression(R[2])) +
      
      theme(legend.position = "none", 
            plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"),
            axis.text = element_text(size=13),
            axis.title=element_text(size=20)) +
      panel_border(colour = "black")
      
  , list(i = i)))
 
  # Add Set 1 supply point only for the last panel and an arrow between the two supply points
  if(i == 4) {temp = temp + geom_point(x = S1, y = S2, size = 3, shape = 15)}
  if(i == 4) {temp = temp + geom_segment(data = supply.seg, aes(x = x1, y = y1, xend = x2, yend = y2), size = 0.5, arrow = arrow(length = unit(0.03, "npc"), angle = 15, type = "closed"))}
  
  impact.subst.ZNGI[[i]] = temp
}


#### Summarize into a four panel figure 
ZNGIs = plot_grid(impact.subst.ZNGI[[1]], impact.subst.ZNGI[[2]], impact.subst.ZNGI[[3]], impact.subst.ZNGI[[4]], nrow = 2, ncol = 2, labels = c("(a)", "(b)", "(c)", "(d)"), label_size=20)




######################################
#### Plot Chesson's plot
######################################

#### Function setting for Chesson's plot
ggplot.MCT.stablizing.logscale = function(YMIN){
  x.1 = seq(1, 0, by=-0.001)
  x.2 = seq(0, -1, by=-0.001)
  y1.1 = 1/((1-x.1))
  y2.1 = 1-x.1
  y1.2 = 1/((1-x.2))
  y2.2 = 1-x.2
  
  cf.df.1 = data.frame(my.x = c(x.1, x.1), 
                       my.y = c(y1.1, y2.1), 
                       whichy = c(rep("y1", times = length(seq(1, 0, by=-0.001))), 
                                  rep("y2", times = length(seq(1, 0, by=-0.001)))))
  rib.dims.1 = data.frame(min.dim = y1.1, max.dim = y2.1, x.dim = x.1)
  
  cf.df.2 = data.frame(my.x = c(x.2, x.2), 
                       my.y = c(y1.2, y2.2), 
                       whichy = c(rep("y1", times = length(seq(0, -1, by=-0.001))), 
                                  rep("y2", times = length(seq(0, -1, by=-0.001)))))
  rib.dims.2 = data.frame(min.dim = y1.2, max.dim = y2.2, x.dim = x.2)
  
  ggplot() + 
    geom_ribbon(data = rib.dims.1, 
                aes(x = x.dim, 
                    ymin = min.dim, 
                    ymax = max.dim), 
                fill = "grey", 
                alpha = 0.3) +
    geom_line(data = cf.df.1, 
              aes(x = my.x, 
                  y = my.y,
                  linetype = whichy,
                  group = whichy), 
              col = "black") +     
    geom_ribbon(data = rib.dims.2, 
                aes(x = x.dim, 
                    ymin = min.dim, 
                    ymax = max.dim), 
                fill = "darkgrey", 
                alpha = 0.6) +
    geom_line(data = cf.df.2, 
              aes(x = my.x, 
                  y = my.y,
                  linetype = whichy,
                  group = whichy), 
              col = "black") + 
    geom_abline(intercept = 0, 
                slope = 0, 
                lty = 2, 
                col = "darkgrey") + 
    geom_vline(xintercept = 0, 
               lty = 2, 
               col = "darkgrey") +
    panel_border(colour = "black") +
    scale_y_log10() +
    scale_linetype_manual(values = c("y2" = "solid", 
                                     "y1" = "dotted")) +
    xlab(expression(paste("Stabilization potential ( 1 - ", rho, " )", sep = ""))) + 
    ylab(expression(paste("Fitness ratio ( ", frac(italic(f[2]), italic(f[1])), " )", sep=""))) +
    theme(legend.position = "none", 
          axis.title.y = element_text(angle = 90), 
          axis.title = element_text(size = 20))
}


#### Setup ploting range with the index 'ri'
#### Set 1 (square, a supply point with S1 > S2)
ri = c(1:(99*2))
K1.plot.region = K1.save[ri]
K2.plot.region = K2.save[ri]
alpha.plot.region = alpha.save[ri]
beta.plot.region = beta.save[ri]
overlap.plot.region = impact.overlap.save[ri]

#### Set 2 (circle, a supply point with S1 < S2)
ri.2 = c(1:(99*2))
K1.plot.region.2 = K1.save.2[ri.2]
K2.plot.region.2 = K2.save.2[ri.2]
alpha.plot.region.2 = alpha.save.2[ri.2]
beta.plot.region.2 = beta.save.2[ri.2]
overlap.plot.region.2 = impact.overlap.save.2[ri.2]


#### Calculate Chesson's fitness and niche difference for each angle
#### Set 1 (square, a supply point with S1 > S2)
niche.overlap = sqrt(alpha.save[ri] * beta.save[ri])
fitness.diffs = (K2.save[ri] / K1.save[ri]) * sqrt(alpha.save[ri] / beta.save[ri])

#### Set 2 (circle, a supply point with S1 < S2)
niche.overlap.2 = sqrt(alpha.save.2[ri.2] * beta.save.2[ri.2])
fitness.diffs.2 = (K2.save.2[ri.2] / K1.save.2[ri.2]) * sqrt(alpha.save.2[ri.2] / beta.save.2[ri.2])


#### Gather data for ploting 
#### Set 1 (square, a supply point with S1 > S2)
df.temp.subst = data.frame(fitness = fitness.diffs, 
                           overlap = overlap.plot.region, 
                           Rho = niche.overlap, 
                           Rho.inv = 1/niche.overlap)
df.gg.subst = gather(df.temp.subst, key = param, value = score, -overlap)
stable.temp = data.frame(AFD = df.gg.subst$score[df.gg.subst$param == "fitness"], 
                         ND = df.gg.subst$score[df.gg.subst$param == "Rho"], 
                         SND = 1 - (df.gg.subst$score[df.gg.subst$param == "Rho"]))

#### Set 2 (circle, a supply point with S1 < S2)
df.temp.subst.2 = data.frame(fitness = fitness.diffs.2, 
                             overlap = overlap.plot.region.2, 
                             Rho = niche.overlap.2, 
                             Rho.inv = 1/niche.overlap.2)
df.gg.subst.2 = gather(df.temp.subst.2, key = param, value = score, -overlap)
stable.temp.2 = data.frame(AFD = df.gg.subst.2$score[df.gg.subst.2$param == "fitness"], 
                           ND = df.gg.subst.2$score[df.gg.subst.2$param == "Rho"], 
                           SND = 1 - (df.gg.subst.2$score[df.gg.subst.2$param == "Rho"]))

#### An arrow between the two supply points
supply.seg.chesson = data.frame(x1 = 1 - (niche.overlap.2[impact.params][4]), y1 = (fitness.diffs.2[4]+0.01),
                                x2 = 1 - (niche.overlap[impact.params][4]), y2 = (fitness.diffs[4])-0.01)


#### Plot setting
impact.params = c(1, 70, 140, 140)
alpha.vec = c(1.0, 0.8, 0.5, 0.5)
curve.vec = c(0.5, 0.5, -0.5, -0.5)
pos.vec = c(0.008, 0.02, 0.048, 0.048)
label.vec = c("theta[1]", "theta[2]", "theta[3]", "theta[3]")
angle.sep = c(29, 29, 23, 23)
mylims.x.Chesson = 0.5
mylims.y.Chesson = 1.7


#### Plot !
impact.subst.Chesson = 
  
  # Plot Chesson framework 
  ggplot.MCT.stablizing.logscale(2) + 
  
  # Plot simulation lines -- Set 2
  geom_line(data = stable.temp.2[1:max(impact.params), ], 
            aes(x=SND, 
                y=AFD), 
            size=0.5,
            color='#454545') +
  
  # Plot angle points -- Set 1 (only the last angle)
  geom_point(aes(x = 1 - (niche.overlap[impact.params][4]), 
                 y = (fitness.diffs[4])), 
             size = 3, 
             shape = 15) + 
  
  # Plot angle points -- Set 2
  geom_point(aes(x = 1 - (niche.overlap.2[impact.params][1]), 
                 y = (fitness.diffs.2[1])), 
             size = 3) +
  geom_point(aes(x = 1 - (niche.overlap.2[impact.params][2]), 
                 y = (fitness.diffs.2[2])), 
             size = 3) +
  geom_point(aes(x = 1 - (niche.overlap.2[impact.params][3]), 
                 y = (fitness.diffs.2[3])), 
             size = 3) + 
  
  # Plot text indicating angles -- Set 1 (only the last angle)
  geom_text(aes(x = 1 - (niche.overlap[impact.params][4]), 
                y = (fitness.diffs[4]) + 0.03), 
            label = "theta[3]", 
            parse = TRUE,  
            size = 5.5) + 

  # Plot text indicating angles -- Set 2
  geom_text(aes(x = 1 - (niche.overlap.2[impact.params][1]), 
                y = (fitness.diffs.2[1]) - 0.03), 
            label = "theta[1]", 
            parse = TRUE,  
            size = 5.5) +
  geom_text(aes(x = 1 - (niche.overlap.2[impact.params][2]), 
                y = (fitness.diffs.2[2]) - 0.03), 
            label = "theta[2]", 
            parse = TRUE,  
            size = 5.5) +
  geom_text(aes(x = 1 - (niche.overlap.2[impact.params][3]), 
                y = (fitness.diffs.2[3]) - 0.03), 
            label = "theta[3]", 
            parse = TRUE,  
            size = 5.5) + 

  # Add arrow between the two supply points
  geom_segment(data = supply.seg.chesson, 
               aes(x = x1, y = y1, xend = x2, yend = y2), 
               size = 0.5, 
               # linetype = 3, 
               arrow = arrow(length = unit(0.02, "npc"),
                             angle = 15, 
                             type = "closed")) + 
  
  # Add text for coexistence/priority effect region:
  geom_text(aes(x = 0.0,
                y = exp(log(1.7)/2)),
            label = "Species 2 wins",
            size = 6.0, 
            fontface = "bold") +  
  geom_text(aes(x = 0.0,
                y = exp(-log(1.7)/2)),
            label = "Species 1 wins",
            size = 6.0, 
            fontface = "bold") +  
  geom_text(aes(x = 0.29,
                y = exp(log(1.7)/18)),
            label = "Coexistence",
            size = 6.0, 
            fontface = "bold") +  
  geom_text(aes(x = -0.29,
                y = exp(log(1.7)/18)),
            label = "Priority effect",
            size = 6.0, 
            fontface = "bold") +
  
  # Set ticks for the log-scaled y-axis
  scale_y_log10(breaks = c(0.6, 1, 1.6), 
                labels = c("0.6", "1", "1.6")) +
  
  # Other settings
  coord_cartesian(expand = c(0, 0), 
                  xlim = c(-mylims.x.Chesson, mylims.x.Chesson), 
                  ylim = c(1/mylims.y.Chesson, mylims.y.Chesson)) +
  theme(legend.position = "none", 
        plot.margin = unit(c(0.8,0.8,0.8,0.8), "lines"), 
        axis.text = element_text(size=13),
        axis.title = element_text(size=20)) +
  panel_border(colour = "black")




######################################
#### Save plots
######################################
dev.off()
CairoPDF(file = "Fig1.pdf", width = 18, height = 9)
plot_grid(ZNGIs, impact.subst.Chesson, labels = c("", "(e)"), nrow = 1, label_size=20)
dev.off()


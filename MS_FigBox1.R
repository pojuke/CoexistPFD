#################################################################################
#### Coexistence theory and the frequency-dependence of priority effects 
#### Ke & Letten (2018) Nature Ecology & Evolution
#### This R script creates the figure in Box 1
#### Coexistence and priority effects in a Lotka-Volterra competition model
#################################################################################




######################################
#### Load packages
######################################
library(tidyr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(grid)
library(gridExtra)




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


#### Plot setting
mylims.x.Chesson = 0.5
mylims.y.Chesson = 1.7


#### Plot !
Chesson.conceptual =
  
  # Plot Chesson framework 
  ggplot.MCT.stablizing.logscale(2) + 
  
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
  
  # Add point within coexistence/priority effect region:
  geom_point(aes(x = 1 - sqrt((0.8*0.96)/(1.4*1.4)), 
                 y = sqrt((1.4*0.96)/(1.4*0.8))), 
             size = 3, shape = 15) + 
  geom_point(aes(x = 1 - sqrt((0.8*0.96)/(0.64*0.64)), 
                 y = sqrt((0.64*0.96)/(0.64*0.8))), 
             size = 3, shape = 16) + 

  # Add mathematical expression
  annotate("text", 
           x=0.29, y=exp(-log(1.7)/8), size=6.0,
           label='paste(italic(rho), " < ", italic(frac(f[2], f[1])), " < ", italic(frac(1, rho)))', 
           parse=T) + 
  annotate("text", 
           x=-0.29, y=exp(-log(1.7)/8), size=6.0,
           label='paste(italic(frac(1, rho)), " < ", italic(frac(f[2], f[1])), " < ", italic(rho))', 
           parse=T) + 
  
  # Other settings
  coord_cartesian(expand = c(0, 0), 
                  xlim = c(-mylims.x.Chesson, mylims.x.Chesson), 
                  ylim = c(1/mylims.y.Chesson, mylims.y.Chesson)) +
  theme(legend.position = "none", 
        # plot.margin = unit(c(0.8, 0.8, 5.0, 0.8), "lines"), 
        axis.text = element_text(size=13),
        axis.title = element_text(size=20)) 

Chesson.conceptual = ggdraw() +
                     draw_plot(Chesson.conceptual, x = 0.25, y = 0, width = .5, height = 1) +
                     draw_plot_label(label = "(b)", size = 18,
                                    x = 0.25, y = 1)




######################################
#### Plot time series for parameter sets
######################################
require('deSolve')

#### The model
LV=function(Time, State, Pars)
{
  with(as.list(c(State, Pars)),
       {
         dN1 = r1*N1*(1 - a11*N1 - a12*N2)
         dN2 = r2*N2*(1 - a21*N1 - a22*N2)
         return(list(c(dN1, dN2)))
       })
}

## The general parameter setting
r1 = 0.1
r2 = 0.1
a11.vec = c(1.4, 0.64)
a21.vec = c(0.8, 0.8)
a12.vec = c(0.96, 0.96)
a22.vec = c(1.4, 0.64)

## The simulation setting (time step and initial condition)
times = seq(0, 300, by = 0.01)
yini.1more = c(N1=1, N2=0.2)
yini.2more = c(N1=0.2, N2=1)

## The saving space for plot 
label.vec = c("(a)", "(c)") 
time.series.plots = list()

## The simulation for coexistence
i = 1
pars.temp = c(r1=0.1, r2=0.1, 
              a11=a11.vec[i], 
              a21=a12.vec[i], 
              a12=a12.vec[i],
              a22=a22.vec[i])
  
out.1more.co = as.data.frame(ode(func = LV, y = yini.1more, parms = pars.temp, times = times, method = "ode45"))
out.2more.co = as.data.frame(ode(func = LV, y = yini.2more, parms = pars.temp, times = times, method = "ode45"))
DATA.1more.co = out.1more.co[seq(0, dim(out.1more.co)[1], by=100) + 1, ]
DATA.2more.co = out.2more.co[seq(0, dim(out.2more.co)[1], by=100) + 1, ]
DATA.1more.co = gather(DATA.1more.co, key=Variable, value=PopSize, -time)
DATA.2more.co = gather(DATA.2more.co, key=Variable, value=PopSize, -time)
  
# Plot for coexistence
temp.1more.co = ggplot(DATA.1more.co[DATA.1more.co$Variable %in% c("N1", "N2"), ], aes(x=time, y=PopSize, color=Variable)) + 
                geom_line(size=1.2) +
                scale_color_manual(values=c("N1"="#d1495b", "N2"="#30638e")) + 
                labs(x=" ", y="Density") +
                # scale_x_continuous(expand=c(0, 5)) +
                scale_y_continuous(limits=c(0, 1.6)) +
                theme(legend.position = "none", 
                      axis.text = element_text(size=13),
                      axis.title = element_text(size=20)) 
            
temp.2more.co = ggplot(DATA.2more.co[DATA.2more.co$Variable %in% c("N1", "N2"), ], aes(x=time, y=PopSize, color=Variable)) + 
                geom_line(size=1.2) +
                scale_color_manual(values=c("N1"="#d1495b", "N2"="#30638e")) + 
                labs(x="Time", y="Density") +
                # scale_x_continuous(expand=c(0, 5)) +
                scale_y_continuous(limits=c(0, 1.6)) +
                theme(legend.position = "none", 
                      axis.text = element_text(size=13),
                      axis.title = element_text(size=20)) 

Coexist = ggdraw() +
          draw_plot(temp.1more.co, x = 0, y = 0.5, width = 1, height = 0.5) +
          draw_plot(temp.2more.co, x = 0, y = 0, width = 1, height = 0.5) +
          draw_plot_label(label = c("(c)", ""), size = 18,
                          x = c(0, 0), y = c(1, 0.5))


## The simulation for priority effect
i = 2
pars.temp = c(r1=0.1, r2=0.1, 
              a11=a11.vec[i], 
              a21=a12.vec[i], 
              a12=a12.vec[i],
              a22=a22.vec[i])
  
out.1more.pe = as.data.frame(ode(func = LV, y = yini.1more, parms = pars.temp, times = times, method = "ode45"))
out.2more.pe = as.data.frame(ode(func = LV, y = yini.2more, parms = pars.temp, times = times, method = "ode45"))
DATA.1more.pe = out.1more.pe[seq(0, dim(out.1more.pe)[1], by=100) + 1, ]
DATA.2more.pe = out.2more.pe[seq(0, dim(out.2more.pe)[1], by=100) + 1, ]
DATA.1more.pe = gather(DATA.1more.pe, key=Variable, value=PopSize, -time)
DATA.2more.pe = gather(DATA.2more.pe, key=Variable, value=PopSize, -time)
  
# Plot for coexistence
temp.1more.pe = ggplot(DATA.1more.pe[DATA.1more.pe$Variable %in% c("N1", "N2"), ], aes(x=time, y=PopSize, color=Variable)) + 
                geom_line(size=1.2) +
                scale_color_manual(values=c("N1"="#d1495b", "N2"="#30638e")) + 
                labs(x=" ", y="Density") +
                # scale_x_continuous(expand=c(0, 5)) +
                scale_y_continuous(limits=c(0, 1.6)) +
                theme(legend.position = "none", 
                      axis.text = element_text(size=13),
                      axis.title = element_text(size=20)) 
            
temp.2more.pe = ggplot(DATA.2more.pe[DATA.2more.pe$Variable %in% c("N1", "N2"), ], aes(x=time, y=PopSize, color=Variable)) + 
                geom_line(size=1.2) +
                scale_color_manual(values=c("N1"="#d1495b", "N2"="#30638e")) + 
                labs(x="Time", y="Density") +
                # scale_x_continuous(expand=c(0, 5)) +
                scale_y_continuous(limits=c(0, 1.6)) +
                theme(legend.position = "none", 
                      axis.text = element_text(size=13),
                      axis.title = element_text(size=20)) 

Priority = ggdraw() +
           draw_plot(temp.1more.pe, x = 0, y = 0.5, width = 1, height = 0.5) +
           draw_plot(temp.2more.pe, x = 0, y = 0, width = 1, height = 0.5) +
           draw_plot_label(label = c("(a)", ""), size = 18,
                          x = c(0, 0), y = c(1, 0.5))

   


######################################
#### Save plots
######################################
dev.off()
CairoPDF(file = "Conceptual_a.pdf", width = 4, height = 9)
Priority
dev.off()


dev.off()
CairoPDF(file = "Conceptual_b.pdf", width = 17, height = 8)
Chesson.conceptual
dev.off()


dev.off()
CairoPDF(file = "Conceptual_c.pdf", width = 4, height = 9)
Coexist
dev.off()



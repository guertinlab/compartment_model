
data.df = read.table("Nilson2017_Kinit_fit_and_deriv.txt", sep="\t", header=T)
ggplot(data.df) + geom_point(aes(x=time_sec, y = kinit_AU, color="Original Data"), size=2) +
         geom_line(aes(x=time_sec, y = fittedPoints, color="Fit"), size=1) +
      scale_y_continuous(
        limits = c(0, 1.2),  # Set lower and upper limits for the y-axis
        breaks = seq(0,1.2,0.2) # Set custom tick positions
       )     +
      scale_x_continuous(
        limits = c(0, 1800),  # Set lower and upper limits for the y-axis
        breaks = seq(0,1800,600) # Set custom tick positions
       )     +
         geom_vline(xintercept = 750, linetype=2) +
         xlab("Time (sec)") +
         ylab("Fraction of Initiation") +
         ggtitle("Fit for Fig. 6B, Nilson et al., 2017") +
         theme_bw() + fontTheme + 
         theme(legend.position="none",
               plot.title = element_text(size=24, face=1.2,hjust=0.5))

ggsave("Nilson2017_Init_vs_time.pdf", width=9.13)
#plot derivative
ggplot(data.df) + 
        geom_point(aes(x=time_sec, y = -1 * derivative, color="Derivative"), size=2, color = "blue") +
         geom_line(aes(x=time_sec, y = -1 * derivative, color="Derivative"), color= "green") +
      #geom_point(aes(x=time_sec, y =  derivative, color="Derivative"), size=2, color = "blue") +
       # geom_line(aes(x=time_sec, y = derivative, color="Derivative"), color= "green") +
        #geom_hline(yintercept = 0, linetype=2) +
         geom_vline(xintercept = 750, linetype=2) +
          scale_y_continuous(trans=log10_trans(),
                breaks=trans_breaks('log10', function(x) 10^x),
                labels=trans_format('log10', math_format(.x)), limits=c(NA,1)) +
 scale_x_continuous(
        limits = c(0, 1800),  # Set lower and upper limits for the y-axis
        breaks = seq(0,1800,600) # Set custom tick positions
       ) + 
         xlab("Time (sec)") +
         #ylab(expression("(Slope of Fit)")) +
                  ggtitle("Fit for Fig. 6B, Nilson et al., 2017") +
         ylab(expression(log[10]~"(-1 x Slope of Fit)")) +
         theme_bw() + fontTheme + 
          theme(plot.title = element_text(size=24, face=1.2,hjust=0.5))
ggsave("Nilson2017_Slope_Init_vs_time.pdf")
ggsave("Nilson2017_log10_Slope_Init_vs_time.pdf", width=8.6)


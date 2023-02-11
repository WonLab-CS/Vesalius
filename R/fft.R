################################################################################
################################   Vesalius      ###############################
################################################################################

#------------------------------/Fourier Transform/-----------------------------#

# internal_fft_reduction <- function(signal,
#     partials = 100) {
#     #--------------------------------------------------------------------------#
#     # Fun fast fourier transform and get revsersed reduced 
#     # TO ADD: input checks
#     #--------------------------------------------------------------------------#
#     signal <- vesalius:::min_max(signal)
#     transformed <- fft(signal)
#     rev <- Re(fft(transformed[seq(1, partials)], inverse = TRUE))
#     rev <- vesalius:::min_max(rev)
#     return(list("reversed" = rev,
#         "reduced_signal" = signal[seq(1, length(signal), l = partials)]))

# }

# partials <- 350
# test <- internal_fft_reduction(sort(vect), partials = partials)
# ylim <- c(min(min(test$reversed), min(test$reduced_signal)),max(max(test$reversed), max(test$reduced_signal)))

# plot(0, type = "n", xlim = c(0, length(vect)), ylim = ylim)
# lines(seq_along(vect), vesalius:::min_max(sort(vect)), lty = 1)
# lines(seq(1,length(vect), l = partials), test$reversed, lty = 1, col = "red")
# lines(seq(1,length(vect), l = partials),test$reduced_signal, lty = 2, col ="blue")





fake_actig <- data.frame(patient_ID = rep(1:3, each = 14*720),
                         Noon_Day = rep(rep(1:14, each = 720), 3),
                         Epoch = rep(seq(1, 720 * 14, 1), 3),
                         Hour = rep(rep(1:24, each = 30), 3 * 14)) %>%
  mutate(Sine = sin(Epoch / round(720 / (2 * pi))),
         Sign = sign(Sine),
         Noise = rnorm(Epoch))

## ------ Test ------
IV_IS_test <- list()

IV_IS_combs <- expand.grid(c("Sine", "Sign", "Noise"), c("moving", "expanding"),
                           c(3, 7)) %>% 
  rename(var = Var1, window = Var2, size = Var3) %>%
  arrange(var, window, size)

for (j in 1L:nrow(IV_IS_combs)) {
  var <- IV_IS_combs$var[j] %>% as.character
  window <- IV_IS_combs$window[j] %>% as.character
  size <- IV_IS_combs$size[j]
  
  IV_title <- paste(var, "IV", window, size, sep = "_")
  IS_title <- paste(var, "IS", window, size, sep = "_")
  
  IV <- lapply(unique(fake_actig$patient_ID), function(ii) {
    calculate_IV(var, ii, window, size, fake_actig)
  }) %>%
    do.call("rbind", .) %>%
    rename_(.dots=setNames(names(.), gsub("IV", IV_title, names(.))))
  
  IS <- lapply(unique(fake_actig$patient_ID), function(ii) {
    calculate_IS(var, ii, window, size, fake_actig)
  }) %>%
    do.call("rbind", .) %>%
    rename_(.dots=setNames(names(.), gsub("IS", IS_title, names(.))))
  
  IV_IS_test[[IV_title]] <- IV
  IV_IS_test[[IS_title]] <- IS
}

IV_IS_test_results <- IV_IS_test[1][[1]]

for (i in 2L:length(IV_IS_test)) {
  IV_IS_test_results <- merge(IV_IS_test_results, IV_IS_test[i][[1]],
        by = c("patient_ID", "Day"))
}

IV_IS_test_results %<>%
  arrange(patient_ID, Day) %>%
  mutate(patient_ID = factor(patient_ID))

## ------ Plot ------

longitudinal_plot <- function(var) {
  p <- ggplot(IV_IS_test_results, aes_string("Day", var, color = "patient_ID")) +
    geom_point(size = 2) +
    geom_smooth() +
    ylim(0, 2)
  
  return(p)
}

test_IV_IS_names <- data.frame(beep = rep(c("IV", "IS"), each = nrow(IV_IS_combs)),
                               rbind(IV_IS_combs, IV_IS_combs)) %>%
  mutate(nam = paste(var, beep, window, size, sep = "_")) %>%
  select(nam) %>%
  unlist

## Plot List
plot_list <- list()
for (j in 1:length(test_IV_IS_names)) {
  plot_list[[j]] <- longitudinal_plot(test_IV_IS_names[j])
}

pdf("Tests/test_IV_IS_plots.pdf")

for (j in 1:length(plot_list)) {
  print(plot_list[j])
}

dev.off()


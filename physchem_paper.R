library(tidyverse)
library(here)
library(readxl)
library(patchwork)
library(tools)

my_theme <-  theme_grey() +
  theme(axis.text = element_text(size = 18)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(axis.title = element_text(size = 22)) +
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_text(size = 20)) +
  theme(strip.text = element_text(face = "bold", size = 15)) +
  theme(legend.position = "right")
theme_set(my_theme)

# Description------------------------------------------------------------------

# Samples "Oct22": 
# series of solvothermal syntheses from urea or thiourea and citric acid
# variable ratio of components
# same temperature, duration, concentration

# material intended for publication on PhysChem



# Functions--------------------------------------------------------------------

# reads a single excel sheet obtained from Varioskan microplate reader
read_sheet <- function(data_file, sheet){
  suppressMessages(
    {
      first_col <- read_excel(path = data_file, 
                          sheet = sheet,
                          col_names = FALSE,
                          range = cell_cols("A"))[[1]]
    }
    )
  
  wl_range <- str_split(first_col[6], ',')[[1]]
  
  exp_type <- ifelse(grepl("Em", wl_range[1]), "emission", "excitation")
  fixed_wl <- as.numeric(gsub(".*?([0-9]+).*", "\\1", wl_range[2]))
  
  first_line <- which(first_col == "Wavelength")
  num_rows <- min(which(is.na(first_col[first_line:length(first_col)])))
  
  num_cols <- ncol(read_excel(path = data_file, 
                              sheet = sheet,
                              col_names = TRUE,
                              range = cell_rows(first_line)))
  
  spectra <- read_excel(path = data_file,
                        sheet = sheet,
                        col_names = TRUE,
                        range = anchored("A11", dim = c(num_rows - 1, num_cols)))
  return(list(spectra = spectra, exp_type = exp_type, fixed_wl = fixed_wl))
}

# reads spectral data from multipage excel file 
# obtained from Varioskan microplate reader
# returns data frame in long format
tidy_data <- function(spectra_lst){
  data <- spectra_lst$spectra %>% 
    pivot_longer(cols = -c("Wavelength"),
                 names_to = "cell",
                 values_to = "intensity") %>%
    rename(var_wavelength = Wavelength) %>%
    mutate(fixed_wavelength = spectra_lst$fixed_wl,
           exp_type = spectra_lst$exp_type,
           cell = str_extract(cell, "(?<=\\().+?(?=\\))"))
  
  return(data)
}

# reads samples data from excel file
# adding samples data to spectral data frame
# returns data frame in long format
read_plate <- function(file_name){
  data_file <- here("data", file_name)
  
  data_sheets <- excel_sheets(path = data_file) %>% 
    .[matches("^data_", vars = .)]
  
  samples <- read_excel(path = data_file, 
                        sheet = "samples",
                        col_names = TRUE)
  
  data <- map(data_sheets, read_sheet, data_file = data_file) %>% 
    map(tidy_data) %>%
    bind_rows() %>%
    left_join(samples, by = "cell")
  
  return(data)
}

# reads a csv file with emission spectrum from Ocean Optics spectrometer
# returns spectral data frame with two columns: 
# wavelength (nm) and intensity (a.u.)
# note: decimal delimiter is ','
read_ocean_fluoro <- function(filename){
  con <- file(description = filename, open = "r", blocking = TRUE)
  n <- 0
  chunk <- ''
  
  while (chunk != '>>>>>Begin Spectral Data<<<<<') {
    n <- n + 1
    chunk <- readLines(con = con, n = 1)
  }
  close(con)
  
  data <- read_delim(file = filename,
                     delim = "\t",
                     locale = locale(decimal_mark = ","),
                     skip = n,
                     col_names = c("wavelength", "intensity"),
                     show_col_types = FALSE)
  data$file <- tools::file_path_sans_ext(basename(filename))
  return(data)
}


# reads the csv file with absorbance spectrum from Ocean Optics spectrometer
# returns spectral data frame with two columns: wavelength (nm) and absorbance
# note: decimal delimiter is ','
read_ocean_abs <- function(filename){
  con <- file(description = filename, open = "r", blocking = TRUE)
  n <- 0
  chunk <- ''
  
  while (chunk != '>>>>>Begin Spectral Data<<<<<') {
    n <- n + 1
    chunk <- readLines(con = con, n = 1)
  }
  close(con)
  
  data <- read_delim(file = filename,
                     delim = "\t",
                     locale = locale(decimal_mark = ","),
                     skip = n,
                     col_names = c("wavelength", "absorbance"),
                     show_col_types = FALSE)
  data$file <- tools::file_path_sans_ext(basename(filename))
  return(data)
}



bl_correction <- function(spectrum, wl_range){
  offset <- spectrum %>%
    filter(wavelength > wl_range[1] & wavelength < wl_range[2]) %>%
    .[[2]] %>%
    mean()
  spectrum[, 2] <- spectrum[, 2] - offset
  return(spectrum)
}


smooth_runmed <- function(spectrum, k = 3) {
  spectrum[, 2] = runmed(x = spectrum[[2]], k = k)
  return(spectrum)
}


integral <- function(spectrum, wl_range){
  intensity <- spectrum %>%
    filter(wavelength > wl_range[1] & wavelength < wl_range[2]) %>%
    .[['intensity']] %>%
    sum()
  return(tibble(integral = intensity / (wl_range[2] - wl_range[1])))
}



# Load data--------------------------------------------------------------------

# example of particles size distribution
U_CA_3_1_size <- read_delim(file = here("data", "urea_3_1_distribution_clean.txt"),
                            delim = "\t",
                            locale = locale(decimal_mark = "."),
                            skip = 1,
                            col_names = c("size", "n", "concentration", "volume", "area"),
                            show_col_types = FALSE,
                            col_select = c(size, n)) %>%
  mutate(cum_n = cumsum(n))




# 1st plate reader experiment, Nov 09, 2022
# urea-citric acid nanoparticles in water
# varied: composition, dilution
plate_1 <- read_plate("cnd_vario_Oct22_1.xlsx")

# 2nd plate reader experiment, Nov 14, 2022
# urea-citric acid nanoparticles in water, Hg
# varied: composition, dilution, presence of Hg
plate_2 <- read_plate("cnd_vario_Oct22_2.xlsx")

# 3rd plate reader experiment, Nov 21, 2022
# maps of thioure:CA emission
# sensitivity of thiourea:CA samples to Hg
# detection limit of Hg with urea:CA samples
plate_3 <- read_plate("cnd_vario_Oct22_3.xlsx")


# reading absorbance data from Ocean Optics
ocean_abs_raw <- here("data", "QY_Oct22", "abs") %>%
  list.files(pattern = "*.txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_abs) 
  

# reading absorbance data from Ocean Optics
ocean_fluoro_raw <- here("data", "QY_Oct22", "fluoro") %>%
  list.files(pattern = "*.txt", full.names = TRUE) %>%
  lapply(FUN = read_ocean_fluoro) 




# Data preparation-------------------------------------------------------------

# combined data from plates
# smoothed spectrum-wise
plate <- bind_rows(list(plate_1, plate_2, plate_3)) %>%
  group_by(plate, cell, exp_type, fixed_wavelength) %>%
  arrange(var_wavelength, .by_group = TRUE) %>%
  mutate(intensity = runmed(x = intensity, k = 5)) %>%
  ungroup()

# baseline correction and smoothing of absorbance spectra from Ocean Optics
# parsing filename into sample metadata
# removing spectral ranges containing artifacts
ocean_abs <- ocean_abs_raw %>%
  lapply(FUN = bl_correction, wl_range = c(800, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  separate(col = file,
           sep = '_',
           into = c("type", NA, "comp1", "comp2", "sample"),
           fill = "right") %>%
  mutate(sample = case_when(sample == "p" ~ "purified",
                            sample == "d"~ "dialysate")) %>%
  mutate(type = case_when(type == "t" ~ "thiourea",
                          type == "u"~ "urea",
                          type == "ref" ~ "reference")) %>%
  mutate(composition = paste(comp1, comp2, sep = ":")) %>%
  filter(wavelength > 260 & wavelength < 900) %>%
  filter(!(wavelength > 655 & wavelength < 659)) %>%
  select(-c("comp1", "comp2")) %>%
  mutate(composition = factor(composition, levels = c("7:1", "5:1", "3:1", "1:1", "1:3", "1:5", "1:7")))


# baseline correction and smoothing of emission spectra from Ocean Optics
# parsing filename into sample metadata
# removing spectral ranges containing artifacts
ocean_fluoro <- ocean_fluoro_raw %>%
  lapply(FUN = bl_correction, wl_range = c(800, 900)) %>%
  lapply(FUN = smooth_runmed, k = 5) %>%
  bind_rows() %>%
  separate(col = file,
           sep = '_',
           into = c("type", NA, "comp1", "comp2", "sample"),
           fill = "right") %>%
  mutate(sample = case_when(sample == "p" ~ "purified",
                            sample == "d"~ "dialysate")) %>%
  mutate(type = case_when(type == "t" ~ "thiourea",
                          type == "u"~ "urea",
                          type == "ref" ~ "reference")) %>%
  mutate(composition = paste(comp1, comp2, sep = ":")) %>%
  filter(wavelength > 300 & wavelength < 900)  %>% 
  select(-c("comp1", "comp2")) %>%
  mutate(composition = factor(composition, levels = c("7:1", "5:1", "3:1", "1:1", "1:3", "1:5", "1:7")))

# calculation of quantum yield
abs_365 <- ocean_abs %>%
  group_by(type, sample, composition) %>%
  filter(wavelength > 360 & wavelength < 370) %>%
  summarize(absorbance = mean(absorbance))

integrals <- ocean_fluoro %>%
  group_by(type, sample, composition) %>%
  group_modify(~integral(spectrum = .x, wl_range = c(390, 800)), .keep = TRUE)

ref <- c(abs = abs_365$absorbance[which(abs_365$type == "reference")],
         int = integrals$integral[which(integrals$type == "reference")])

QY <- abs_365 %>%
  left_join(integrals, by = c("type", "sample", "composition")) %>%
  mutate(qy_pct = 54 * integral / ref["int"] * ref["abs"] / absorbance) %>%
  mutate(qy_pct = round(qy_pct, digits = 1))


# Final figures for paper------------------------------------------------------
# particles size distribution
# urea:CA = 3:1
factor <- max(U_CA_3_1_size$n) / sum(U_CA_3_1_size$n)

fig_3 <- ggplot(data = U_CA_3_1_size, aes(x = size)) +
  geom_line(aes(y = cum_n * max(n) / sum(n)), size = 2, color = "red") +
  geom_line(aes(y = n), size = 2) +
  scale_x_log10(limits = c(10, 1000), name = "Size/ nm") +
  scale_y_continuous(name = 'Number of particles', limits = c(0, max(U_CA_3_1_size$n)), 
                     sec.axis = sec_axis(trans = ~ . / max(U_CA_3_1_size$n) * 100, 
                                         name = "Cummulative fraction, %\n"))
ggsave(here("plots", "fig_3.tiff"), plot = fig_3)  


# absorbance spectra of the purified samples
fig_4a <- ggplot(data = ocean_abs %>% 
                   filter(type == "urea" & sample == "purified")) +
  geom_line(aes(x = wavelength, 
                y = absorbance, 
                color = composition), 
            size = 2) +
  labs(x = "Wavelength, nm", y = "Absorbance", color = "urea/citric acid") +
  theme(legend.position = c(0.75, 0.7)) +
  scale_color_brewer(palette = "Dark2")
ggsave(here("plots", "fig_4a.tiff"), plot = fig_4a)  

fig_4b <- ggplot(data = ocean_abs %>% 
                   filter(type == "thiourea" & sample == "purified")) +
  geom_line(aes(x = wavelength, 
                y = absorbance, 
                color = composition), 
            size = 2) +
  labs(x = "Wavelength, nm", y = "Absorbance", color = "thiourea/citric acid") +
  theme(legend.position = c(0.75, 0.7)) +
  scale_color_brewer(palette = "Dark2")
ggsave(here("plots", "fig_4b.tiff"), plot = fig_4b)

# comparison of absorbance spectra of dialysates and purified samples
fig_5a <- ggplot(data = ocean_abs %>% 
                   filter(type == "urea" & composition == "1:1")) +
  geom_line(aes(x = wavelength, y = absorbance, color = sample), size = 2) +
  labs(x = "Wavelength, nm", y = "Absorbance", color = "Sample") +
  theme(legend.position = c(0.85, 0.8)) +
  scale_color_brewer(palette = "Dark2")
ggsave(here("plots", "fig_5a.tiff"), plot = fig_5a)  

fig_5b <- ggplot(data = ocean_abs %>% 
                   filter(type == "thiourea" & composition == "1:1")) +
  geom_line(aes(x = wavelength, y = absorbance, color = sample), size = 2) +
  labs(x = "Wavelength, nm", y = "Absorbance", color = "Sample") +
  theme(legend.position = c(0.85, 0.8)) +
  scale_color_brewer(palette = "Dark2")
ggsave(here("plots", "fig_5b.tiff"), plot = fig_5b)  


# emission spectra of the purified samples (ex 365 nm)
fig_6a <- ggplot(data = ocean_fluoro %>% 
                   filter(type == "urea" & sample == "purified") %>%
                   filter(wavelength > 400)) +
  geom_line(aes(x = wavelength, 
                y = intensity, 
                color = composition), 
            size = 2) +
  labs(x = "Wavelength, nm", y = "Intensity", color = "urea/citric acid") +
  theme(legend.position = c(0.75, 0.7)) +
  scale_color_brewer(palette = "Dark2")
ggsave(here("plots", "fig_6a.tiff"), plot = fig_6a)  

fig_6b <- ggplot(data = ocean_fluoro %>% 
                   filter(type == "thiourea" & sample == "purified") %>%
                   filter(wavelength > 400)) +
  geom_line(aes(x = wavelength, 
                y = intensity, 
                color = composition), 
            size = 2) +
  labs(x = "Wavelength, nm", y = "Intensity", color = "thiourea/citric acid") +
  theme(legend.position = c(0.75, 0.7)) +
  scale_color_brewer(palette = "Dark2")
ggsave(here("plots", "fig_6b.tiff"), plot = fig_6b)

# comparison of emission spectra of dialysates and purified samples
fig_6c <- ggplot(data = ocean_fluoro %>% 
                   filter(type == "urea" & composition == "1:1") %>%
                   filter(wavelength > 400)) +
  geom_line(aes(x = wavelength, 
                y = intensity, 
                color = sample), size = 2) +
  labs(x = "Wavelength, nm", y = "Intensity", color = "Sample") +
  theme(legend.position = c(0.85, 0.8)) +
  scale_color_brewer(palette = "Dark2")
ggsave(here("plots", "fig_6c.tiff"), plot = fig_6c)  

fig_6d <- ggplot(data = ocean_fluoro %>% 
                   filter(type == "thiourea" & composition == "1:1") %>%
                   filter(wavelength > 400)) +
  geom_line(aes(x = wavelength, y = intensity, color = sample), size = 2) +
  labs(x = "Wavelength, nm", y = "Intensity", color = "Sample") +
  theme(legend.position = c(0.85, 0.8)) +
  scale_color_brewer(palette = "Dark2")
ggsave(here("plots", "fig_6d.tiff"), plot = fig_6d)  


# emission at different dilutions urea:CA 1:1
fig_7 <- ggplot(data = plate %>%
                  filter(plate == 2) %>%
                  filter(exp_type == "emission") %>%
                  filter(type == "urea") %>%
                  filter(dilution != 2) %>%
                  filter(composition == "1_1" & hg == 0) %>%
                  filter(var_wavelength <= 700)) +
  geom_contour_filled(aes(x = var_wavelength, y = fixed_wavelength, z = intensity), bins = 20) +
  labs(x = "Emission, nm", y = "Excitation, nm") +
  facet_wrap(~factor(dilution), ncol = 3) +
  theme(legend.position = "none")
ggsave(here("plots", "fig_7.tiff"), plot = fig_7)


# comparison of emission maps for different urea:CA samples 
urea_dilutions <- data.frame(composition = c("7_1", "5_1", "3_1", "1_1", "1_3", "1_5", "1_7"),
                             dilution = c(2, 5, 10, 20, 5, 2, 2))

urea_maps <- plate %>%
  filter(plate == 2) %>%
  filter(exp_type == "emission") %>%
  filter(type == "urea") %>%
  filter(hg == 0) %>%
  filter(var_wavelength <= 700) %>%
  merge(urea_dilutions) %>%
  group_by(composition, dilution) %>%
  mutate(intensity = intensity / max(intensity))


fig_8 <- ggplot(data = urea_maps) +
  geom_contour_filled(aes(x = var_wavelength, y = fixed_wavelength, z = intensity), bins = 20) +
  labs(x = "Emission, nm", y = "Excitation, nm") +
  facet_wrap(~factor(composition, 
                     levels = c("7_1", "5_1", "3_1", "1_1", "1_3",  "1_5", "1_7"),
                     labels = c("7:1", "5:1", "3:1", "1:1", "1:3",  "1:5", "1:7")), 
             ncol = 4) +
  theme(legend.position = "none")
ggsave(here("plots", "fig_8.tiff"), plot = fig_8)


# comparison of emission maps for different thiourea:CA samples 
thiourea_dilutions <- data.frame(composition = c("7_1", "5_1", "3_1", "1_1", "1_3", "1_5"),
                             dilution = c(2, 2, 2, 2, 2, 2))

thiourea_maps <- plate %>%
  filter(plate == 3) %>%
  filter(exp_type == "emission") %>%
  filter(type == "thiourea") %>%
  filter(hg == 0) %>%
  filter(var_wavelength <= 700) %>%
  filter(fixed_wavelength != 380) %>%
  merge(thiourea_dilutions) %>%
  group_by(composition, dilution) %>%
  mutate(intensity = intensity / max(intensity))


fig_9 <- ggplot(data = thiourea_maps) +
  geom_contour_filled(aes(x = var_wavelength, y = fixed_wavelength, z = intensity), bins = 20) +
  labs(x = "Emission, nm", y = "Excitation, nm") +
  facet_wrap(~factor(composition, 
                     levels = c("7_1", "5_1", "3_1", "1_1", "1_3",  "1_5"),
                     labels = c("7:1", "5:1", "3:1", "1:1", "1:3",  "1:5")), 
             ncol = 3) +
  theme(legend.position = "none")
ggsave(here("plots", "fig_9.tiff"), plot = fig_9)


# effect of Hg on emission of urea:CA samples of selected composition
urea_hg_dilutions <- data.frame(composition = c("7_1", "1_1", "1_7"),
                                dilution = c(2, 20, 2))

urea_hg_maps <- plate %>%
  filter(plate == 2) %>%
  filter(exp_type == "emission") %>%
  filter(type == "urea") %>%
  filter(var_wavelength <= 700) %>%
  merge(urea_hg_dilutions) %>%
  group_by(composition) %>%
  mutate(intensity = intensity / max(intensity)) %>%
  mutate(composition = factor(composition,
                              levels = c("7_1", "1_1", "1_7"),
                              labels = c("7:1", "1:1", "1:7"))) %>%
  mutate(hg = factor(hg,
                     levels = c(0, 1),
                     labels = c("no Hg", "with Hg")))

fig_10 <- ggplot(data = urea_hg_maps) +
  geom_contour_filled(aes(x = var_wavelength, y = fixed_wavelength, z = intensity), bins = 20) +
  labs(x = "Emission, nm", y = "Excitation, nm") +
  facet_grid(composition ~ hg) +
  theme(legend.position = "none")
ggsave(here("plots", "fig_10.tiff"), plot = fig_10)




# estimation of detection limit
detection_wl_u71 <- data.frame(composition = c("7_1", "7_1", "7_1"),
                               fixed_wavelength = c(540, 440, 360),
                               var_wavelength = c(630, 550, 460))

detection_u71 <- plate %>%
  filter(plate == 3) %>%
  filter(type == "urea") %>%
  filter(exp_type == "emission") %>%
  filter(var_wavelength <= 700) %>%
  merge(detection_wl_u71) %>%
  group_by(fixed_wavelength) %>%
  mutate(rel_intensity = intensity / max(intensity)) %>%
  filter(hg_conc > 0) %>%
  mutate(fixed_wavelength = factor(fixed_wavelength,
                                   levels = c(360, 440, 540),
                                   labels = c("360 / 460", "440 / 550", "540 / 630")))

fig_11a <- ggplot(data = detection_u71, aes(x = hg_conc, y = rel_intensity, color = factor(fixed_wavelength))) +
  geom_point(size = 3) +
  geom_smooth(se = F, method = 'loess', size = 2) +
  scale_x_log10(limits = c(2e-8, 2e-4)) +
  theme(legend.position = c(0.85, 0.8)) +
  labs(x = "Hg(II) concentration, mol/L", y = "Relative intensity", color = "Ex / Em, nm")
ggsave(here("plots", "fig_11a.tiff"), plot = fig_11a)





detection_wl <- data.frame(composition = c("7_1", "1_1", "1_3", "1_7"),
                           fixed_wavelength = c(560, 460,460, 340),
                           var_wavelength = c(630, 550, 550, 420))

detection <- plate %>%
  filter(plate == 3) %>%
  filter(type == "urea") %>%
  filter(exp_type == "emission") %>%
  filter(var_wavelength <= 700) %>%
  merge(detection_wl) %>%
  group_by(composition) %>%
  mutate(rel_intensity = intensity / max(intensity)) %>%
  filter(hg_conc > 0) %>%
  mutate(composition = factor(composition,
                              levels = c("7_1", "1_1", "1_3", "1_7"),
                              labels = c("7:1", "1:1", "1:3", "1:7")))

fig_11b <- ggplot(data = detection, aes(x = hg_conc, y = rel_intensity, color = composition)) +
  geom_point(size = 3) +
  geom_smooth(se = F, method = 'loess', size = 2) +
  scale_x_log10(limits = c(2e-8, 5e-4)) +
  theme(legend.position = c(0.86, 0.8)) +
  labs(x = "Hg(II) concentration, mol/L", y = "Relative intensity", color = "urea/citric acid")
ggsave(here("plots", "fig_11b.tiff"), plot = fig_11b)




# effect of Hg on emission of thiourea:CA samples of selected composition
thiourea_hg_dilutions <- data.frame(composition = c("3_1", "1_3"),
                                    dilution = c(2, 2))

thiourea_hg_maps <- plate %>%
  filter(plate == 3) %>%
  filter(exp_type == "emission") %>%
  filter(type == "thiourea") %>%
  filter(var_wavelength <= 700) %>%
  filter(fixed_wavelength != 380) %>%
  merge(thiourea_hg_dilutions) %>%
  group_by(composition) %>%
  mutate(intensity = intensity / max(intensity)) %>%
  mutate(composition = factor(composition,
                              levels = c("3_1", "1_3"),
                              labels = c("3:1", "1:3"))) %>%
  mutate(hg = factor(hg,
                     levels = c(0, 1),
                     labels = c("no Hg", "with Hg")))

fig_12 <- ggplot(data = thiourea_hg_maps) +
  geom_contour_filled(aes(x = var_wavelength, y = fixed_wavelength, z = intensity), bins = 20) +
  labs(x = "Emission, nm", y = "Excitation, nm") +
  facet_grid(composition ~ hg) +
  theme(legend.position = "none")
ggsave(here("plots", "fig_12.tiff"), plot = fig_12)
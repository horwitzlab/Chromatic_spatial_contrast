# Chromatic spatial contrast
This repository contains supplementary information and Matlab code accompanying the manuscript by De and Horwitz, "Coding of chromatic spatial contrast by macaque V1 neurons"

# monitor details
 Matlab structure containing the following fields:
    background_rgb: intensity of the red, green, and blue primaries to the background of the display in normalized units (0 = minimum intensity, 1 = maximum intensity).
    wavelengths: lattice of wavelengths at which irradiance was measured (in nm).
    red_phosphor_spd: spectral power distribution of the red monitor phosphor measured at maximum intensity. Irradiance is measured in "watts/sr•m^2•bin” for each 5 nm-wide bin.
    green_phosphor_spd: as above, but for the green phosphor.
    blue_phosphor_spd: as above, but for the blue phosphor.
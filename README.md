# CASTOR

CASTOR Plots

I've written a script for the CASTOR Plots based on Dr. Matt's email. I've gotten the code up and running and it's giving me a plot and the fluxes—but the flux for each band seems way too low (it is many orders below the total flux or total luminosity). I've included the script in the Python file and some typical plots with values for each band below. I can't seem to trace any physics/math errors in my code.

The temperature I used in this script is a median temperature (about 32,000K) from the 500 temperatures I calculated using the formula in the Arnett paper coupled with the blackbody equation, given the inputs of masses ranging from 8–20 solar masses, at an average value of 750 solar radii. I've varied the radius of the progenitor from 8 to 20 solar masses (typical range). The plots seem to mostly be identical in shape but the flux values vary. I've attached them below.

---

### 8M☉

![8M☉](https://github.com/user-attachments/assets/9402a7d2-e902-4645-ace3-39fb7e891605)

| Band | Flux (W/m²) | Luminosity (W) |
|------|-------------|----------------|
| UV   | 1.361e+05   | 4.659e+29      |
| U    | 1.123e+04   | 3.846e+28      |
| G    | 4.000e+03   | 1.370e+28      |

**Bolometric flux:** 5.657e+10 W/m²  
**Bolometric luminosity:** 1.937e+35 W

---

### 12M☉

![12M☉](https://github.com/user-attachments/assets/f96fffc3-088f-437d-a5c0-67c7a4706f1f)

| Band | Flux (W/m²) | Luminosity (W) |
|------|-------------|----------------|
| UV   | 1.229e+05   | 4.210e+29      |
| U    | 1.015e+04   | 3.475e+28      |
| G    | 3.614e+03   | 1.238e+28      |

**Bolometric flux:** 5.657e+10 W/m²  
**Bolometric luminosity:** 1.937e+35 W

---

### 16M☉

![16M☉](https://github.com/user-attachments/assets/7415a844-ea25-4373-b32b-aa67c05ac6e4)

| Band | Flux (W/m²) | Luminosity (W) |
|------|-------------|----------------|
| UV   | 1.144e+05   | 3.918e+29      |
| U    | 9.445e+03   | 3.234e+28      |
| G    | 3.363e+03   | 1.152e+28      |

**Bolometric flux:** 5.657e+10 W/m²  
**Bolometric luminosity:** 1.937e+35 W

---

### 20M☉

![20M☉](https://github.com/user-attachments/assets/a2bd86dd-ba38-47e3-8ae6-f6cc093484a7)

| Band | Flux (W/m²) | Luminosity (W) |
|------|-------------|----------------|
| UV   | 1.082e+05   | 3.705e+29      |
| U    | 8.933e+03   | 3.059e+28      |
| G    | 3.181e+03   | 1.089e+28      |

**Bolometric flux:** 5.657e+10 W/m²  
**Bolometric luminosity:** 1.937e+35 W

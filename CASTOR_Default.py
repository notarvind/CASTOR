#script with default values of 10M☉ mass, 1e51 erg/s explosion energy, and 1e13 cm radius, as per Dr. Matt's email

import numpy as np
import matplotlib.pyplot as plt

c = 3.0e8
E_exp = 10**51 #ergs
kappa = 0.4 #photospheric opacity in cm^2/g
K_0 = 5.2* 10**43 #avg K value, in paper
radius = 10**13 #cm
mass = 10 #solar masses
sigma = 5.67*10**-8
h = 6.62607015*10**-34 #plancks constant
k = 1.380649 * 10**-23 #boltzmann constant

lum2 = K_0*(radius/(10**14))*(E_exp/10**51)*(0.4/kappa)*(1/mass) #using the Arnett's paper formula

radius_m = radius/100 #converting to m

lum2_watts = lum2 * 1e-7 #watts

temp_default = (lum2_watts/(4*np.pi* (radius_m**2) * sigma))**(1/4) #blackbody

#initialising bands
uv_band = np.linspace(150, 300, 1000)*1e-9
u_band = np.linspace(300, 400, 1000)*1e-9
g_band = np.linspace(400, 500, 1000)*1e-9
all_bands = np.linspace(150, 500, 1000)*1e-9

#planck function
def planck_function(wavelengths, temperatures):
    exponent = (h * c) / (wavelengths * k * temperatures)
    B = (2.0 * h * c**2) / (wavelengths**5) / (np.exp(exponent) - 1)
    return B

#flux calculation
B_lambda_uv = planck_function(uv_band, temp_default)
B_lambda_u = planck_function(u_band, temp_default)
B_lambda_g = planck_function(g_band, temp_default)
B_lambda_all = planck_function(all_bands, temp_default)

plt.figure(figsize=(8, 5))
plt.plot(uv_band, B_lambda_uv, color="purple", label="uv band")
plt.plot(u_band, B_lambda_u, color="blue", label="u band")
plt.plot(g_band, B_lambda_g, color="green", label="g band")
plt.legend()
plt.xlabel("Wavelengths of uv band (m)")
plt.ylabel("Spectral Radiance $B_\\lambda$ (W/m²·sr·m)")
plt.title("Blackbody Spectrum across CASTOR Bands at 10M☉, 1e51 erg/s, 1e13 cm")
plt.grid(True)
plt.show()

#assume isotropic - flux same in all directions (sphere)

total_flux_density_uv = B_lambda_uv * np.pi
total_flux_density_u = B_lambda_u * np.pi
total_flux_density_g = B_lambda_g * np.pi

#integrating
uv_flux = np.trapezoid(total_flux_density_uv, uv_band)
u_flux = np.trapezoid(total_flux_density_u, u_band)
g_flux = np.trapezoid(total_flux_density_g, g_band)

print("Spectral uv radiance: ", B_lambda_uv[250])
print(f"Total uv flux: {uv_flux:.3e} W/m²")
print(f"Total u flux: {u_flux:.3e} W/m²")
print(f"Total g flux: {g_flux:.3e} W/m²")

#bolometric flux from Stefan-Boltzmann
F_bol = sigma * temp_default**4
print(f"Bolometric flux: {F_bol:.3e} W/m²")
print(f"UV flux fraction: {uv_flux / F_bol:.3%}")


#total luminosities
L_uv = uv_flux * 4 * np.pi * radius_m**2
L_u = u_flux * 4 * np.pi * radius_m**2
L_g = g_flux * 4 * np.pi * radius_m**2
L_bol = lum2_watts

print("Peak wavelength: ", (2.898*10**-3)/temp_default) #wiens displacement law 
print("Luminosity in watts: ", lum2_watts)
print("Temperature: ", temp_default)
print(f"UV Luminosity: {L_uv:.3e} W")
print(f"U Luminosity: {L_u:.3e} W")
print(f"G Luminosity: {L_g:.3e} W")
print(f"Bolometric Luminosity: {L_bol:.3e} W")
print(f"UV % of Bolometric: {L_uv / L_bol:.3%}")


#table
print("\n=== Luminosity Summary ===")
print(f"{'UV':<10} {uv_flux:<15.3e} {L_uv:<20.3e} {(L_uv/L_bol)*100:<10.2f}")
print(f"{'U':<10} {u_flux:<15.3e} {L_u:<20.3e} {(L_u/L_bol)*100:<10.2f}")
print(f"{'G':<10} {g_flux:<15.3e} {L_g:<20.3e} {(L_g/L_bol)*100:<10.2f}")

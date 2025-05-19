import numpy as np
import matplotlib.pyplot as plt

c = 3.0e8
E_exp = 10**51 #ergs
kappa = 0.4 #photospheric opacity in cm^2/g
K_0 = 5.2* 10**43 #avg K value, in paper
low_mass = 8 #lowest mass, in solar masses  ---> this is what i vary.
avg_radius = 750 #avg solar radii
sigma = 5.67*10**-8
h = 6.62607015*10**-34 #plancks constant
k = 1.380649 * 10**-23 #boltzmann constant


R_sun_cm = 6.96e10  #solar radius, cm
radii = np.linspace(500, 1000, num=500) #radii in solar radii
masses = np.linspace(8, 20, num=500) #masses in solar masses

lum1 = np.zeros(500)
lum2 = np.zeros(500)

radii_cm = radii * R_sun_cm
avg_radius_cm = avg_radius * R_sun_cm

for a in range(500):
    lum1[a] = K_0*(radii_cm[a]/(10**14))*(E_exp/10**51)*(0.4/kappa)*(1/low_mass)

temp_mconstant = np.zeros(500)
lum1_watts = lum1 * 1e-7
radii_m = radii_cm / 100  #convert to m

for b in range(500):
    temp_mconstant[b] = (lum1_watts[b]/(4*np.pi* (radii_m[b]**2) * sigma))**(1/4)

for c in range(500):
    lum2[c] = K_0*(avg_radius_cm/(10**14))*(E_exp/10**51)*(0.4/kappa)*(1/masses[c])

temp_rconstant = np.zeros(500)
lum2_watts = lum2 * 1e-7
avg_radius_m = avg_radius_cm / 100  #convert to m

for d in range(500):
    temp_rconstant[d] = (lum2_watts[d]/(4*np.pi* (avg_radius_m**2) * sigma))**(1/4)

uv_band = np.linspace(150, 300)*1e-9
u_band = np.linspace(300, 400)*1e-9
g_band = np.linspace(400, 500)*1e-9
all_bands = np.linspace(150, 500)*1e-9

def planck_function(wavelengths, temperatures):
    exponent = (h * c) / (wavelengths * k * temperatures)
    B = (2.0 * h * c**2) / (wavelengths**5) / (np.exp(exponent) - 1)
    return B

B_lambda_uv = planck_function(uv_band, temp_mconstant[250])
B_lambda_u = planck_function(u_band, temp_mconstant[250])
B_lambda_g = planck_function(g_band, temp_mconstant[250])
B_lambda_all = planck_function(all_bands, temp_mconstant[250])

plt.figure(figsize=(8, 5))
plt.plot(uv_band, B_lambda_uv, color="purple", label="uv band")
plt.plot(u_band, B_lambda_u, color="blue", label="u band")
plt.plot(g_band, B_lambda_g, color="green", label="g band")
plt.legend()
plt.xlabel("Wavelengths of uv band (nm)")
plt.ylabel("Spectral Radiance $B_\\lambda$ (W/m²·sr·m)")
plt.title("Blackbody Spectrum across CASTOR Bands at Median Temperature at 8 M☉")
plt.grid(True)
plt.show()

#range of peak wavelengths:
print((2.898*10**-3)/temp_rconstant[0]) #wiens displacement law min
print((2.898*10**-3)/temp_rconstant[499]) #wiens displacement law max

#assume isotropic - flux same in all directions (sphere)

total_flux_density_uv = B_lambda_uv * np.pi
total_flux_density_u = B_lambda_u * np.pi
total_flux_density_g = B_lambda_g * np.pi

#integrating
uv_flux = np.trapezoid(total_flux_density_uv, uv_band)
u_flux = np.trapezoid(total_flux_density_u, u_band)
g_flux = np.trapezoid(total_flux_density_g, g_band)

print(f"Total uv flux: {uv_flux:.3e} W/m²")
print(f"Total u flux: {u_flux:.3e} W/m²")
print(f"Total g flux: {g_flux:.3e} W/m²")


#bolometric flux from Stefan-Boltzmann
F_bol = sigma * temp_rconstant[250]**4
print(f"Bolometric flux: {F_bol:.3e} W/m²")
print(f"UV flux fraction: {uv_flux / F_bol:.3%}")


#total luminosities
L_uv = uv_flux * 4 * np.pi * avg_radius_m**2
L_u = u_flux * 4 * np.pi * avg_radius_m**2
L_g = g_flux * 4 * np.pi * avg_radius_m**2
L_bol = lum2_watts[250]

print(f"UV Luminosity: {L_uv:.3e} W")
print(f"U Luminosity: {L_u:.3e} W")
print(f"G Luminosity: {L_g:.3e} W")
print(f"Bolometric Luminosity: {L_bol:.3e} W")
print(f"UV % of Bolometric: {L_uv / L_bol:.3%}")

#table
print("\n=== Luminosity Summary ===")
print(f"{'Band':<10} {'Flux (W/m²)':<15} {'Luminosity (W)':<20} {'% of L_bol':<10}")
print(f"{'UV':<10} {uv_flux:<15.3e} {L_uv:<20.3e} {(L_uv/L_bol)*100:<10.2f}")
print(f"{'U':<10} {u_flux:<15.3e} {L_u:<20.3e} {(L_u/L_bol)*100:<10.2f}")
print(f"{'G':<10} {g_flux:<15.3e} {L_g:<20.3e} {(L_g/L_bol)*100:<10.2f}")





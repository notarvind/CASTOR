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


lum2 = K_0*(radius/(10**14))*(E_exp/10**51)*(0.4/kappa)*(1/mass) 
radius_m = radius/100 #converting to m
lum2_watts = lum2 * 1e-7 #watts


temp_default = (lum2_watts/(4*np.pi* (radius_m**2) * sigma))**(1/4) 


uv_band = np.linspace(150, 300, 1000)*1e-9
u_band = np.linspace(300, 400, 1000)*1e-9
g_band = np.linspace(400, 500, 1000)*1e-9
all_bands = np.linspace(150, 500, 1000)*1e-9

#planck function
def planck_function(wavelengths, temperatures):
    temperatures = np.maximum(temperatures, 1e-10) 
    exponent = (h * c) / (wavelengths * k * temperatures)
    B = (2.0 * h * c**2) / (wavelengths**5) / (np.exp(exponent) - 1)
    return B


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
plt.ylabel("Spectral Flux Density $B_\\lambda$ (W/m²·sr·m)")
plt.title("Blackbody Spectrum across CASTOR Bands with at 20M☉, 1e51 ergs, 1e13 cm")
plt.grid(True)
plt.show()

#assume isotropic - flux same in all directions 
total_flux_density_uv = B_lambda_uv * np.pi 
total_flux_density_u = B_lambda_u * np.pi
total_flux_density_g = B_lambda_g * np.pi

#integrating
uv_flux = np.trapz(total_flux_density_uv, uv_band)
u_flux = np.trapz(total_flux_density_u, u_band)
g_flux = np.trapz(total_flux_density_g, g_band)

print(f"Spectral uv radiance: {B_lambda_uv[250]:.3e}")
print(f"Total uv flux: {uv_flux:.3e} W/m²")
print(f"Total u flux: {u_flux:.3e} W/m²")
print(f"Total g flux: {g_flux:.3e} W/m²")

#bolometric flux from Stefan-Boltzmann
F_bol = sigma * temp_default**4
print(f"Bolometric flux: {F_bol:.3e} W/m²")
print(f"UV flux fraction: {uv_flux / F_bol:.3%}")

#total luminosities in bands
L_uv = uv_flux * 4 * np.pi * radius_m**2
L_u = u_flux * 4 * np.pi * radius_m**2
L_g = g_flux * 4 * np.pi * radius_m**2
L_bol = lum2_watts

print(f"Peak wavelength: {(2.898*10**-3)/temp_default:.3e} m") #wiens displacement law 
print(f"Luminosity in watts: {lum2_watts:.3e}")
print(f"Temperature: {temp_default:.2f} K")
print(f"UV Luminosity: {L_uv:.3e} W")
print(f"U Luminosity: {L_u:.3e} W")
print(f"G Luminosity: {L_g:.3e} W")
print(f"Bolometric Luminosity: {L_bol:.3e} W")
print(f"UV % of Bolometric: {L_uv / L_bol:.3%}")


print("\n=== Luminosity Summary ===")
print(f"{'Band':<10} {'Flux (W/m²)':<15} {'Luminosity (W)':<20} {'% of Bolometric':<15}")
print(f"{'UV':<10} {uv_flux:<15.3e} {L_uv:<20.3e} {(L_uv/L_bol)*100:<15.2f}")
print(f"{'U':<10} {u_flux:<15.3e} {L_u:<20.3e} {(L_u/L_bol)*100:<15.2f}")
print(f"{'G':<10} {g_flux:<15.3e} {L_g:<20.3e} {(L_g/L_bol)*100:<15.2f}")

#Mass vs Radius plot
mass_range = np.linspace(5, 20, 50)
radius_range = np.logspace(12, 14, 50) 
radius_range_m = radius_range/100 #converting to m
energy_range = np.logspace(50, 52, 50) 

mass_grid, radius_grid = np.meshgrid(mass_range, radius_range)
total_mr_lum_grid = np.zeros_like(mass_grid, dtype=float)

for a in range(mass_grid.shape[0]):
    for b in range(mass_grid.shape[1]):
        mass_val = mass_grid[a, b]
        radius_val = radius_grid[a, b] 
        radius_m_mr = radius_val/100 
        
        lum1 = K_0*(radius_val/(10**14))*(E_exp/10**51)*(0.4/kappa)*(1/mass_val)
        lum1_watts = lum1 * 1e-7
        
        temp_mr = (lum1_watts/(4*np.pi* (radius_m_mr**2) * sigma))**(1/4)
        
        if not np.isfinite(temp_mr) or temp_mr <= 0:
            total_mr_lum_grid[a, b] = 0
            continue

        B_lambda_uv_1 = planck_function(uv_band, temp_mr)*np.pi
        B_lambda_u_1 = planck_function(u_band, temp_mr)*np.pi
        B_lambda_g_1 = planck_function(g_band, temp_mr)*np.pi
        
        L_uv_mr = np.trapz(B_lambda_uv_1, uv_band)*4 * np.pi * radius_m_mr**2
        L_u_mr = np.trapz(B_lambda_u_1, u_band)*4 * np.pi * radius_m_mr**2
        L_g_mr = np.trapz(B_lambda_g_1, g_band)*4 * np.pi * radius_m_mr**2
        
        total_mr_lum_grid[a, b] = L_u_mr + L_uv_mr + L_g_mr

normalised_lum_mr = total_mr_lum_grid / np.max(total_mr_lum_grid)    
plt.yscale("log")   
contour_mr = plt.contourf(mass_range, radius_range_m, normalised_lum_mr, levels=100, cmap="plasma")
cbar = plt.colorbar(contour_mr)     
cbar.set_label("Normalised Total Luminosity") 
plt.xlabel("Mass (M☉)")
plt.ylabel("Radius (m)")
plt.title("Normalised Total Luminosity in uv, u, and g Bands (Mass vs. Radius)")
plt.tight_layout()
plt.show()

#Radius vs Energy plot
radius_grid2, energy_grid = np.meshgrid(radius_range, energy_range) 
total_re_lum_grid = np.zeros_like(radius_grid2, dtype=float) 

for i in range(radius_grid2.shape[0]):
    for j in range(radius_grid2.shape[1]):
        energy_val = energy_grid[i, j]
        radius_val2 = radius_grid2[i, j] 
        radius_m_re = radius_val2/100 #in m
        
        lum3 = K_0*(radius_val2/(10**14))*(energy_val/10**51)*(0.4/kappa)*(1/mass) 
        lum3_watts = lum3 * 1e-7
        
        temp_re = (lum3_watts/(4*np.pi* (radius_m_re**2) * sigma))**(1/4)

        B_lambda_uv_2 = planck_function(uv_band, temp_re)*np.pi
        B_lambda_u_2 = planck_function(u_band, temp_re)*np.pi
        B_lambda_g_2 = planck_function(g_band, temp_re)*np.pi
        
        L_uv_re = np.trapz(B_lambda_uv_2, uv_band)*4 * np.pi * radius_m_re**2
        L_u_re = np.trapz(B_lambda_u_2, u_band)*4 * np.pi * radius_m_re**2
        L_g_re = np.trapz(B_lambda_g_2, g_band)*4 * np.pi * radius_m_re**2
        
        total_re_lum_grid[i, j] = L_u_re + L_uv_re + L_g_re


normalised_lum_re = total_re_lum_grid / np.max(total_re_lum_grid)
plt.xscale("log")    
plt.yscale("log")   
contour_re = plt.contourf(radius_range_m, energy_range, normalised_lum_re, levels=100, cmap="plasma")
cbar = plt.colorbar(contour_re)     
cbar.set_label("Normalised Total Luminosity") 
plt.ylabel("Energy (ergs)")
plt.xlabel("Radius (m)") 
plt.title("Normalised Total Luminosity in uv, u, and g Bands")
plt.tight_layout()
plt.show()

#Mass vs Energy plot
mass_grid2, energy_grid2 = np.meshgrid(mass_range, energy_range) 
total_me_lum_grid = np.zeros_like(mass_grid2, dtype=float)

for y in range(mass_grid2.shape[0]):
    for z in range(mass_grid2.shape[1]):
        mass_val2 = mass_grid[y, z]
        energy_val2 = energy_grid2[y, z]
        
        lum4 = K_0*(radius/(10**14))*(energy_val2/10**51)*(0.4/kappa)*(1/mass_val2) 
        lum4_watts = lum4 * 1e-7
        
        temp_me = (lum4_watts/(4*np.pi* ((radius/100)**2) * sigma))**(1/4)

        B_lambda_uv_3 = planck_function(uv_band, temp_me)*np.pi
        B_lambda_u_3 = planck_function(u_band, temp_me)*np.pi
        B_lambda_g_3 = planck_function(g_band, temp_me)*np.pi
        
        L_uv_me = np.trapz(B_lambda_uv_3, uv_band)*4 * np.pi * radius_m_re**2
        L_u_me = np.trapz(B_lambda_u_3, u_band)*4 * np.pi * radius_m_re**2
        L_g_me = np.trapz(B_lambda_g_3, g_band)*4 * np.pi * radius_m_re**2
        
        total_me_lum_grid[y, z] = L_uv_me + L_u_me + L_g_me


normalised_lum_me = total_me_lum_grid / np.max(total_me_lum_grid)
plt.yscale("log")   
contour_me = plt.contourf(mass_range, energy_range, normalised_lum_me, levels=100, cmap="plasma")
cbar = plt.colorbar(contour_me)     
cbar.set_label("Normalised Total Luminosity") 
plt.ylabel("Energy (ergs)")
plt.xlabel("Mass (M☉)")
plt.title("Normalised Total Luminosity in uv, u, and g Bands")
plt.tight_layout()
plt.show()



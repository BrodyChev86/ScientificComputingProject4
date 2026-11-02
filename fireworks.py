#!/usr/bin/env python3
"""
fireworks.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.widgets import CheckButtons



# Representative spectral lines in nanometers (used claude to list all the points in (wavelength nm, intensity) format)
# For the following elements, intensitiies are relative to the element's strongest line so different thresholds are used for each element
sodium = np.array([(411.370, 150),(432.462, 1),(439.003, 1),(439.281, 600),(439.334, 1),(440.512, 400),(445.523, 500),(449.087, 400),(449.418, 1),
                   (449.766, 1),(466.4811, 2),(466.8560, 2),(497.8541, 2),(498.2813, 5),(514.8838, 1),(515.3402, 1),(541.455, 50),(568.2633, 4),
                   (568.8193, 1),(568.8205, 7),(588.9950, 1000),(589.5924, 500),(615.4225, 2),(616.0747, 3),(651.421, 50),(653.070, 70),(654.404, 70),
                   (654.575, 70),]) #https://physics.nist.gov/PhysRefData/Handbook/Tables/sodiumtable2.htm

copper = np.array([(406.264, 200),(465.112, 200),(490.9734, 700),(493.1698, 600),(495.3724, 500),(505.1793, 600),(510.554, 150),(515.324, 200),
                      (521.820, 250),(529.252, 150),(570.024, 150),(578.213, 150),(600.0120, 400),(615.0384, 400),(615.4222, 500),(621.6939, 500),
                      (621.9844, 500),(627.3349, 700),(630.1009, 600),(637.7840, 500),(642.3884, 600),(644.8559, 500),(647.0168, 600),(648.1437, 500),
                      (662.4292, 500),(664.1396, 500),(740.4354, 700),]) #https://physics.nist.gov/PhysRefData/Handbook/Tables/coppertable2.htm 

strontium = np.array([(403.038, 20),(407.771, 1000),(416.180, 4),(421.552, 700),(430.545, 7),(460.733, 1000),(472.228, 50),(474.192, 30),(478.432, 20),
                         (481.188, 70),(483.208, 60),(487.249, 50),(487.632, 30),(496.226, 120),(496.794, 20),(522.220, 20),(522.511, 30),(522.927, 30),
                         (523.855, 40),(525.690, 70),(545.084, 25),(548.084, 110),(550.417, 50),(552.183, 40),(553.481, 30),(554.005, 30),(640.847, 140),
                         (650.400, 80),(655.026, 25),(661.726, 50),(679.105, 30),(687.838, 70),(689.259, 20),(707.010, 80),(730.941, 40)]) # https://physics.nist.gov/PhysRefData/Handbook/Tables/strontiumtable2.htm 

barium = np.array([(389.1779, 20),(413.0649, 25),(413.2427, 9),(416.6001, 3),(428.3097, 30),(452.4926, 2),(455.4033, 1000),(489.9927, 6),(493.4077, 300),
                   (553.5481, 1000),(577.7619, 80),(585.3675, 40),(599.7087, 80),(611.0783, 250),(614.1713, 300),(649.6898, 200),(649.8760, 250),(652.7312, 110),
                   (659.5325, 100),(667.5270, 50),(669.3842, 50),(686.5686, 9),(705.9943, 200),(712.0331, 30),(719.5230, 9),(728.0296, 150),(739.2405, 15),
                   (748.8075, 25),]) #https://physics.nist.gov/PhysRefData/Handbook/Tables/bariumtable2.htm 

lithium = np.array([(413.2624, 2),(413.2625, 2),(415.5946, 10),(432.5419, 15),(432.5471, 15),(460.2831, 15),(460.2898, 30),(460.2902, 30),(467.1651, 30),
                    (467.1705, 10),(467.8056, 80),(467.8290, 25),(478.8363, 80),(488.1320, 15),(488.1386, 15),(488.1490, 5),(497.1665, 5),(497.1748, 10),
                    (503.7915, 15),(548.356, 300),(548.450, 500),(548.511, 400),(610.3542, 300),(610.3654, 400),(610.3667, 400),(670.7775, 500),(670.7926, 1000),]) #https://physics.nist.gov/PhysRefData/Handbook/Tables/lithiumtable2.htm 

calcium = np.array([(393.36614, 1000),(396.84673, 1000),(422.6727, 1000),(430.2527, 500),(430.7741, 500),(442.5441, 500),(443.4960, 500),(443.5688, 500),
                    (445.4781, 600),(445.5887, 600),(445.6605, 400),(487.8132, 500),(501.9971, 300),(518.8848, 500),(526.5557, 500),(527.0270, 500),(534.9472, 500),
                    (558.1971, 500),(558.8757, 500),(559.0120, 500),(559.4468, 500),(559.8487, 500),(585.7452, 600),(610.2722, 500),(612.2219, 600),(616.2172, 600),
                    (616.9055, 500),(616.9559, 600),(643.9073, 700),(644.9810, 600),(645.6874, 300),(646.2566, 700),(647.1660, 600),(649.3780, 600),(649.9649, 600),
                    (657.2777, 500),(671.7685, 600),(714.8147, 700),(720.2194, 600),(732.6146, 700),]) #https://physics.nist.gov/PhysRefData/Handbook/Tables/calciumtable2.htm 

# ELEMENT DATABASE
elements = {"Sodium": sodium, "Copper": copper, "Strontium": strontium, "Barium": barium, "Lithium": lithium, "Calcium": calcium}

# Spectral Line ranges for compound salts
molecular_bands = {
    "Na": (584, 594), #Na expaned to make visible in plot normally 589
    "Li": (664, 674), #Li expanded to make visible in plot normally 671
    "SrCl": (636, 662), #SrCl
    "CaCl": (593, 622), #CaCl
    "BaCl": (507, 553), #BaCl
    "CuCl": (410, 500) #CuCl
}

element_bands = {
    "Sodium": "Na",
    "Lithium": "Li",
    "Strontium": "SrCl",
    "Calcium": "CaCl",
    "Barium": "BaCl",
    "Copper": "CuCl",
}



fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharey=True)

def plot_spectrum(ax, wavelengths, band_name, title):
    ax.set_title(title)
    ax.set_xlim(380, 720)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([380, 440, 500, 560, 620, 680, 720])
    ax.set_xticklabels(["380", "440", "500", "560", "620", "680", "720"], fontsize=5)
    ax.set_xlabel("Wavelength (nm)").set_fontsize(5)
    band = molecular_bands[band_name]
    ax.axvspan(band[0], band[1], color="red", alpha=0.2)
    for wavelength in wavelengths:
        ax.axvline(wavelength[0], color="black", alpha=wavelength[1]/1000)  
    
# Plots the spectra of the different elements 
plot_spectrum(axes[0, 0], sodium, "Na", "Sodium")
plot_spectrum(axes[0, 1], lithium, "Li", "Lithium")
plot_spectrum(axes[0, 2], strontium, "SrCl", "Strontium")
plot_spectrum(axes[1, 0], calcium, "CaCl", "Calcium")
plot_spectrum(axes[1, 1], barium, "BaCl", "Barium")
plot_spectrum(axes[1, 2], copper, "CuCl", "Copper")

# Plots the firework color
axes[2, 0].remove()
firework_ax = fig.add_subplot(3, 3, 7)
firework_ax.set_title("Firework color")
firework_ax.set_axis_off()
firework_circle = Circle((0, 0), 1, color="black")
firework_ax.add_patch(firework_circle)
firework_ax.set_xlim(-1.1, 1.1)
firework_ax.set_ylim(-1.1, 1.1)
firework_ax.set_aspect("equal")

# Plots the check buttons for selecting elements
check_ax = axes[2, 1]
check_ax.set_title("Select elements")
check_buttons = CheckButtons(check_ax, list(elements))


def wavelength_to_rgb(band):
    wavelength = (band[0] + band[1]) / 2
    if wavelength < 380 or wavelength > 780: #outside visible range
        return (0.0, 0.0, 0.0)
    if wavelength < 440: #violet/blue
        return ((440 - wavelength) / 60, 0.0, 1.0)
    if wavelength < 490: #blue/cyan
        return (0.0, (wavelength - 440) / 50, 1.0)
    if wavelength < 510: #cyan/green
        return (0.0, 1.0, (510 - wavelength) / 20)
    if wavelength < 580: #green/yellow
        return ((wavelength - 510) / 70, 1.0, 0.0)
    if wavelength < 645: #yellow/red
        return (1.0, (645 - wavelength) / 65, 0.0)
    return (1.0, 0.0, 0.0)

# Updates the firework color based on the selected elements generated using claude 
def update_firework(element):
    selected_colors = []
    for name, selected in zip(elements, check_buttons.get_status()):
        if selected:
            band = molecular_bands[element_bands[name]]
            selected_colors.append(wavelength_to_rgb(band))
    if selected_colors:
        firework_circle.set_facecolor(np.mean(selected_colors, axis=0))
    else:
        firework_circle.set_facecolor("black")
    fig.canvas.draw_idle()


check_buttons.on_clicked(update_firework)
axes[2, 2].set_visible(False)
update_firework(None)

plt.show()
"""
Author: Tayeb Kakeshpour
Description: This script generates figure 3
"""
import numpy as np
import matplotlib.pyplot as plt
from bin.tsitools import ParticleData
from scipy.optimize import curve_fit
from PIL import Image

def exponential_decay(t,lamda):
    return np.exp(-lamda * t) 

def fraction(T,lamda):
    alpha = 1.0/7
    return (alpha/lamda)*(1-np.exp(-T*lamda))

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
time_points = np.array([0,1,2,3,4,5])
fine_time_points = np.linspace(-0.5, 5.5, 100)
p1 = ParticleData('data',str(123))
p2 = ParticleData('data',str(124))
p3 = ParticleData('data',str(125))
bins = [
        ['bin1','0.3-2.7 $\mu$m',[1,11]],
        ['bin2','2.7-4.2 $\mu$m',[11,13]],
        ['bin3','4.2-6.5 $\mu$m',[13,15]],
        ['bin4','6.5-10 $\mu$m',[15,17]],
        ['bin5','$\geq$ 10 $\mu$m',[17,-1]],
        ]

# organize the data
for i in range(len(bins)):
    val1 = []
    val2 = []
    val3 = []

    for time in [60,120,180,240,300,360]:
        c1 = p1.get_particle_counts_at_time_np(time)[bins[i][2][0]-1:bins[i][2][1]-1]
        c2 = p2.get_particle_counts_at_time_np(time)[bins[i][2][0]-1:bins[i][2][1]-1]
        c3 = p3.get_particle_counts_at_time_np(time)[bins[i][2][0]-1:bins[i][2][1]-1]

        if i == 4:
            c1 = p1.get_particle_counts_at_time_np(time)[-1]
            c2 = p2.get_particle_counts_at_time_np(time)[-1]
            c3 = p3.get_particle_counts_at_time_np(time)[-1]            
        
        sum1 = np.sum(c1)
        sum2 = np.sum(c2)
        sum3 = np.sum(c3)

        if time == 60:
            constant1 = sum1
            constant2 = sum2
            constant3 = sum3
        
        val1.append(sum1)
        val2.append(sum2)
        val3.append(sum3)

    bins[i].append(val1)
    bins[i].append(val2)
    bins[i].append(val3)

fig, axs = plt.subplots(3, 1, figsize=(4.5, 3*3), dpi=300) 

rates = []
rates_e = []
##########################################################
### PANEL A ##############################################
##########################################################
tsi_time_points = np.linspace(0,6,20)
axs[0].errorbar(tsi_time_points , exponential_decay(t=tsi_time_points,lamda=1/7), fmt='--', capsize=5, color = 'k', label='loss to OPS')

for i in range(len(bins)):
    sum_values = np.sum(np.stack((bins[i][3], bins[i][4], bins[i][5]), axis=0), axis=0)
    sum_values = sum_values / sum_values[0]

    popt, pcov = curve_fit(exponential_decay, time_points, sum_values, p0=[0.4])
    rates.append(popt[0])

    predicted_values = exponential_decay(time_points, *popt)
    residuals = sum_values - predicted_values
    ss_res = np.sum(residuals**2)
    n = len(time_points)
    p = len(popt)
    standard_error = np.sqrt(ss_res / (n - p))
    rates_e.append(standard_error)

    fitted_curve = exponential_decay(fine_time_points, *popt)
    
    color = colors[i % len(colors)]
    norm_factor = exponential_decay(-0.5, *popt)

    axs[0].errorbar(time_points+0.5, sum_values/norm_factor, yerr=None, fmt='o', capsize=5, color=color, label=bins[i][1])
    axs[0].plot(fine_time_points+0.5, fitted_curve/norm_factor, color=color, linestyle='--')

axs[0].set_ylabel('$C_{\mathrm{N,D}}(t)/C_{\mathrm{N,D}}(0)$')
axs[0].set_xlabel('t (min)')
axs[0].set_xlim([0,6.2])
axs[0].legend(ncol=2, frameon=False, fontsize = 9)
##########################################################
### PANEL B ##############################################
##########################################################
bin = np.array([2.156,3.343,5.182,8.031])
bin_sq = bin*bin
rates = np.array(rates)
rates_fit = rates[:-1]
rates_e = np.array(rates_e)

d_loss, alpha = np.polyfit(bin_sq, rates_fit, 1)

bin_sq_t = np.linspace(0,120,100)
predicted_rates = np.poly1d((d_loss, alpha))(bin_sq_t)

axs[1].errorbar(bin_sq, rates[:-1], fmt='ok', yerr=rates_e[:-1], capsize = 3, label='rates obtained from panel a')
axs[1].errorbar(106.0, rates[-1], fmt='o',color= colors[4], yerr=rates_e[-1], capsize = 3,label= 'estimated from loss rate')
axs[1].plot(bin_sq_t, predicted_rates, '--k',label='fitted line')
axs[1].set_ylabel('λ (min⁻¹)')
axs[1].set_xlabel('D$^2$ ($\mu$m$^2$)')
axs[1].set_xlim([0,125])

##########################################################
### PANEL C ##############################################
##########################################################
diameter = np.array([0.300,0.374,0.465,0.579,0.721,0.897,1.117,1.391,\
                     1.732,2.156,2.685,3.343,4.162,5.182,6.451,8.031,10.300])
lamdas = np.poly1d((d_loss, alpha))(diameter*diameter)
measured_fractions = fraction(T=3,lamda=lamdas)

axs[2].plot(diameter,measured_fractions , '--k',label='fitted line')
axs[2].set_ylabel('$f$')
axs[2].set_xlabel('D ($\mu$m)')

##########################################################
### PANEL SETUP ##########################################
##########################################################
legend_font_size = 8.5
axs[0].text(5.9,0.97,'a',fontsize=13)
axs[0].legend(borderpad=0.5, labelspacing=0.5 , handletextpad=0.3,\
               ncols=2,frameon= False,fontsize=legend_font_size,
               bbox_to_anchor=(0.3, 1), loc='upper left', borderaxespad=0.)

axs[1].text(124,0.84,'b',fontsize=13)
axs[1].set_xlim([0,130])
axs[1].set_ylim([-0.1,0.93])
axs[2].text(10.3,0.315,'c',fontsize=13)
axs[2].set_yticks([0.2,0.3])

plt.tight_layout()
plt.savefig('figure3.pdf')
print('measured fractions = ', np.round(measured_fractions, 2).tolist())

fig.canvas.draw()
data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
image = Image.fromarray(data)
image = image.quantize(colors=256)
image.save('figure3.tiff', 'TIFF', compression='tiff_lzw')
plt.close()
exit()


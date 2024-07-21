"""
Author: Tayeb Kakeshpour
Description: This script generates figure 4
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from bin.tsitools import ParticleData
from bin.rhtools import calcRH,calcVol
from scipy.stats import gmean
import numpy.ma as ma
from PIL import Image

def addFlatten(python_list, np_array):
    python_list.extend(np_array.tolist())
    return python_list

def custom_mean(arr, axis=0):
    arithmetic_means = np.mean(arr, axis=axis)
    arr_nonzero = np.where(arr == 0, arithmetic_means, arr)
    geometric_means = gmean(arr_nonzero, axis=axis)
    logs = np.log(arr_nonzero)
    log_std = np.std(logs, axis=axis)
    multiplier = np.exp(log_std)
    lower_bounds = geometric_means / multiplier
    upper_bounds = geometric_means * multiplier
    masked_result = np.ma.masked_equal(geometric_means, 0)
    error_bars = np.array([geometric_means - lower_bounds, upper_bounds - geometric_means])
    return masked_result, error_bars

Chigh = 33.0 # mg/L at 83.4% RH and 35 C
Clow  = 18.54 # mg/L at 100% RH and 21.2 C
Vlow  = 7.0  # L - the volume of the cylinder
Vin   = 5.09 # input voltage 
Texp  = 21.0 # temperature
df = pd.read_csv('data.csv')

fig = plt.figure(figsize=(12, 12), dpi=300)
gs = gridspec.GridSpec(4, 3, width_ratios=[1, 1, 0.5])
axes = [[None for _ in range(3)] for _ in range(4)]
for row in range(4):
    for col in range(3):
        axes[row][col] = plt.subplot(gs[row, col])

measured_fractions =  np.array([0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.32, 0.31, 0.31, 0.3, 0.29, 0.27, 0.24, 0.21, 0.17])

Activities = ['breath', 'AAH','popeye','counting']
Samplings = ['direct','cold-short-funnel','warm-short-funnel','cold-long-funnel','warm-long-funnel']
colors = ['purple' ,'blue', 'red', 'green', 'orange','grey']
location = [0,-1.5,-0.5,+0.5,+1.5]
texts    = ["breath","'aah'", "'popeye'",'counting']
RHs = []
capsize = 2
alpha_value = 0.4

for i in range(len(Activities)):
    activity = Activities[i]
    RHs.append([Activities[i]])
    VLs = []
    for j in range(len(Samplings)):
        sampling = Samplings[j]
        dN_dlogD = []
        dV_dlogD = []
        RH_current = []
        VL_current = []
        for k in range(len(df)):
            if df.iloc[k]['Activity'] == Activities[i] and df.iloc[k]['Sampling'] == Samplings[j]:

                RH1 = calcRH(float(df.iloc[k]['RH1']),Vin,Texp)
                RH2 = calcRH(float(df.iloc[k]['RH2']),Vin,Texp)
                RH_change = RH2-RH1
                RH_current.append(RH_change)

                if Samplings[j] == 'direct':
                    volume =  calcVol(Chigh,Clow,Vlow,RH_change)
                else:
                    volume = 1.

                TSI = ParticleData('data',str(df.iloc[k]['TSI']))
                # multiply by 3 the function averages and sums are needed
                # divide by 0.43 for total volume correction
                dN_dlogD_current = TSI.get_dN_dlogD_np()/measured_fractions/volume*3/0.43
                dV_dlogD_current = TSI.get_dV_dlogD_np()/measured_fractions/volume*3/0.43

                dN_dlogD.append(dN_dlogD_current)
                dV_dlogD.append(dV_dlogD_current)
                
                VL_current.append(np.sum(dV_dlogD_current*TSI.get_dlogD_np()))

        RH_avg, RH_std = np.average(np.array(RH_current)),np.std(np.array(RH_current))
        RHs[i].append([Samplings[j].replace('-', ' '),RH_avg,RH_std])
        VL_gmean = gmean(np.array(VL_current))
        log_std = np.std(np.log(VL_current))
        lower_error = VL_gmean / np.exp(log_std)
        upper_error = VL_gmean * np.exp(log_std)        
        VL_errors = [[VL_gmean - lower_error], [upper_error - VL_gmean]]
        VLs.append([Samplings[j].replace('-', ' '),VL_gmean,VL_errors])   

        dN_dlogD = np.array(dN_dlogD)
        dV_dlogD = np.array(dV_dlogD)

        dN_dlogD_gmean, dN_dlogD_error = custom_mean(dN_dlogD)

        dV_dlogD_gmean, dV_dlogD_error  = custom_mean(dV_dlogD)

        meanDiameters, binWidths = TSI.get_mean_diameters_np(),  TSI.get_bin_widths_np()
        error_properties = {'ecolor': colors[j], 'capsize': capsize, 'alpha': alpha_value}

        axes[i][0].errorbar(meanDiameters+(location[j])*0.2*binWidths,dN_dlogD_gmean,yerr=dN_dlogD_error,color=colors[j],**error_properties)
        axes[i][1].errorbar(meanDiameters+(location[j])*0.2*binWidths,dV_dlogD_gmean,yerr=dV_dlogD_error,color=colors[j],**error_properties)

        axes[i][0].plot(meanDiameters+(location[j])*0.2*binWidths,dN_dlogD_gmean,color=colors[j],marker='o',markersize=3,label=Samplings[j].replace('-', ' '))
        axes[i][1].plot(meanDiameters+(location[j])*0.2*binWidths,dV_dlogD_gmean,color=colors[j],marker='o',markersize=3,label=Samplings[j].replace('-', ' '))

        if j in [0,1]:
            axes[i][j].set_xscale('log')
            axes[i][j].set_yscale('log')

        axes[i][2].set_yscale('log')

    for w, (subcategory, value, std) in enumerate(VLs):
        axes[i][2].bar(w, value, color=colors[w], yerr=std, capsize=5)
        axes[i][2].set_xticks([]) 
    
    axes[i][0].set_ylim(1,500000)
    axes[i][1].set_ylim(1,5000000)
    axes[i][2].set_ylim(10,1000000)
    axes[i][0].set_xlim(0.25,14)
    axes[i][1].set_xlim(0.25,14)
    axes[i][0].text(.28, 2, '$%s$' % texts[i], fontsize=11, ha='left')


#############################################################################################
### extract the data for the background (BS) and plot them on one of the panels #############
#############################################################################################
dN_dlogD_BG = []
dV_dlogD_BG = []
VL_BG = []
volume = 1.
for k in range(len(df)):
    if df.iloc[k]['Activity'] == 'hold':
        TSI = ParticleData('data',str(df.iloc[k]['TSI']))
        dN_dlogD_BG.append(TSI.get_dN_dlogD_np()/measured_fractions/volume*3/0.43)
        dV_dlogD_BG.append(TSI.get_dV_dlogD_np()/measured_fractions/volume*3/0.43)
        VL_BG.append(np.sum((TSI.get_dV_dlogD_np()/measured_fractions/volume*3/0.43)*TSI.get_dlogD_np()))

dN_dlogD_BG_gmean,dN_dlogD_BG_error = custom_mean(np.array(dN_dlogD_BG))
dV_dlogD_BG_gmean,dV_dlogD_BG_error = custom_mean(np.array(dV_dlogD_BG))


error_properties = {'ecolor': colors[-1], 'capsize': capsize, 'alpha': alpha_value}
axes[0][0].errorbar(meanDiameters+(0)*0.2*binWidths, dN_dlogD_BG_gmean, yerr=dN_dlogD_BG_error,**error_properties)
axes[0][0].plot(meanDiameters+(0)*0.2*binWidths, dN_dlogD_BG_gmean,'o-',markersize=3,color = colors[-1],label='background')
axes[0][1].errorbar(meanDiameters+(0)*0.2*binWidths, dV_dlogD_BG_gmean, yerr=dV_dlogD_BG_error,**error_properties)
axes[0][1].plot(meanDiameters+(0)*0.2*binWidths, dV_dlogD_BG_gmean,'o-',markersize=3,color = colors[-1],label='background')
#############################################################################################
### set the axis labels #####################################################################
#############################################################################################
[axes[q][0].set_ylabel('dN/dlogD (L$^{-1}$)',fontsize=11) for q in range(4)]
[axes[q][1].set_ylabel('dV/dlogD (fL$\cdot$L$^{-1}$)',fontsize=11) for q in range(4)]
[axes[q][2].set_ylabel('total (fL$\cdot$L$^{-1}$)',fontsize=11) for q in range(4)]
[axes[3][q].set_xlabel('D (\u03BCm)',fontsize=11) for q in range(2)]

labels = ['a','b','c','d','e','f','g','h','i','j','k','l']
[axes[q][0].text(11.5, 150000, labels[3*q], fontsize=16, fontfamily='Times New Roman') for q in range(4)]
[axes[q][1].text(11.5, 1000000, labels[(3*q)+1], fontsize=16, fontfamily='Times New Roman') for q in range(4)]
[axes[q][2].text(4.15, 350000, labels[(3*q)+2], fontsize=16, fontfamily='Times New Roman') for q in range(4)]

axes[0][0].legend(frameon=False,ncols=2, fontsize=9, loc='upper left')
plt.tight_layout()
plt.savefig('figure4.pdf')

fig.canvas.draw()
data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
image = Image.fromarray(data)
image = image.quantize(colors=256)
image.save('figure4.tiff', 'TIFF', compression='tiff_lzw')
plt.close()
exit()
 


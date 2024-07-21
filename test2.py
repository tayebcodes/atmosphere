"""
Author: Tayeb Kakeshpour
Description: This script generates figure 2
"""
from bin.rhtools import calcRH,calcVol
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from PIL import Image

Chigh = 33.0 # mg/L at 83.4% RH and 35 C
Clow  = 18.54 # mg/L at 100% RH and 21.2 C
Vlow  = 7.0  # L - the volume of the cylinder
Vin   = 5.09 # input voltage
Texp  = 21.0 # temperature

df = pd.read_csv('./data.csv')

Activities = ['breath', 'AAH','popeye','counting']
Samplings = ['direct','cold-short-funnel','warm-short-funnel','cold-long-funnel','warm-long-funnel']
colors = ['purple' ,'blue', 'red', 'green', 'orange','grey']
texts    = ["breath","'aah'", "'popeye'",'counting']
RHs = []
     
samplings_dict = {
    'direct': 'D',
    'cold-short-funnel': 'UFST',
    'warm-short-funnel': 'HFST',
    'cold-long-funnel': 'UFLT',
    'warm-long-funnel': 'HFLT'
}

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(3.5, 5), dpi=300)
###########################################################################
### PANEL A ###############################################################
###########################################################################
activity = 'breath'
samplings = ['cold-short-funnel','warm-short-funnel','cold-long-funnel','warm-long-funnel']
categories = []
values = []
errors = []
data = []

for i in range(len(samplings)):
    sampling = samplings[i]
    RH_current = []
    for j in range(len(df)):
        if df.iloc[j]['Activity'] == activity and df.iloc[j]['Sampling'] == sampling:
            RH1 = calcRH(float(df.iloc[j]['RH1']),Vin,Texp)
            RH2 = calcRH(float(df.iloc[j]['RH2']),Vin,Texp)
            RH_change = RH2-RH1
            RH_current.append(RH_change)
    categories.append(samplings_dict[sampling])
    values.append(np.mean(RH_current))
    errors.append(np.std(RH_current))

axes[0,0].bar(categories, values, yerr=errors, capsize=5, color='skyblue')
axes[0,0].set_ylabel('ΔRH (%)')
axes[0,0].set_title('a')
axes[0,0].xaxis.set_major_locator(FixedLocator(range(len(categories))))
axes[0,0].set_xticklabels(categories, rotation=45, ha='right')

###########################################################################
### PANEL B ###############################################################
###########################################################################
sampling = 'warm-short-funnel'
activities = ['breath', 'AAH','popeye','counting']
categories = []
values = []
errors = []
data = []

for i in range(len(activities)):
    activity = activities[i]
    RH_current = []
    for j in range(len(df)):
        if df.iloc[j]['Activity'] == activity and df.iloc[j]['Sampling'] == sampling:
            RH1 = calcRH(float(df.iloc[j]['RH1']),Vin,Texp)
            RH2 = calcRH(float(df.iloc[j]['RH2']),Vin,Texp)
            RH_change = RH2-RH1
            RH_current.append(RH_change)
    categories.append(activity)
    values.append(np.mean(RH_current))
    errors.append(np.std(RH_current))

axes[0,1].bar(categories, values, yerr=errors, capsize=5, color='skyblue')
axes[0,1].set_ylabel('ΔRH (%)')
axes[0,1].set_title('b')
axes[0,1].xaxis.set_major_locator(FixedLocator(range(len(categories))))
axes[0,1].set_xticklabels(categories, rotation=45, ha='right')

###########################################################################
### PANEL C ###############################################################
###########################################################################
sampling = 'direct'
activities = ['breath', 'AAH','popeye','counting']
categories = []
values = []
errors = []
data = []

for i in range(len(activities)):
    activity = activities[i]
    RH_current = []
    for j in range(len(df)):
        if df.iloc[j]['Activity'] == activity and df.iloc[j]['Sampling'] == sampling:
            RH1 = calcRH(float(df.iloc[j]['RH1']),Vin,Texp)
            RH2 = calcRH(float(df.iloc[j]['RH2']),Vin,Texp)
            RH_change = RH2-RH1
            RH_current.append(RH_change)
    categories.append(activity)
    values.append(np.mean(RH_current))
    errors.append(np.std(RH_current))

axes[1,0].bar(categories, values, yerr=errors, capsize=5, color='skyblue')
axes[1,0].set_ylabel('ΔRH (%)')
axes[1,0].set_title('c')
axes[1,0].xaxis.set_major_locator(FixedLocator(range(len(categories))))
axes[1,0].set_xticklabels(categories, rotation=45, ha='right')

###########################################################################
### PANEL D ###############################################################
###########################################################################
sampling = 'direct'
activities = ['breath', 'AAH','popeye','counting']
categories = []
values = []
errors = []
data = []

for i in range(len(activities)):
    activity = activities[i]
    RH_current = []
    for j in range(len(df)):
        if df.iloc[j]['Activity'] == activity and df.iloc[j]['Sampling'] == sampling:
            RH1 = calcRH(float(df.iloc[j]['RH1']),Vin,Texp)
            RH2 = calcRH(float(df.iloc[j]['RH2']),Vin,Texp)
            RH_change = RH2-RH1
            volume =  calcVol(Chigh,Clow,Vlow,RH_change)
            RH_current.append(volume)
    categories.append(activity)
    values.append(np.mean(RH_current))
    errors.append(np.std(RH_current))

axes[1,1].bar(categories, values, yerr=errors, capsize=5, color='peachpuff')
axes[1,1].set_ylabel('V$_{c}$ (L)')
axes[1,1].set_title('d')
axes[1,1].xaxis.set_major_locator(FixedLocator(range(len(categories))))
axes[1,1].set_xticklabels(categories, rotation=45, ha='right')


plt.tight_layout()
plt.savefig('figure2.pdf')

fig.canvas.draw()
data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
image = Image.fromarray(data)
image = image.quantize(colors=256)
image.save('figure2.tiff', 'TIFF', compression='tiff_lzw')

plt.close()
exit()

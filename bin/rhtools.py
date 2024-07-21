"""
Author: Tayeb Kakeshpour
Description: Tools for calculating humidity
"""

def calcRH(Vout,Vin,T):
	sensorRH = ((Vout/Vin)- 0.1515)/0.00636
	trueRH = sensorRH/(1.0546-0.00216*T)
	return trueRH

def calcVol(Chigh,Clow,Vlow,RHlow):
	Vhigh = (Clow * Vlow * (RHlow/100))/Chigh
	return Vhigh

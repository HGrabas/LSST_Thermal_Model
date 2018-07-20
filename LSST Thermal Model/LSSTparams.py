import math
#Parameters used for the SPICE circuit simulation

#Parameters used for the model
CCD_DC_current = 3.84
CCD_temp = 173
Cryo_temp = 145

R1V = 0.127
R3V = 0.192
R4V = 0.162
R5V = 0.1020
R6V = 90
R7V = 30.1
R15V = 199
R16V = 0.013
R18V = 181
R19V = 558

#PI1 params
KG1 = 15
KIG1 = 1/500
AWG1 = 4
MING1 = 0
MAXG1 = 2

#PI2 params
KG2 = 0.2
KIG2 = 1/1000
AWG2 = 4
MING2 = -8
MAXG2 = 0

#Lense parameters
Vair = 1
TA = 298
AT = 0.0144

EN = 0.65
EL = 0.91
EE = 0.91
ELL = 1/(1/EN+1/EL-1)
EAL = 1/(1/EE+1/EL-1)

Rback = 10000
Rfront = 21.9/math.sqrt(Vair)
IG = 5.67e-8

LensGainFront = '{:0.3e}'.format(abs(AT*IG*EAL))
LensGainBack  = '{:0.3e}'.format(abs(AT*IG*ELL))
SurfLimMin = 100
SurfLimMax = 400

Heater_gain= 4
Cryo_gain = -2

Cryo_table = (('0V','0'),('126.7V','15.2'),('135.7V','17.4'),('195.9V','19.7'),('220V','28.8'),('289V','71.7'))

AG = 0.022
EAF = 0.029
EAG = 0.03
ShroudGainFront = '{:0.3e}'.format(IG*AG*EAF)
ShroudGainBack  = '{:0.3e}'.format(IG*AG*EAG)

EB = 0.05
ER = 0.05
ERB = EB*ER/(EB+ER-EB*ER)
BB = 0.092
RebGain = '{:0.3e}'.format(IG*ERB*BB)
RebRes = 46

EP = 0.05
PP = 0.03
GridGain = '{:0.3e}'.format(IG*EP*PP)

#Cold circuit model
CurrSrccurrent = 17.28
ASPICcurrent   = 4.35
CABACcurrent   = 9
ADCcurrent     = 0.42
DiffAmpcurrent = 3.63
Heatercurrent  = 9.2
FPGAcurrent    = 9.45

#Cold plate PI regulator
KG3   = 10
KIG3  = 1/100
AWG3  = 4
MING3 = 0
MAXG3 = 20
ColdTemp = 251
Cold_transcd = 6
ColdPlateCap = 385

Cold_table = (('140','0'),('228','29'),('233','45'),('238','63'),('243','85'),('248','111'),('253','142'),('258','177'),('300','200'))

RShroud = 100000

GridR = 200



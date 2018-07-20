#!/usr/bin/env python

"""
The cryo circuit starts with the sensors and the heat sources
connected to the sensors and runs through the raft base plates,
the thermal straps, the REC box sides to the cryo plate,
where the heat is removed by the cryo refrigerator evaporator
loops. The schematic for this circuit is shown in Fig. 1.
Material is represented by the transmission lines, TLURC8.
TLURC stands for lumped RC transmission. The 8 indicates that
the overall capacitance and resistance has been divided into
8 pieces. Most of the resistors represent thermal-contact
impedances. R7 represents the flex cables (all of them in one
RTM) connecting the sensors to the REBs (all three of them in
an RTM).  PID and PID2 are the proportional-integral
controllers, which regulate the temperatures of the sensors
and the cryo plate respectively by controlling heaters G1
and G2. 

"""

__author__ = "Hervé Grabas"
__email__ = "herve.grabas@gmail.com"

import os
import csv
import numpy as np
import matplotlib.pyplot as plt

import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()

from PySpice.Probe.Plot import plot
from PySpice.Spice.Netlist import Circuit
from PySpice.Spice.Netlist import SubCircuitFactory
from PySpice.Spice.NgSpice.Shared import NgSpiceShared
from PySpice.Unit import *

import LSSTparams as par

os.environ['PySpiceLogLevel']="ERROR"


# Some proto code to use callbacks in PySpice - work in progress
# Would be usefull to do a software (Python) PID block
#####
class MyNgSpiceShared(NgSpiceShared):
  def __init__(self, gain, awgain, timeConstant, 
               minOutput, maxOutput, smoothTime,
               setPoint, tolerance, **kwargs):
      super().__init__(**kwargs)
      self.gain = gain
      self.awgain = awgain
      self.timeConstant = timeConstant
      self.minOutput = maxOutput
      self.smoothTime = smoothTime
      self.setPoint = setPoint
      self.tolerance = tolerance
      
      self.actualtime = 0
      
  #Java code for the PID block   
  #def compute_Pi(self,voltage, time):
  #  measTime = time
  #  if (lastTime == 0.) {
  #      lastTime = measTime - 1.;
  #      smoothInput = aveInput;
  #  }                            # initialize the filter
  #
  #  # perform the filter
  #  smoothInput = ((measTime - lastTime) * aveInput
  #      + (smoothTime - measTime + lastTime) * smoothInput) / smoothTime;
  #
  #  # find the "error"
  #  setError = (setpoint - smoothInput);
  #
  #  # apply the overall loop gain
  #  smoothOutput = gain * setError;
  #
  #  # add the proportional and integral terms
  #  propInt = smoothOutput + errorIntegral;
  #
  #  # Make sure the final result is within bounds, truncate otherwise
  #  output = (propInt > maxOutput) ? maxOutput :
  #           (propInt < minOutput) ? minOutput : propInt;
  #
  #  # calculate the amount of output truncation
  #  outputTrunc = awgain * (output - propInt);
  #
  #  # integrate the error less amount proportinal to the truncation
  #  errorIntegral +=
  #      (smoothOutput + outputTrunc) * (measTime - lastTime) / timeConst;
  #
  #  lastTime = measTime;
  #
  #  # add base output
  #  return output + baseOutput;
  
  def get_vsrc_data(self, voltage, time, node, ngspice_id):
    self._logger.debug('ngspice_id-{} get_vsrc_data @{} node {}'.format(ngspice_id, time, node))
    voltage[0] = self._amplitude * math.sin(self._pulsation * time)
    return 0

  def get_isrc_data(self, current, time, node, ngspice_id):
    self._logger.debug('ngspice_id-{} get_isrc_data @{} node {}'.format(ngspice_id, time, node))
    current[0] = 1.
    return 0

  def send_data(self, actual_vector_values, number_of_vectors, ngspice_id):
    for i in range(int(number_of_vectors)):
        actual_vector_value = data.vecsa[i]
        print("    Vector: {} {} +i {}".format(ffi.string(actual_vector_value.name),
                                               actual_vector_value.creal,
                                               actual_vector_value.cimag))
    
    return 0
  def get_sync_data(time, deltatime, olddelta, ngspice_id, id_number, user_data):
    id_number = self._ngspice_id
    

  def send_init_data_callback(data,  ngspice_id, user_data):
      print(">>> send_init_data")
      number_of_vectors = data.veccount
      for i in range(number_of_vectors):
          print("  Vector:", ffi.string(data.vecs[i].vecname))

      return 0
###### End of proto code

#PI Subcircuit
class PI(SubCircuitFactory):
  __name__ = 'PI'
  __nodes__ = ('In+', 'In-','Out', 'K=1','KI=1','AWG=1','MINV=1','MAXV=1')
  
  def __init__(self):
    super().__init__()
    self.BehavioralSource(1,'prop'    ,self.gnd,voltage_expression ='{K}*(V(In+)-V(In-))')
    self.BehavioralSource(2,'sum'     ,self.gnd,voltage_expression ='V(prop)+V(integ)')
    self.BehavioralSource(3,'overfl'  ,self.gnd,voltage_expression ='{AWG}*(V(Out)-V(sum))')
    self.BehavioralSource(4,'antiwind',self.gnd,voltage_expression ='V(overfl)+V(prop)')
    self.BehavioralSource(5,'Out'     ,self.gnd,voltage_expression ='V(sum) < {MINV} ? {MINV} : V(sum) > {MAXV} ? {MAXV} : V(sum)')
    
    #integrator part - fixme with int model from NGSpice?
    #circuit.raw_spice += 'R2 out 0 1kOhm'
    self.VoltageControlledCurrentSource(1,self.gnd,'integ','antiwind',self.gnd,transconductance='{KI}')
    self.C(1,'integ',self.gnd,1)
    self.R(1,'integ',self.gnd,10000000) 

#BoltzmannSurface subcircuit
class BoltzmannSurface(SubCircuitFactory):
  __name__ = 'BoltzmannSurface'
  __nodes__ = ('In','Out','Gain=1','Rin=1','VMIN=0','VMAX=1')
  
  def __init__(self):
    super().__init__()
    self.R(1,'In','Out','{Rin}')
    self.BehavioralSource(1,'V+' ,self.gnd,voltage_expression ='V(In) < {VMIN} ? {VMIN} : V(In) > {VMAX} ? {VMAX} : V(In)')
    self.BehavioralSource(2,'V-' ,self.gnd,voltage_expression ='V(Out) < {VMIN} ? {VMIN} : V(Out)> {VMAX} ? {VMAX} : V(Out)')
    self.BehavioralSource(3,'vp+',self.gnd,voltage_expression ='POW(V(V+),4)')
    self.BehavioralSource(4,'vp-',self.gnd,voltage_expression ='POW(V(V-),4)')
    self.VoltageControlledCurrentSource(1,'In','Out','vp+','vp-',transconductance='{Gain}')


#Definition of the LSST thermal model circuit - One RAFT
circuit = Circuit("MainCryo")

#PI regulators
circuit.subcircuit(PI())
circuit.X(1,'PI','CCDSetPoint'  ,'L3_CCD'   ,'HeaterCtrl'    ,K=par.KG1,KI=par.KIG1,AWG=par.AWG1,MINV=par.MING1,MAXV=par.MAXG1)
circuit.X(2,'PI','PlateSetPoint','CryoPlate','PlateCtrl'     ,K=par.KG2,KI=par.KIG2,AWG=par.AWG2,MINV=par.MING2,MAXV=par.MAXG2)
circuit.X(3,'PI','ColdSetPoint' ,'ColdBAR'  ,'ColdHeaterCtrl',K=par.KG3,KI=par.KIG3,AWG=par.AWG3,MINV=par.MING3,MAXV=par.MAXG3)

#Heater functions
circuit.VoltageControlledCurrentSource(1,circuit.gnd,'HeaterOut','HeaterCtrl'    ,circuit.gnd,transconductance=par.Heater_gain)
circuit.VoltageControlledCurrentSource(2,circuit.gnd,'PlateIn'  ,'PlateCtrl'     ,circuit.gnd,transconductance=par.Cryo_gain)
circuit.VoltageControlledCurrentSource(6,circuit.gnd,'ColdPlate','ColdHeaterCtrl',circuit.gnd,transconductance=par.Cold_transcd)

#Behavioral models
circuit.BehavioralSource(1,circuit.gnd,'Heatercurrent',current_expression='4*SQRT(MAX(V(HeaterCtrl)*(2.01-V(HeaterCtrl)),0))')

circuit.NonLinearCurrentSource(3,'PlateIn'  ,circuit.gnd,expression='V(PlateIn, 0)'  ,table=par.Cryo_table)
circuit.NonLinearCurrentSource(7,'ColdPlate',circuit.gnd,expression='V(ColdPlate, 0)',table=par.Cold_table)

#Emissive surface models
circuit.subcircuit(BoltzmannSurface())
circuit.X(4,'BoltzmannSurface','L3_CCD'    ,'LensBack'      ,Gain=par.LensGainBack   ,Rin=par.Rback  ,VMIN=par.SurfLimMin,VMAX=par.SurfLimMax)
circuit.X(5,'BoltzmannSurface','Ambient'   ,'LensFront'     ,Gain=par.LensGainFront  ,Rin=par.Rfront ,VMIN=par.SurfLimMin,VMAX=par.SurfLimMax)
circuit.X(6,'BoltzmannSurface','Ambient'   ,'ShroudFront'   ,Gain=par.ShroudGainFront,Rin=par.RShroud,VMIN=par.SurfLimMin,VMAX=par.SurfLimMax)
circuit.X(7,'BoltzmannSurface','ShroudBack','CryoGridShroud',Gain=par.ShroudGainBack ,Rin=par.RShroud,VMIN=par.SurfLimMin,VMAX=par.SurfLimMax)
circuit.X(8,'BoltzmannSurface','REBRAD'    ,'BoxRAD'        ,Gain=par.RebGain        ,Rin=par.RebRes ,VMIN=par.SurfLimMin,VMAX=par.SurfLimMax)
circuit.X(9,'BoltzmannSurface','ColdPlate' ,'CryoPlate'     ,Gain=par.GridGain       ,Rin=par.GridR  ,VMIN=par.SurfLimMin,VMAX=par.SurfLimMax)

#Lumped resistivity elements
circuit.UniformDistributedRCLine('T1' ,'CryoPlate'     ,'PlateIn'       ,circuit.gnd,model='T1' ,length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T2' ,'T28out'        ,'CurrentSource' ,circuit.gnd,model='T2' ,length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T3' ,'BoxIn'         ,'BoxRAD'        ,circuit.gnd,model='T3' ,length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T4' ,'StrapIn'       ,'StrapOut'      ,circuit.gnd,model='T4' ,length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T5' ,'Flex_lumped'   ,'HeaterOut'     ,circuit.gnd,model='T5' ,length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T6' ,'L3_CCD'        ,'Flex_entry'    ,circuit.gnd,model='T6' ,length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T7' ,'LensBack'      ,'LensFront'     ,circuit.gnd,model='T7' ,length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T8' ,'CurrentSource' ,'T29out'        ,circuit.gnd,model='T8' ,length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T9' ,'T29out'        ,'CurrentASPIC'  ,circuit.gnd,model='T9' ,length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T10','CurrentASPIC'  ,'REBRAD'        ,circuit.gnd,model='T10',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T11','REBRAD'        ,'CABACcurrent'  ,circuit.gnd,model='T11',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T12','CABACcurrent'  ,'PassCryoIn'    ,circuit.gnd,model='T12',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T13','PassCryoIn'    ,'ColdBAR'       ,circuit.gnd,model='T13',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T14','T29out'        ,'ColdBAR'       ,circuit.gnd,model='T14',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T15','ColdBAR'       ,'ADCcurrent'    ,circuit.gnd,model='T15',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T16','FPGAIn'        ,'ColdBAR'       ,circuit.gnd,model='T16',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T17','ADCcurrent'    ,'DiffAmpcurrent',circuit.gnd,model='T17',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T18','ColdStrapIn'   ,'ColdPlate'     ,circuit.gnd,model='T18',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T20','DiffAmpcurrent','Heatercurrent' ,circuit.gnd,model='T20',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T21','Heatercurrent' ,'FPGAIn'        ,circuit.gnd,model='T21',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T22','FPGAIn'        ,'FPGAcurrent'   ,circuit.gnd,model='T22',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T23','FPGAcurrent'   ,'REBDAQ'        ,circuit.gnd,model='T23',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T24','BoxRAD'        ,'BoxOut'        ,circuit.gnd,model='T24',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T25','ShroudBack'    ,'ShroudFront'   ,circuit.gnd,model='T25',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T26','CryoGridShroud','CryoPlate'     ,circuit.gnd,model='T26',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T27','GridFront'     ,'GridBack'      ,circuit.gnd,model='T27',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T28','REBCCD'        ,'T28out'        ,circuit.gnd,model='T28',length=1,number_of_lumps=8)
circuit.UniformDistributedRCLine('T29','T28out'        ,'T29out'        ,circuit.gnd,model='T29',length=1,number_of_lumps=8)

#Heat generators
circuit.I(1,circuit.gnd,'CurrentSource' ,par.CurrSrccurrent)
circuit.I(2,circuit.gnd,'CurrentASPIC'  ,par.ASPICcurrent)
circuit.I(3,circuit.gnd,'CABACcurrent'  ,par.CABACcurrent)
circuit.I(4,circuit.gnd,'ADCcurrent'    ,par.ADCcurrent)
circuit.I(5,circuit.gnd,'DiffAmpcurrent',par.DiffAmpcurrent)
circuit.I(6,circuit.gnd,'Heatercurrent' ,par.Heatercurrent)
circuit.I(7,circuit.gnd,'FPGAcurrent'   ,par.FPGAcurrent)
circuit.I(8,circuit.gnd,'L3_CCD'        ,par.CCD_DC_current)

#Temperature sources and set points
circuit.V(1,'Ambient'      ,circuit.gnd,par.TA)
circuit.V(2,'PlateSetPoint',circuit.gnd,par.Cryo_temp)
circuit.V(3,'ColdSetPoint' ,circuit.gnd,par.ColdTemp)
circuit.V(4,'CCDSetPoint'  ,circuit.gnd,par.CCD_temp)

#Contact resistors
circuit.R(1 ,'BoxOut'    ,'CryoPlate'  ,'{F*'+str(par.R1V)+'}')
circuit.R(3 ,'StrapOut'  ,'BoxIn'      ,'{F*'+str(par.R3V)+'}')
circuit.R(4 ,'HeaterOut' ,'StrapIn'    ,'{F*'+str(par.R4V)+'}')
circuit.R(5 ,'Flex_entry','Flex_lumped','{FC*'+str(par.R5V)+'}')
circuit.R(6 ,'REBDAQ'    ,'Ambient'    ,par.R6V)
circuit.R(7 ,'Flex_entry','REBCCD'     ,par.R7V)
circuit.R(15,'Ambient'   ,'GridFront'  ,par.R15V)
circuit.R(16,'ColdBAR'   ,'ColdStrapIn',par.R16V)
circuit.R(18,'HeaterOut' ,'REBCCD'     ,par.R18V)
circuit.R(19,'Ambient'   ,'CryoPlate'  ,par.R19V)

#Capacitor
circuit.C(1,'ColdPlate',circuit.gnd,par.ColdPlateCap)

circuit.G1.output_minus.add_current_probe(circuit)
circuit.G2.output_minus.add_current_probe(circuit)
circuit.G6.output_minus.add_current_probe(circuit)

circuit.parameter('F=1','F=2') ##TODO report bad implementation 
circuit.parameter('FC={SQRT(F)}','F={SQRT(F)}') ##TODO report bad implementation 

#Lumped models created from .csv file - good practice to avoids errors
with open('LumpedValues.csv') as csvfile:
  reader = csv.DictReader(csvfile, delimiter='\t')
  for row in reader:
    circuit.model(row['Reference'],'URC',RPERL=row['R(Ohms)'],CPERL=row['C(F)'])
#End of circuit definition


#Definition of simulator
simulator = circuit.simulator(temperature=25, nominal_temperature=25)
#simulator.options(GMINSTEPS=10,GMIN=1e-18,PIVTOL=1e-15)

#Defining two simulation
analysis = simulator.dc(V1=slice(268,300,1))

#analysis2 = simulator.dc(F=slice(0.1,10,1))
analysis2 = simulator.operating_point()

simulator = circuit.simulator(temperature=25, nominal_temperature=25,F=2)

#Saving the SPICE netlist into a file
f1 = open('./maincryo.sp', 'w')
print(str(circuit), file=f1)
f1.close

print("Branches names:")
print(analysis.branches.keys())

#Plotting thr simualtion result
figure = plt.figure(1, (20, 10))
axe = plt.subplot(111)
plt.grid()
axe.plot(analysis.ambient, -analysis.vg1_output_minus, 'r-+')
axe.plot(analysis.ambient, -analysis.vg2_output_minus, 'b-+')
axe.set_xlabel('Ambient temperature [T K]')
axe.set_ylabel('Heat [W]')
plt.legend(('Raft heat', 'Cryo plate heat'))
plt.tight_layout()
plt.show()

#Printing the DC node voltages = temperatures in K
for node in analysis2.nodes.values():
    print('Node {}: {:4.1f} V'.format(str(node), float(node)))

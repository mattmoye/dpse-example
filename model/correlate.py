#####################################################
#
#  7 January 2010
#  Bryan A. Toth
#  University of California, San Diego
#  btoth@physics.ucsd.edu
#
#  This script correlates couplings with their
#  appropriate equation for use in determining the R-
#  calue of the optimization.
#
#  This script has been developed as part of a suite of
#  python scripts to define a dynamic parameter estimation
#  problem using the optimization software IPOPT, but is
#  generally applicable to any application needing
#  discretized derivatives of a vector field.
#
######################################################

import sympy as sym

import discretize

nY = discretize.nY
nP = discretize.nP
nU = discretize.nU
nI = discretize.nI

Lvars = discretize.Lvars
Lparams = discretize.Lparams
Ldata = discretize.Ldata
cup = discretize.Lcouple
Lstimuli = discretize.Lstimuli

#Make symbols using sympy module
Sv = []
Sp = []
Sk = []
Sd = []
Si = []
for i in range(len(Lvars)):
  Sv.append(sym.Symbol(Lvars[i]))
for i in range(len(Lparams)):
  Sp.append(sym.Symbol(Lparams[i]))
for i in range(nU):
  Sd.append(sym.Symbol(Ldata[i]))
  Sk.append(sym.Symbol(cup[2*i]))
  Sk.append(sym.Symbol(cup[2*i+1]))
  # Sk includes the coupling and the derivative of the coupling: k1,k1d,k2,k2d,etc ...
for i in range(nI):
  Si.append(sym.Symbol(Lstimuli[i]))


eqns = discretize.Feqnstr

correl = []
for i in range(nY):
  correl.append(0)
func = []
synch = []
for j in range(nU):
   for i in range(nY):
     foo = eqns[i].find(cup[2*j])
     if foo == -1:
       continue
     else:
        correl[i] = 1
        temp = eqns[i].partition(cup[2*j])
        temp1 = temp[0]
        func.append(temp1.rpartition('+')[0])
        synch.append(temp[1]+temp[2])
#Feqns = []
#for k in range(nU):
#  sTemp1 = func[k]
#  for i in range(len(Lvars)):
#    sTemp2 = "Sv[%d]" % i
#    sTemp1 = sTemp1.replace(Lvars[i],sTemp2)
#  for i in range(len(Lparams)):
#    sTemp2 = "Sp[%d]" % i
#    sTemp1 = sTemp1.replace(Lparams[i],sTemp2)
#  for i in range(len(cup)):
#   sTemp2 = "Sk[%d]" % i
#    sTemp1 = sTemp1.replace(cup[i],sTemp2)
#  for i in range(nU):
#    sTemp2 = "Sd[%d]" % i
#    sTemp1 = sTemp1.replace(Ldata[i],sTemp2)
#  for i in range(nI):
#    sTemp2 = "Si[%d]" % i
#    sTemp1 = sTemp1.replace(Lstimuli[i],sTemp2)
#  sTemp2 = "Feqns.append("
#  sTemp2 = sTemp2 + sTemp1 + ")"
#  exec sTemp2
#print func
#print Feqns

def subvars(mystr):
    mytemp = mystr
    mytemp = sym.sympify(mytemp)
    mytemp = sym.ccode(mytemp)

    for j in range(len(Sv)):
      Srep = "Xval[%d]" % j
      Sfind = Lvars[j]
      mytemp = mytemp.replace(Sfind,Srep)

    for j in range(len(Sp)):
      Srep = "Pval[%d]" % j
      Sfind = Lparams[j]
      mytemp = mytemp.replace(Sfind,Srep)

    for j in range(len(Sk)):
      Srep = "K11val[%d]" % (j/2)
      Sfind = cup[j]
      mytemp = mytemp.replace(Sfind,Srep)

    for j in range(len(Sd)):
      Srep = "Xdval[%d]" % j
      Sfind = Ldata[j]
      mytemp = mytemp.replace(Sfind,Srep)

    for j in range(len(Si)):
      Srep = "Ival[%d]" % j
      Sfind = Lstimuli[j]
      mytemp = mytemp.replace(Sfind,Srep)

    return mytemp

strEqns = []
for i in range(len(func)):
  Stemp = str(func[i])
  Stemp = subvars(Stemp)
  strEqns.append(Stemp)

cupEqns = []
for i in range(len(synch)):
  Stemp = str(synch[i])
  Stemp = subvars(Stemp)
  cupEqns.append(Stemp)

#print strEqns
#print cupEqns

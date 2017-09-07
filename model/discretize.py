#####################################################
#
#  20 October 2009
#  Bryan A. Toth
#  University of California, San Diego
#  btoth@physics.ucsd.edu
#
#  This script performs symbolic Hermite-Simpson
#  integration on a vector field given in the text file
#  equations.txt, takes the Jacobian and Hessian of the
#  vector field, and stores the results in arrays that
#  are used by other python scripts in this directory.
#
#  This script has been developed as part of a suite of
#  python scripts to define a dynamic parameter estimation
#  problem using the optimization software IPOPT, but is
#  generally applicable to any application needing
#  discretized derivatives of a vector field.
#
######################################################

import sympy as sym
from sympy import *
import re

#  Opening and reading the text file with the vector field information

file = open('equations.txt','r')
temp=[]  # Array to hold equations.txt information
for line in file:
  if line.startswith('#'): # Pound used as comment in text file
     continue
  elif line.startswith('\\'): # In case file has UTF-8 markers
     continue
  else:
     temp.append(line)
file.close()

h=[]  # Array to hold unformatted equations.txt information
for i in range(len(temp)):
  temp1=temp[i].rstrip( )
  h.append(temp1)

# Initialize problem variables
nY=0
nP=0
nU=0
nI=0
nF=0

# Problem name
Problem = h[0]

# Problem variables
a=h[1].split(',')
nY=int(a[0])
nP=int(a[1])
nU=int(a[2])
nI=int(a[3])
if len(a) > 4:
  nF=int(a[4])

# Import equations as strings
Feqnstr = []
for k in range(nY):
   Feqnstr.append(h[k+2])

# Import objective function as string
Fobjstr = []
Fobjstr.append(h[nY+2])

# Import variable, parameter, control, data, and stimuli names as strings
Lvars = []
for k in range(nY):
   Lvars.append(h[k+3+nY])

Lparams = []
for k in range(nP):
   Lparams.append(h[k+3+nY+nY])

Lcouple = []
Ldata = []
for k in range(nU):
    Lcouple.append(h[k+3+nY+nY+nP])
    Lcouple.append('d'+ h[k+3+nY+nY+nP])
    Ldata.append(h[k+3+nY+nY+nP+nU])
Lstimuli = []
for k in range(nI):
    Lstimuli.append(h[k+3+2*nY+nP+2*nU])

# Import function names as strings
Funcstr = []
Funcarg = []
for k in range(nF):
    temp = h[k+3+2*nY+nP+2*nU+nI].split(',')
    Funcstr.append(temp[0])
    Funcarg.append(int(temp[1]))
#Lvars.reverse()
#Lcouple.reverse()
Fdim = len(Feqnstr)
Pdim = len(Lparams)

print ("Making symbols using sympy module:")
Sv = []
Sp = []
Sk = []
Sd = []
Si = []
print ("  variables...")
for i in range(len(Lvars)):
  Sv.append(sym.Symbol(Lvars[i]))
print ("  parameters...")
for i in range(len(Lparams)):
  Sp.append(sym.Symbol(Lparams[i]))
print ("  data...")
for i in range(nU):
  Sd.append(sym.Symbol(Ldata[i]))
  Sk.append(sym.Symbol(Lcouple[2*i]))
  Sk.append(sym.Symbol(Lcouple[2*i+1]))
# Sk includes the coupling and the derivative of the coupling: k1,k1d,k2,k2d,etc ...
for i in range(nI):
  Si.append(sym.Symbol(Lstimuli[i]))

Sall = Sv + Sk + Sp

print ("  functions...")
Sf = []
for i in range(nF):
  Sf.append(sym.Function(Funcstr[i]))

hstep = sym.Symbol("hstep")
# Define symbolic vector field
Feqns = []
for k in range(Fdim):
  sTemp1 = Feqnstr[k]
  for i in range(len(Lvars)):
    sTemp2 = "Sv[%d]" % i
    sTemp1 = sTemp1.replace(Lvars[i],sTemp2)
  for i in range(len(Lparams)):
    sTemp2 = "Sp[%d]" % i
    sTemp1 = sTemp1.replace(Lparams[i],sTemp2)
  for i in range(len(Lcouple)):
    sTemp2 = "Sk[%d]" % i
    sTemp1 = sTemp1.replace(Lcouple[i],sTemp2)
  for i in range(nU):
    sTemp2 = "Sd[%d]" % i
    sTemp1 = sTemp1.replace(Ldata[i],sTemp2)
  for i in range(nI):
    sTemp2 = "Si[%d]" % i
    sTemp1 = sTemp1.replace(Lstimuli[i],sTemp2)
  for i in range(nF):
    sTemp2 = "Sf[%d]" % i
    sTemp1 = sTemp1.replace(Funcstr[i],sTemp2)
  sTemp2 = "Feqns.append("
  sTemp2 = sTemp2 + sTemp1 + ")"
  exec (sTemp2)

# Define symbolic objective function
Fobj = []
sTemp1 = Fobjstr[0]
for i in range(len(Lvars)):
  sTemp2 = "Sv[%d]" % i
  sTemp1 = sTemp1.replace(Lvars[i],sTemp2)
for i in range(len(Lparams)):
  sTemp2 = "Sp[%d]" % i
  sTemp1 = sTemp1.replace(Lparams[i],sTemp2)
for i in range(len(Lcouple)):
  sTemp2 = "Sk[%d]" % i
  sTemp1 = sTemp1.replace(Lcouple[i],sTemp2)
for i in range(nU):
  sTemp2 = "Sd[%d]" % i
  sTemp1 = sTemp1.replace(Ldata[i],sTemp2)
sTemp2 = "Fobj.append("
sTemp2 = sTemp2 + sTemp1 + ")"
exec (sTemp2)


# For a continuous version of the Jacobian and Hessian of the
#  vector field, the following can be printed.  Otherwise, these
#  are not needed in the script.

#--------Jacobian---------------
#J = diff(Feqns,Sall[:]) # Jacobian
#---------Hessian---------------
#H = diff(J,Sall[:]) # Hessian


# Perform Hermite-Simpson discretization and label these as constraints
#  in the optimization problem.

#  The format of these constraints is one element in AllCon for each
#   dynamical variable in the vector field, with each entry containing
#   three entries for the discretization, corresponding to the current
#   time-step, next time-step, and midpoint value, as defined by the
#   integration rule.  For other integration rules, this part must change.

AllCon = []
# Simpson Constraints
for k in range(len(Sv)):
  tCon = []
  tCon.append(Sv[k] + (hstep/6.0)*Feqns[k])
  tCon.append(-Sv[k] + (hstep/6.0)*Feqns[k])
  tCon.append((2.0*hstep/3.0)*Feqns[k])
  AllCon.append(tCon)

# Hermite Constraints
for k in range(len(Sv)):
  tCon = []
  tCon.append(0.5*Sv[k] + (hstep/8.0)*Feqns[k])
  tCon.append(0.5*Sv[k] - (hstep/8.0)*Feqns[k])
  tCon.append(-Sv[k])
  AllCon.append(tCon)

# Hermite Control Constraint
for k in range(nU):
  tCon = []
  tCon.append(0.5*Sk[2*k] + (hstep/8.0)*Sk[2*k+1])
  tCon.append(0.5*Sk[2*k] - (hstep/8.0)*Sk[2*k+1])
  tCon.append(-Sk[2*k])
  AllCon.append(tCon)

# Add objective function
AllObj = []
AllObj.append(Fobj[0])
AllObj.append(0)
AllObj.append(Fobj[0])

#  The following dictionaries are for a subsequent substitution in the
#   function subvars.  These are the actual variable array names that will
#   be used in the eventual c++ IPOPT program.


dict1 = {0:"Xval",1:"Xvalp1",2:"Xval2"}
dict2 = {0:"K11val",1:"K11valp1",2:"K11val2"}
dict3 = {0:"dK11val",1:"dK11valp1",2:""}
dict4 = {0:"Xdval",1:"Xdvalp1",2:"Xdval2"}
dict5 = {0:"Ival",1:"Ivalp1",2:"Ival2"}

#  The following function performs a variable substitution, and
#   is called later in the code.

def subvars(mystr,myi):
   mytemp = mystr
#  The following two lines are very important
#  Sympy converts all inputs into a simplified form
#  for calculations.  Specifically, this involves
#  any exponentials (a^b) put in the form a**b.
#  Whereas this form is acceptable for fortran outputs,
#  this needs to change to pow(a,b) for C++ outputs.
#  Sympify and ccode combine to make this transformation.
#  This transformation is done at this point in the code,
#  since sympify will not operate on an
#  expression that includes brackets - which are added
#  in this function.
   mytemp = sym.sympify(mytemp)
   mytemp = sym.ccode(mytemp)

   for j in range(len(Sv)):
#      Srep = dict1[n] + "[%d]" % (len(Sv)-j-1)
      Srep = dict1[n] + "[%d]" % j
      Sfind = Lvars[j]
      mytemp = mytemp.replace(Sfind,Srep)

   for j in range(len(Sp)):
      Srep = "Pval[%d]" % j
      Sfind = Lparams[j]
      mytemp = mytemp.replace(Sfind,Srep)

#   for j in range(len(Sk)):
#      Srep = dict2[n] + "[%d]" % (j/2)
#      Sfind = Lcouple[j]
#      mytemp = mytemp.replace(Sfind,Srep)
   for j in range(len(Sk)):
      if (j % 2 == 0):
#        Srep = dict2[n] + "[%d]" % ((len(Sk)-j-1)/2)
        Srep = dict2[n] + "[%d]" % (j/2)
        Sfind = Lcouple[j]
        mytemp = mytemp.replace(Sfind,Srep)

      elif (j % 2 == 1):
#        Srep = dict3[n] + "[%d]" % ((len(Sk)-j-1)/2)
        Srep = dict3[n] + "[%d]" % (j/2)
        Sfind = Lcouple[j]
        mytemp = mytemp.replace(Sfind,Srep)

   for j in range(len(Sd)):
      Srep = dict4[n] + "[%d]" % j
      Sfind = Ldata[j]
      mytemp = mytemp.replace(Sfind,Srep)

   for j in range(len(Si)):
      Srep = dict5[n] + "[%d]" % j
      Sfind = Lstimuli[j]
      mytemp = mytemp.replace(Sfind,Srep)

   return mytemp
# END subvars

def subfunc(mystr,myi):
   mytemp = mystr
   Dsearch = re.findall('D\(', mytemp)
   for num in range(len(Dsearch)):
     jacsearch = re.search('D\(([a-z]+)\((-?[0-9]+(\.[0-9]+)?, |[A-Za-z0-9]+\[[0-9]+\], )+(-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\), (-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\)', mytemp)
     hessearch = re.search('D\(([a-z]+)\((-?[0-9]+(\.[0-9]+)?, |[A-Za-z0-9]+\[[0-9]+\], )+(-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\), (-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\]), (-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\)', mytemp)
     if jacsearch:
       dervar = jacsearch.group(6)
       var = re.findall('[A-za-z0-9]+\[[0-9]+\]|-?[0-9]*\.?[0-9]+', jacsearch.group(0))
       for i in range(len(var)-1):
         if var[i] == dervar:
           jacnum = i+1
       rep = jacsearch.group(1) + 'jac('
       for i in range(len(var)-1):
         temp = var[i] + ','
         rep += temp
       rep += str(jacnum) +')'
       mytemp = re.sub('D\(([a-z]+)\((-?[0-9]+(\.[0-9]+)?, |[A-Za-z0-9]+\[[0-9]+\], )+(-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\), (-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\)', rep,mytemp,1)
     if hessearch:
       dervar1 = hessearch.group(6)
       dervar2 = hessearch.group(8)
       hesnum1 = 0
       hesnum2 = 0
       var = re.findall('[A-za-z0-9]+\[[0-9]+\]|-?[0-9]*\.?[0-9]+', hessearch.group(0))
       for i in range(len(var)-2):
         if var[i] == dervar1:
           hesnum1 = i+1
         if var[i] == dervar2:
           hesnum2 = i+1
       rep = hessearch.group(1) + 'hes('
       for i in range(len(var)-2):
         temp = var[i] + ','
         rep += temp
       rep += str(hesnum1) +','+ str(hesnum2)+ ')'
       mytemp = re.sub('D\(([a-z]+)\((-?[0-9]+(\.[0-9]+)?, |[A-Za-z0-9]+\[[0-9]+\], )+(-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\), (-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\]), (-?[0-9]+(\.[0-9]+)?|[A-Za-z0-9]+\[[0-9]+\])\)', rep, mytemp, 1)
   return mytemp
#end subfunc

print ("Building constraint equation strings...")
strAllCon = []
for icon in range(len(AllCon)):
   temp1 = []
   for n in [0,1,2]:
      Stemp = str(AllCon[icon][n])
      Stemp = subvars(Stemp,n)
      Stemp = subfunc(Stemp,n)
      temp1.append(Stemp)
   strAllCon.append(temp1)

print ("Building objective function strings...")
strObj = []
temp1 = []
for n in [0,1,2]:
    Stemp = str(AllObj[n])
    Stemp = subvars(Stemp,n)
    Stemp = subfunc(Stemp,n)
    temp1.append(Stemp)
strObj.append(temp1)

print ("Building Jacobian...")
sJac = []
for icon in range(len(AllCon)):
   temp1 = []
   for jvar in range(len(Sall)):
      temp2 = []
      for n in [0,1,2]:
         import pdb; pdb.set_trace() ## DEBUG ##
         Stemp = str(sym.diff(AllCon[icon][n],Sall[jvar]))
         Stemp = subvars(Stemp,n)
         Stemp = subfunc(Stemp,n)
         temp2.append(Stemp)
      temp1.append(temp2)
   sJac.append(temp1)

print ("Building objective gradient strings...")
sObj = []
for jvar in range(len(Sall)):
   temp2 = []
   for n in [0,1,2]:
      Stemp = str(sym.diff(AllObj[n],Sall[jvar]))
      Stemp = subvars(Stemp,n)
      Stemp = subfunc(Stemp,n)
      temp2.append(Stemp)
   sObj.append(temp2)

print ("Building Hessian...")
sHes = []
for icon in range(len(AllCon)):
   temp1 = []
   for jvar in range(len(Sall)):
      temp2 = []
      for kvar in range(len(Sall)):
         temp3 = []
         for n in [0,1,2]:
            Stemp = str(sym.diff(sym.diff(AllCon[icon][n],Sall[jvar]),Sall[kvar]))
            Stemp = subvars(Stemp,n)
            Stemp = subfunc(Stemp,n)
            temp3.append(Stemp)
         temp2.append(temp3)
      temp1.append(temp2)
   sHes.append(temp1)

print ("Adding Hessian for the objective function...")
temp1=[]
for jvar in range(len(Sall)):
    temp2 = []
    for kvar in range(len(Sall)):
        temp3 = []
        for n in [0,1,2]:
            Stemp = str(sym.diff(sym.diff(AllObj[n],Sall[jvar]),Sall[kvar]))
            Stemp = subvars(Stemp,n)
            Stemp = subfunc(Stemp,n)
            temp3.append(Stemp)
        temp2.append(temp3)
    temp1.append(temp2)
sHes.append(temp1)

# Fill out Jacobian vector
# Jacobian includes all constraints, but not objective
# VJac will include all non-zero elements of the Jacobian
#  with the format: (value, row, column, time) where time
#  refers to whether the element is current time, next time
#  or midpoint time in discretized form

VJac = []
for i in range(len(AllCon)):
   for j in range(len(Sall)):
         if j < (len(Sv) + len(Sk)): # Same as nY + 2*nU
             for n in [0,1,2]:
                F = sJac[i][j][n]
                if F != '0':
                    temp1 = []
                    temp1.append(F)
                    temp1.append(i)
                    temp1.append(j)
                    temp1.append(n)
                    VJac.append(temp1)
# Distinction is made between variable/control entries and parameter entries based
#  on how the sJac matrix is set up: derivatives with respect to the parameters
#  must take into account all time information - adding the current time, next
#  time, and midpoint time derivatives together into a single entry.

         else:
             F = sJac[i][j][0]
             if F != '0':
                temp2 = []
                temp1 = ''
                for n in [0,1,2]:
                    F = sJac[i][j][n]
                    if n == 0:
                        temp1 = temp1 + F
                    else:
                        temp1 = temp1 + '+' + F
                temp2.append(temp1)
                temp2.append(i)
                temp2.append(j)
                VJac.append(temp2)

# Fill out objective gradient
VObj = []
for j in range(len(Sall)):
    if j < (len(Sv) + len(Sk)):
        for n in [0,1,2]:
            F = sObj[j][n]
            if F != '0':
                temp1 = []
                temp1.append(F)
                temp1.append(j)
                temp1.append(n)
                VObj.append(temp1)
    else:
        F = sObj[j][0]
        if F != '0':
            temp2 = []
            temp1 = ''
            for n in [0,1,2]:
                F = sObj[j][n]
                if n == 0:
                    temp1 = temp1 + F
                else:
                    temp1 = temp1 + '+' + F
            temp2.append(temp1)
            temp2.append(j)
            VObj.append(temp2)

# Fill out Hessian vector
# This vector has the format [index, constraint, row, column, time, value] for each entry
# Index is a counter that refers to a specific row/column combination, so that the same
#  combination is not used more than once.
# Constraint refers to which constraint that the Hessian element is taking derivatives of

VHes = []

index = 0
oddball = 0

# index tracks how many entries there are, oddball tracks the parameter/parameter entries.
#  This distinction is necessary since each time step will have its own entries, but
#  parameter/parameter derivatives only have one Hessian entry.

# Fill out for constraint 1 first
# Note that this is a symmetrical matrix - filling out lower diagonal-half only.

for j in range(len(Sall)):
# Sall is Sv + Sk + Sp, the symbolic representations of the variables, controls, and parameters.
   for k in range(j+1):
     if k < (len(Sv) + len(Sk)):
       for n in [0,2]:
         H = sHes[0][j][k][n]
         if H != '0':
             temp1 = []
             temp1.append(index)
             temp1.append(0)
             temp1.append(j)
             temp1.append(k)
             temp1.append(n)
             temp1.append(H)
             temp1.append(1)
             if n == 0:
                temp1.append(sHes[0][j][k][1])
             VHes.append(temp1)
             index = index + 1
     else:  # These are the parameter/parameter derivatives
       H = sHes[0][j][k][0]
       if H != '0':
         temp2 = []
         temp1 = ''
         for n in [0,1,2]:
           H = sHes[0][j][k][n]
           if n == 0:
              temp1 = temp1 + H
           else:
              temp1 = temp1 + '+' + H
         temp2.append(index)
         temp2.append(0)
         temp2.append(j)
         temp2.append(k)
         temp2.append(-1)
         temp2.append(temp1)
         temp2.append(1)
         VHes.append(temp2)
         index = index + 1
         oddball = oddball + 1


# Fill out additional constraints, checking to see if row/column has been indexed already
# Need to make distinction for constraint/constraint derivatives
for i in range(len(AllCon)):
   for j in range(len(Sall)):
      for k in range(j+1):
        if k < (len(Sv) + len(Sk)):
          for n in [0,2]:
            H = sHes[i+1][j][k][n]
            new = 0
            if H != '0':
                # Check to see if row/column has been used before
                for u in range(len(VHes)):
                     if j == VHes[u][2]:
                         if k == VHes[u][3]:
                           if n == VHes[u][4]:
                               oldindex = VHes[u][0]
                               new = 1
                               temp1 = []
                               temp1.append(oldindex)
                               temp1.append(i+1)
                               temp1.append(j)
                               temp1.append(k)
                               temp1.append(n)
                               temp1.append(H)
                               temp1.append(0)
                               if n == 0:
                                temp1.append(sHes[i+1][j][k][1])
                if new == 1:
                   VHes.append(temp1) # Use old index
                elif new == 0: # Assign new index to new row/column combination
                   temp1 = []
                   temp1.append(index)
                   temp1.append(i+1)
                   temp1.append(j)
                   temp1.append(k)
                   temp1.append(n)
                   temp1.append(H)
                   temp1.append(1)
                   if n==0:
                    temp1.append(sHes[i+1][j][k][1])
                   VHes.append(temp1)
                   index = index + 1
        else: # Take care of the oddballs: parameter/parameter combinations
          H = sHes[i+1][j][k][0]
          new = 0
          if H != '0':
            for u in range(len(VHes)):
                if j == VHes[u][2]:
                    if k == VHes[u][3]:
                        oldindex = VHes[u][0]
                        new = 1
                        temp2 = []
                        temp1 = ''
                        for n in [0,1,2]:
                            H = sHes[i+1][j][k][n]
                            if n ==0:
                                temp1 = temp1 + H
                            else:
                                temp1 = temp1 + '+' + H
                        temp2.append(oldindex)
                        temp2.append(i+1)
                        temp2.append(j)
                        temp2.append(k)
                        temp2.append(-1)
                        temp2.append(temp1)
                        temp2.append(0)
            if new == 1:
                VHes.append(temp2)
            elif new == 0:
                temp2 = []
                temp1 = ''
                for n in [0,1,2]:
                    H = sHes[i+1][j][k][n]
                    if n == 0:
                        temp1 = temp1 + H
                    else:
                        temp1 = temp1 + '+' + H
                temp2.append(index)
                temp2.append(i+1)
                temp2.append(j)
                temp2.append(k)
                temp2.append(-1)
                temp2.append(temp1)
                temp2.append(1)
                VHes.append(temp2)
                index = index + 1
                oddball = oddball + 1

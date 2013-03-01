import sys
import os


ARGLIST = sys.argv
ARGLEN = len(ARGLIST)

derp = "(options) output_file_name integration_steps\n        -h for options\n        -u for usage and disclaimers\n"
options = "    -P path  : specify a standard PATH (including end /) for all files assuming standard names which is overwritten by other options\n\
    -e path  : specify an equations.txt to use.  Otherwise looks in current directory\n\
    -s path  : specify a specs file to use to find NSKIP.  Otherwise looks in current directory\n\
    -d path  : specify the IPOPT output path file to use (data.dat).  Otherwise looks in current directory\n\
    -p path  : specify the IPOPT output parameter file to use (param.dat).  Otherwise looks in current directory\n\
    -c path  : specify a directory to find stimulus files (including end /).  Otherwise looks in current directory\n\
    -r int   : specify the number of integration steps per model time step.  Default is currently 10\n"
usage = "This code is for forward integration of the estimates provided by IPOPT.  This code takes the end of the path estimate and parameter estimates as the starting point and integrates these forward using RK4.  This code requires an output file name and the number of time steps to integrate as input.  There are as of yet not much in the way of error messages, but most errors are due to files with an improper number of or order of entries.  If you successfully did a run with IPOPT the files used for that run and created by that run will work.  If a segfault occurs, the most common reason is stimulus files are not where you think they are or are not long enough.\n"


equations_file = "equations.txt"
specs_file = "specs.txt"
data_file = "data.dat"
param_file = "param.dat"
cpath = ''
PATH = ''
RESO = '10'

ispath=0
isequa=0
isspec=0
isdata=0
ispara=0
iscpat=0
isreso=0

argnum = 1

if ARGLEN==1:
   print derp
   sys.exit(0)

if ARGLIST[1] == '-h':
   print options
   sys.exit(0)

if ARGLIST[1] == '-u':
   print usage
   sys.exit(0)


while argnum < (ARGLEN -1):

   if ARGLIST[argnum] == '-e':
      equations_file_temp = ARGLIST[argnum+1]
      isequa=1
      argnum=argnum+2

   elif ARGLIST[argnum] == '-s':
      specs_file_temp = ARGLIST[argnum+1]
      isspec=1
      argnum=argnum+2

   elif ARGLIST[argnum] == '-d':
      data_file_temp = ARGLIST[argnum+1]
      isdata=1
      argnum=argnum+2

   elif ARGLIST[argnum] == '-p':
      param_file_temp = ARGLIST[argnum+1]
      ispara=1
      argnum=argnum+2

   elif ARGLIST[argnum] == '-P':
      PATH = ARGLIST[argnum+1]
      ispath=1
      argnum=argnum+2

   elif ARGLIST[argnum] == '-c':
      cpath_temp = ARGLIST[argnum+1]
      argnum=argnum+2

   elif ARGLIST[argnum] == '-r':
      RESO = ARGLIST[argnum+1]
      argnum=argnum+2

   else:
      break

outfile = ARGLIST[argnum]
NTIME = ARGLIST[argnum+1]

if ispath:
   equations_file = PATH+"equations.txt"
   specs_file = PATH+"specs.txt"
   data_file = PATH+"data.dat"
   param_file = PATH+"param.dat"
   cpath = PATH

if isequa:
   equations_file = equations_file_temp

if isdata:
   data_file = data_temp_file

if isspec:
   specs_file = specs_file_temp

if ispara:
   param_file = param_file_temp

if iscpat:
   cpath = cpath_temp



file = open(equations_file,"r")

temp=[]
for line in file:
  if line.startswith('#'): # Pound used as comment in text file
     continue
  elif line.startswith('\\'): # In case file has UTF-8 markers
     continue
  elif line=='':
     continue
  else:
     temp.append(line)

file.close()

h=[]
for i in range(len(temp)):
    temp1=temp[i].rstrip( )
    if temp1!='':
      h.append(temp1)



a=h[1].split(',')
NSTAT = int(a[0])
NPARA = int(a[1])
NCOUP = int(a[2])
NSTIM = int(a[3])
NFUNC = int(a[4])

if len(a)>5:
  NMEAS = int(a[5])
else:
  NMEAS=NCOUP



if len(h)!=(3+2*NSTAT+NPARA+NCOUP+NSTIM+NFUNC+NMEAS):
    print "improper equations.txt - check the declared number of variables with what are given\n"
    print 'number of non-comment lines in equations.txt       : '+str(len(h))+'\n'
    print 'number expected based on declared variable numbers : '+str(3+2*NSTAT+NPARA+NCOUP+NSTIM+NFUNC+NMEAS)+'\n'
    sys.exit()

FEQNSTR = []
VARSTR = []
for i in range(NSTAT):
    FEQNSTR.append(h[2+i])
    VARSTR.append(h[3+NSTAT+i])

OBJSTR = h[2+NSTAT]

PRMSTR = []
for i in range(NPARA):
    PRMSTR.append(h[3+2*NSTAT+i])

CPLSTR = []
for i in range(NCOUP):
    CPLSTR.append(h[3+2*NSTAT+NPARA+i])

MSRSTR = []
for i in range(NMEAS):
    MSRSTR.append(h[3+2*NSTAT+NPARA+NCOUP+i])

STMSTR = []
for i in range(NSTIM):
    STMSTR.append(h[3+2*NSTAT+NPARA+NCOUP+NMEAS+i])


OUTFUNC = []
for i in range(NSTAT):
    temp = FEQNSTR[i]
    for j in range(NSTAT):
        temp = temp.replace(VARSTR[j],"x["+str(j)+"]")
    for j in range(NPARA):
        temp = temp.replace(PRMSTR[j],"p["+str(j)+"]")
    for j in range(NCOUP):
        temp = temp.replace(CPLSTR[j],'0')
    for j in range(NMEAS):
        temp = temp.replace(MSRSTR[j],'0')
    for j in range(NSTIM):
        temp = temp.replace(STMSTR[j],'stim['+str(j)+']')
    OUTFUNC.append(temp)

file = open(param_file,"r")

temp=[]
for line in file:
  if line.startswith('#'): # Pound used as comment in text file
     continue
  elif line.startswith('\\'): # In case file has UTF-8 markers
     continue
  else:
     temp.append(line)

file.close()

PLIST=[]
for i in range(len(temp)):
    temp1=temp[i].rstrip( )
    a=temp1.split('\t')
    if temp1!='':
      PLIST.append(a[3])
      print a[3]

if len(PLIST)!=NPARA:
  print 'error : number of parameter values in param.dat does not equal NP in equations.txt\n\n   aborting\n'
  sys.exit(0)

file = open(data_file,'r')

#lazy way to get last line
time=0
for line in file:
    temp = line
    time = time+1

#print time

file.close()

temp = temp.rstrip( )
a = temp.split() #no argument in split assumes any amount of whitespace is a delimiter

if len(a)!=(NSTAT+NCOUP+NMEAS+1):
   print 'warning : data file appears to be from a source other than that generated by IPOPT\n'
   print '          be sure the data file you are using corresponds to the correct model.\n'

XLIST = []
for i in range(NSTAT):
   XLIST.append(a[i+1])


file = open(specs_file,'r')
temp=[]
for line in file:
  if line.startswith('#'): # Pound used as comment in text file
     continue
  elif line.startswith('\\'): # In case file has UTF-8 markers
     continue
  else:
     temp.append(line)

file.close()

h=[]
for i in range(len(temp)):
  h.append(temp[i].rstrip( ))

NSKIP = h[1]
time = time+int(NSKIP)
DT = str(float(h[2])*0.5)

ISTIMFILE=[]
for i in range(NSTIM):
   ISTIMFILE.append(h[3+NMEAS+i])


#uncomment these if errors occurm it may help finding them
   
#print XLIST
#print PLIST
#print ISTIMFILE
#print time
#sys.exit(0)

########################################
#######  make vector field      ########
########################################

file = open("RKTEMP.cpp","w")

file.write(
'#include <cmath>\n'
'#include <iostream>\n'
'#include <fstream>\n'
'#include <string>\n'
'#include <iomanip>\n'
'#include <sstream>\n'
)

if NFUNC!=0:
    file.write(
'#include \"'+PATH+'myfunctions.cpp\"\n'
)

file.write(
'#define RESO '+RESO+'\n'
'#define NSTAT '+str(NSTAT)+'\n'
'#define NPARA '+str(NPARA)+'\n'
'#define NSTIM '+str(NSTIM)+'\n'
'#define NTIME RESO*'+NTIME+'\n'
'#define NMEAS '+str(NMEAS)+'\n'
'#define DT '+DT+'\n'
'#define dt DT/((double) RESO)\n'

'using namespace std;\n'

'void Func(double *dx, double *x, double *p, double *stim);\n\n'
'int main(int argc, char **argv){\n\n'
'    double *K1 = new double[NSTAT];\n'
'    double *K2 = new double[NSTAT];\n'
'    double *K3 = new double[NSTAT];\n'
'    double *K4 = new double[NSTAT];\n'
'    double *X = new double[NSTAT];\n'
'    double *STIM = new double[2*NSTIM];\n'
'    double *Stimstep = new double[NSTIM];\n'
'    double *STIMTEMP = new double[3*NSTIM];\n'
'    double *P = new double[NPARA];\n'
'    int i,j,k;\n'
'    double Temp[NSTAT];\n'
'    const int WIDTH=13;\n'
'    const int PRECISION=6;\n\n'
)

j= 0
for i in PLIST:
    file.write(
'    P['+str(j)+'] = '+i+';\n'
    )
    j+=1

j= 0
for i in XLIST:
    file.write(
'    X['+str(j)+'] = '+i+';\n'
    )
    j+=1


file.write(
'    ifstream stim[NSTIM];\n'
)
for i in range(len(ISTIMFILE)):
     file.write(
'    stim['+str(i)+'].open(\"'+cpath+ISTIMFILE[i]+'\");\n'
#'    cout<<\"'+cpath+ISTIMFILE[i]+'\"<<endl;\n'
)
file.write(
'    ofstream *outfile=new ofstream;\n'
'    outfile->open(\"'+outfile+'\");\n'
'    (*outfile)<<setw(WIDTH)<<setprecision(PRECISION)<<'+str(time)+';\n'
'    (*outfile)<<setw(WIDTH)<<setprecision(PRECISION)<<'+str(time)+'*DT;\n'
'    for(i=0;i<NSTAT;i++){\n'
'        (*outfile)<<setw(WIDTH)<<setprecision(PRECISION)<<X[i];\n'
'    }\n\n'
'    (*outfile)<<endl;\n'
)
if NSTIM !=0:
    file.write(
#'    cout<<NSTIM*('+str(time)+')<<endl;\n'
'    for(j=0;j<NSTIM*('+str(time)+');j++){\n' #this should skip until the right place in the stimulus files
'        STIM[j%NSTIM+NSTIM] = STIM[j%NSTIM];\n'
'        stim[j%NSTIM]>>STIM[j%NSTIM];\n\n'
'    }\n'
'    for(j=0;j<NSTIM;j++){\n'
'        Stimstep[j] = (STIM[j+NSTIM]-STIM[j])/((double) RESO);\n'
'        STIMTEMP[j] = STIM[j];\n'
'        STIMTEMP[j+NSTIM] = STIM[j] + .5*Stimstep[j];\n'
'        STIMTEMP[j+2*NSTIM] = STIM[j] + Stimstep[j];\n'
'    }\n\n'
    )
file.write(
'    for(i=0;i<NTIME;i++){\n'
'        if(X[0] != X[0]){\n'
'            cout<<"NAN error"<<endl;\n'
'            break;\n'
'        }\n'
'        Func(K1,X,P,&STIMTEMP[0]);\n'
'        for(j=0;j<NSTAT;j++)\n'
'            Temp[j] = X[j]+.5*dt*K1[j];\n\n'
'        Func(K2,Temp,P,&STIMTEMP[NSTIM]);\n'
'        for(j=0;j<NSTAT;j++)\n'
'            Temp[j] = X[j]+.5*dt*K2[j];\n\n'
'        Func(K3,Temp,P,&STIMTEMP[NSTIM]);\n'
'        for(j=0;j<NSTAT;j++)\n'
'            Temp[j] = X[j]+dt*K3[j];\n\n'
'        Func(K4,Temp,P,&STIMTEMP[2*NSTIM]);\n'
'        for(j=0;j<NSTAT;j++)\n'
'            X[j] = X[j]+1.0/6.0*dt*(K1[j]+2.0*K2[j]+2.0*K3[j]+K4[j]);\n\n'
'        if((i+1)%RESO==0){\n'
)
#if NSTIM != 0:
file.write(
'            for(j=0;j<NSTIM;j++){\n'
'                STIM[j]=STIM[j+NSTIM];\n'
'                stim[j]>>STIM[j+NSTIM];\n'
'                Stimstep[j] = (STIM[j+NSTIM]-STIM[j])/((double) RESO);\n'
'            }\n'
    )
file.write(
'            (*outfile)<<setw(WIDTH)<<setprecision(PRECISION)<<'+str(time)+'+i/RESO;\n'
'            (*outfile)<<setw(WIDTH)<<setprecision(PRECISION)<<('+str(time)+'+i/RESO)*DT;\n'
'            for(j=0;j<NSTAT;j++)\n'
'                (*outfile)<<setw(WIDTH)<<setprecision(PRECISION)<<X[j];\n\n'
'                //cout<<(i+1)/RESO<<endl;\n'
'            (*outfile)<<endl;\n'
'        }\n'
)
if NSTIM != 0:
    file.write(
'        for(j=0;j<NSTIM;j++){\n'
'            STIMTEMP[j]=STIMTEMP[j+2*NSTIM];\n'
'            STIMTEMP[j+NSTIM]=STIMTEMP[j]+.5*Stimstep[j];\n'
'            STIMTEMP[j+2*NSTIM]=STIMTEMP[j]+Stimstep[j];\n'
'        }\n'
    )
file.write(
'    }\n'
'    return 0;\n'
'}\n'
)

file.write(
'void Func(double *dx, double *x, double *p, double *stim){\n'
)
for i in range(NSTAT):
    file.write(
'    dx['+str(i)+'] = '+OUTFUNC[i]+';\n'
    )
file.write('}\n\n')
file.close()



os.system("g++ -o theproblem -lm RKTEMP.cpp")
os.system("./theproblem")
os.system("rm RKTEMP.cpp")
os.system("rm theproblem")


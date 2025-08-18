import subprocess
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import os
import sys

rc('mathtext', default='regular')
data = []
temps = []
isps = []
ox_perc = []
MWs = []
O_Fs=[]
c_stars=[]
gammas=[]

for i in np.linspace(0.01,100,200):
    ox_perc.append(i)
    inp = """
problem   
    rocket  frozen  nfz=1 
  p,bar=400,
  pi/p=400,
react  
  fuel=CH4(L) wt=$CH4$  t,k=112  
  oxid=O2(L) wt=$O2$  t,k=90  
output  
    plot t isp 
end
    """
    inp = inp.replace("$CH4$", str(i))
    inp = inp.replace("$O2$", str(100-i))
    with open("testcomp.inp", "w") as f:
        f.write(inp)
    process = subprocess.Popen("FCEA2.exe", stdin=subprocess.PIPE, shell=True)
    time.sleep(0.1)
    process.stdin.write(b"testcomp\n")
    process.stdin.flush()
    stdout, stderr = process.communicate()

    with open("testcomp.out", "r") as f:
        lines = f.readlines()
        temp = None
        isp = None
        gamma = None
        cstar = None
        MolW = None
        o_f=None
        for line in lines:
            if ("T, K" in line):
                words = line.split()
                temp = float(words[2])
            elif ("Isp," in line):
                words = line.split()
                try:
                    isp = float(words[3])
                except:
                    pass
            elif ("GAMMAs" in line):
                words = line.split()
                gamma = float(words[1])
            elif ("CSTAR" in line):
                words = line.split()
                try:
                    cstar = float(words[2])
                except:
                    pass
            elif ("M, (1/n)" in line):
                words = line.split()
                try:
                    MolW = float(words[3])
                except:
                    pass
            elif ("O/F" in line):
                words = line.split()
                try:
                    o_f = float(words[1])
                except:
                    pass
        # data.append([i, temp, isp])
        temps.append(temp)
        MWs.append(MolW)
        O_Fs.append(o_f)
        c_stars.append(cstar)
        gammas.append(gamma)
        if isp is not None:
            isps.append(isp)
        else:
            try:
                isps.append((2*gamma**2/(gamma-1)*((2/(gamma+1))**((gamma+1)/(gamma-1)))*(1-0.01**(1-1/gamma)))**0.5*cstar)
            except:
                isps.append(None)
print(O_Fs)
print(temps)
print(isps)
print(MWs)
print(c_stars)
'''
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ox_perc, isps, '-', label = 'Isp')
# ax.plot(time, Rn, '-', label = 'Rn')
ax2 = ax.twinx()
ax2.plot(ox_perc, temps, '-r', label = 'Tf')
ax.grid()
ax.set_xlabel("Oxidiser percentage")
ax.set_ylabel(r"Isp (m/s)")
ax2.set_ylabel(r"Temperature (K)")
ax2.set_ylim(min(x for x in temps if x is not None)-100, max(x for x in temps if x is not None)+100)
ax.set_ylim(min(x for x in isps if x is not None)-100, max(x for x in isps if x is not None)+100)
ax.legend(loc='upper left')
ax2.legend(loc='upper right')
plt.show()
'''
#os.makedirs("TRC/LOXMETHANE.txt", exist_ok=True)
#assert os.path.isfile("TRC/LOXMETHANE.txt")
sys.stdout = open("LOXMETHANE400bar.txt", 'w')
print(f"O/F\ttemp\tisp\tmw\tcstar\tgamma")
for i in range(len(O_Fs)):
    print(f"{O_Fs[i]}\t{temps[i]}\t{isps[i]}\t{MWs[i]}\t{c_stars[i]}\t{gammas[i]}")
sys.stdout.close()
sys.stdout = sys.__stdout__
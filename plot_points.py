import matplotlib.pyplot as plt
import csv

x = []
y = []

with open('points.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))
data= ""
#plt.figure(figsize=(8, 6), dpi=400)
plt.scatter(x,y,alpha=0.5,s=0.5)
plt.xlabel(r'$\tau$')
plt.ylabel('d (dose size)')
plt.title('Aceptable vs unaceptable area')
#plt.legend(loc=1)
data= "IC={H[0]==10^-4 /0.005,CTL[0]==10^-4/0.005,M[0]==0.1,Den[0]==0,IL2[0]==0}\nUpper bound = 0.8 * Tumor carrying capacity "
plt.figtext(0.5, 0.07, data, horizontalalignment='center',
            fontsize=9, multialignment='left',
            bbox=dict(boxstyle="round", facecolor='#D8D8D8',
                      ec="0.5", pad=0.5, alpha=1), fontweight='bold')
plt.show()
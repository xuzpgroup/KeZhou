import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

# plt.rcParams["font.family"]="Arial"

file_name="gas.tar.dat"
file_gr="G.bck.xyz"

# file_name="tt.dat"
matrix=np.genfromtxt(file_name)
matrix_OX=np.genfromtxt(file_gr)

Lx=63.00000000000000
Ly=111.5440000000000
x=matrix[:, 1]
y=matrix[:, 2]

xlo=min(x)
ylo=min(y)
xhi=max(x)
yhi=max(y)


gx=matrix_OX[:, 1]
gx1=matrix_OX[:, 1]-Lx
gy=matrix_OX[:, 2]


# Plot data
# fig1 = plt.figure()
# plt.plot(x,y,'.r')
# plt.xlabel('x')
# plt.ylabel('y')

# make a figure twice as tall as it is wide
w, h = plt.figaspect(0.6)
fig2 = plt.figure(figsize=(w,h))
#fig2 = plt.figure()

# plt.hist2d(x, y, bins=100, norm=LogNorm(), vmin=1, vmax=1.0E4)
plt.hist2d(x, y, bins=200, norm=LogNorm(), cmap=plt.cm.Greens)
# plt.hexbin(x, y, gridsize=128, norm=LogNorm(), cmap=plt.cm.Greens)
 
# plot background 
for i in range(-2,2):
	for j in range(-2,2):
		plt.plot(gx+Lx*i,gy+Ly*j,'or', ms=2)


# plt.plot(gx,gy,'or', ms=2)
# plt.plot(gx1,gy,'or', ms=2)

# plt.plot(ox1,oy,'og', ms=8)
# plt.colorbar()
plt.xlim(xlo-5, xhi+5)
plt.ylim(ylo-5, yhi+5)

plt.xlabel('x', fontsize=15)
plt.ylabel('y', fontsize=15)
cbar = plt.colorbar()
# cbar.set_clim(-2.0, 2.0)
cbar.ax.set_ylabel('Counts')
# plt.savefig('2Dhis.eps', format='eps', dpi=300)
plt.savefig('2Dhis.png', format='png', dpi=300)
plt.show()



# Estimate the 2D histogram
# nbins = 50
# H, xedges, yedges = np.histogram2d(x,y,bins=nbins)

# # H needs to be rotated and flipped
# H = np.rot90(H)
# H = np.flipud(H)

# # Mask zeros
# Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero

# # Plot 2D histogram using pcolor
# fig2 = plt.figure()
# plt.pcolormesh(xedges,yedges,Hmasked)
# plt.xlabel('x')
# plt.ylabel('y')
# cbar = plt.colorbar()
# cbar.ax.set_ylabel('Counts')
# plt.show()
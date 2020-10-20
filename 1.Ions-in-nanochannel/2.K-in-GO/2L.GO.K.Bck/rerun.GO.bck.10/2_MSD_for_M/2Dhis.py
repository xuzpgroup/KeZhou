import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

# Create some random numbers

file_name="gas.tar.dat"
file_OX1="OX1"
file_OX2="OX2"

# file_name="tt.dat"
matrix=np.genfromtxt(file_name)
matrix_OX1=np.genfromtxt(file_OX1)
matrix_OX2=np.genfromtxt(file_OX2)


x=matrix[:, 1]
y=matrix[:, 2]

ox12=matrix_OX1[:, 1]-63
ox11=matrix_OX1[:, 1]
oy11=matrix_OX1[:, 2]
oy12=matrix_OX1[:, 2]-111.544


ox22=matrix_OX2[:, 1]-63
ox21=matrix_OX2[:, 1]
oy21=matrix_OX2[:, 2]
oy22=matrix_OX2[:, 2]-111.544


# Plot data
# fig1 = plt.figure()
# plt.plot(x,y,'.r')
# plt.xlabel('x')
# plt.ylabel('y')

# make a figure twice as tall as it is wide
w, h = plt.figaspect(0.81)
fig2 = plt.figure(figsize=(w,h))
#fig2 = plt.figure()

# plt.hist2d(x, y, bins=100, norm=LogNorm(), vmin=1, vmax=1.0E4)
plt.hist2d(x, y, bins=100, norm=LogNorm())
plt.plot(ox11,oy11,'og', ms=8)
plt.plot(ox11,oy12,'og', ms=8)
plt.plot(ox12,oy12,'og', ms=8)
plt.plot(ox12,oy11,'og', ms=8)

plt.plot(ox21,oy21,'ob', ms=8)
plt.plot(ox21,oy22,'ob', ms=8)
plt.plot(ox22,oy22,'ob', ms=8)
plt.plot(ox22,oy21,'ob', ms=8)

# plt.colorbar()
plt.xlim(20.0, 70.0)
plt.ylim(0.0, 60)
plt.xlabel('x', fontsize=15)
plt.ylabel('y', fontsize=15)
cbar = plt.colorbar()
# cbar.set_clim(-2.0, 2.0)
cbar.ax.set_ylabel('Counts')
plt.savefig('2Dhis.K.png', format='png')
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
plt.show()

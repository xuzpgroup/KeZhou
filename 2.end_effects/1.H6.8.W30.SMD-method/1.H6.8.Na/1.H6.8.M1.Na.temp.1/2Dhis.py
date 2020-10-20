import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

# Create some random numbers

file_name="pos.data"
# file_OX="OX"

# file_name="tt.dat"
matrix=np.genfromtxt(file_name)
# matrix_OX=np.genfromtxt(file_OX)

x=matrix[:, 1]
y=matrix[:, 3]

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

# plt.colorbar()
plt.xlim(15.0, 40.0)
plt.ylim(0.0, 20.0)
plt.xlabel('x', fontsize=15)
plt.ylabel('y', fontsize=15)
cbar = plt.colorbar()
# cbar.set_clim(-2.0, 2.0)
cbar.ax.set_ylabel('Counts')
plt.savefig('2Dhis.ions.tar.png', format='png')
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

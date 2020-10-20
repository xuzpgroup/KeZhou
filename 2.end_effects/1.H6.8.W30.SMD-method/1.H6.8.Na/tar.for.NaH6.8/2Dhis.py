import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.ticker import  MultipleLocator
from matplotlib.ticker import  FormatStrFormatter
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
w, h = plt.figaspect(0.45)
fig2 = plt.figure(figsize=(w,h))
#fig2 = plt.figure()
xmajorLocator = MultipleLocator(5.0) 
xmajorFormatter = FormatStrFormatter('%.0f') 
ymajorLocator = MultipleLocator(5.0) 
ymajorFormatter = FormatStrFormatter('%.0f')


# plt.hist2d(x, y, bins=100, norm=LogNorm(), vmin=1, vmax=1.0E4)
ax=plt.subplot(111)
plt.hist2d(x, y, bins=100, norm=LogNorm())
# plt.hist2d(x, 20.0-y, bins=100, norm=LogNorm())


ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_major_formatter(xmajorFormatter)
ax.yaxis.set_major_locator(ymajorLocator)
ax.yaxis.set_major_formatter(ymajorFormatter)

# plt.colorbar()
# plt.xlim(15.0, 50.0)
# plt.ylim(0.0, 20.0)
plt.xlabel('x', fontsize=15)
plt.ylabel('z', fontsize=15)
cbar = plt.colorbar()
# cbar.set_clim(-2.0, 2.0)
cbar.ax.set_ylabel('Counts')
plt.savefig('2Dhis.ions.tar.png', format='png')
# plt.show()



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

from spectral import *
import matplotlib.pyplot as plt
import numpy as np

# Replace 'your_envi_file' with the actual path to your ENVI header file
envi_header_file1 = '/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2/ATBD/atbd.hdr'
envi_data_file1 = '/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2/ATBD/atbd.envi'

envi_header_file2 = '/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2/ATBD/atbd_combine.hdr'
envi_data_file2 = '/home/corroyez/Documents/NC_Full/03_RESULTS/Mormal/LAI_PAI_Masks/LAI_S2/ATBD/atbd_combine.envi'

# Read the ENVI file
image1 = envi.open(envi_header_file1, image=envi_data_file1)
image2 = envi.open(envi_header_file2, image=envi_data_file2)

# Accessing spectral data
spectral_data1 = image1.load()
spectral_data2 = image2.load()

# Flatten the 3D array to a 1D array for histogram
flat_data1 = spectral_data1.flatten()
flat_data2 = spectral_data2.flatten()

# Remove any NaN or NoData values
flat_data1 = flat_data1[~np.isnan(flat_data1)]
flat_data2 = flat_data2[~np.isnan(flat_data2)]

# Plot histogram
plt.hist(flat_data1, bins=200, color='green', alpha=0.5)
plt.hist(flat_data2, bins=200, color='red', alpha=0.5)

plt.title('Histogram of ENVI File')
plt.xlabel('Pixel Values')
plt.ylabel('Frequency')
plt.legend()
plt.show()
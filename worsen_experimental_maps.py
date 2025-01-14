import numpy as np
import matplotlib.pyplot as plt
import mrcfile
import os
from utils import read_density_data_mrc, delete_extension_from_filename

experimental_data = os.path.join(os.getcwd(), "experimental_data")
input_density_name = "6ply_Ligand_Map_scaled_boxed_pixel_0.5_rescaled_For_Elisei.mrc"
input_density_path = os.path.join(experimental_data, input_density_name)

with mrcfile.open(input_density_path, "r") as mrc_file:
    input_density = mrc_file.data
    input_header = mrc_file.header
    input_voxel_size = mrc_file.voxel_size

# std = np.std(input_density)
# print(f"Std: {std}")
std = 2.0

data1 = np.random.normal(loc=1.5*std, scale=std, size=10000)
data2 = np.random.normal(loc=-1.5*std, scale=std, size=10000)
data = np.concatenate([data1, data2])
hist, bins = np.histogram(data, bins=200)

bin_midpoints = bins[:-1] + np.diff(bins)/2
cdf = np.cumsum(hist)
cdf = cdf / cdf[-1]
values = np.random.rand(input_density.shape[0] * input_density.shape[1] * input_density.shape[2])
value_bins = np.searchsorted(cdf, values)
noise_values = bin_midpoints[value_bins]
noise = np.array(noise_values.reshape(input_density.shape), dtype=input_density.dtype)
noise[input_density <= 0.0] = 0.0

output_density_name = f"{delete_extension_from_filename(input_density_name)}_worsen.mrc"
output_density_path = os.path.join(experimental_data, output_density_name)
with mrcfile.new(output_density_path, overwrite=True) as mrc_file:
        mrc_file.set_data(input_density + noise)
        mrc_file.voxel_size = input_voxel_size
        mrc_file.header.origin = input_header.origin



plt.subplot(121)
plt.hist(data, 100)
plt.subplot(122)
plt.hist(noise_values, 100)
plt.show()



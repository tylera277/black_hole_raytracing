
#


# Imports
from PIL import Image
from matplotlib import pyplot as plt

import pandas as pd
import numpy as np


##############################################
# Main Code Body

img = Image.open('../yosemite.jpeg')
img_array = np.array(img).astype(int)

img_height = img_array.shape[0]
img_width = img_array.shape[1]

print("PICTURE HEIGHT: ", img_height)
print("PICTURE WIDTH: ", img_width)


angle_data = pd.read_csv('../outputs/angles_rkf4_openmp_JULY.csv')

#angle_data = pd.read_csv('../angles.csv')

angle_data = angle_data.apply(pd.to_numeric)

# Shift the data points to within the measurements respective angles
angle_data['theta_prime_adjusted'] = angle_data[' theta_prime'] % 3.14159
angle_data['phi_prime_adjusted'] = angle_data[' phi_prime'] % (2*3.14159)


angle_data['camera_pixel_from_bottom'] = angle_data['theta_cs'] / (3.14159/img_height)
angle_data['camera_pixel_from_left'] = angle_data[' phi_cs'] / ((2*3.14159)/img_width)


# angle_data['camera_pixel_from_left'] = angle_data['camera_pixel_from_left'] - 2050
# angle_data['camera_pixel_from_bottom'] = angle_data['camera_pixel_from_bottom'] - 2080

angle_data['warped_pixel_from_bottom'] = angle_data['theta_prime_adjusted'] / (3.14159/img_height)
angle_data['warped_pixel_from_left'] = angle_data['phi_prime_adjusted'] / ((2*3.14159)/img_width)


# 1025, 2080
angle_data_cleaned = angle_data.dropna(axis=0)

print(angle_data.columns)


angle_data_cleaned = angle_data_cleaned.astype({"warped_pixel_from_left": int, "warped_pixel_from_bottom": int,
                        "camera_pixel_from_left":int, "camera_pixel_from_bottom":int})

print(angle_data_cleaned.head)

warped_picture = np.zeros((img_height, img_width, 3))

print("MAX1:", angle_data_cleaned['camera_pixel_from_left'].max())
print("Min1:", angle_data_cleaned['camera_pixel_from_left'].min())

print("MAX2:", angle_data_cleaned['camera_pixel_from_bottom'].max())
print("Min2:", angle_data_cleaned['camera_pixel_from_bottom'].min())

print(angle_data_cleaned.columns)

# angle_data_cleaned = angle_data_cleaned.to_csv("../angles_rkf4_openmp_output.csv")

print(angle_data_cleaned.head())
angle_data_cleaned = angle_data_cleaned.to_numpy().astype(int)



vertical_theta = 3.14159/2.0
horizontal_phi = 3.14159/2.0
counter = 0
max_value = 255
horizontal_start = 1025
vertical_start = 2081
warped_picture[:][:] = [255, 255, 255]
print("GERANIMO:", img_array.shape)
for i in range(1, 200*200, 1):

    if (angle_data_cleaned[i][2] == 0) and (angle_data_cleaned[i][3] == 0):
        print("YEPPER ", counter)
        counter += 1
        warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][0] = 0

        warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][1] = 0

        warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][2] = 0

    else:
        warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][0] = img_array[angle_data_cleaned[i][8]][angle_data_cleaned[i][9]][0]

        warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][1] = img_array[angle_data_cleaned[i][8]][angle_data_cleaned[i][9]][1]

        warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][2] = img_array[angle_data_cleaned[i][8]][angle_data_cleaned[i][9]][2]


warped_picture = warped_picture.astype(int)

plt.imshow(warped_picture, interpolation='nearest')
plt.show()




from PIL import Image
from matplotlib import pyplot as plt


import pandas as pd

import numpy as np


img = Image.open('../star_background.tiff')
img_array = np.array(img).astype(int)
print(img_array[1][15])
img_width = img_array.shape[0]
img_height = img_array.shape[1]

print(img_width, img_height)


angle_data = pd.read_csv('../angles_rkf4.csv')
#angle_data = pd.read_csv('../angles.csv')

print(angle_data.columns)

angle_data['theta_prime_adjusted'] = angle_data[' theta_prime'] % 3.14159
angle_data['phi_prime_adjusted'] = angle_data[' phi_prime'] % (2*3.14159)

angle_data['camera_pixel_from_left'] = angle_data[' phi_cs'] / (2*3.14159/4101)
angle_data['camera_pixel_from_bottom'] = angle_data['theta_cs'] / (3.14159/4162)

# angle_data['camera_pixel_from_left'] = angle_data['camera_pixel_from_left'] - 1025
#angle_data['camera_pixel_from_bottom'] = angle_data['camera_pixel_from_bottom'] - 2080

angle_data['warped_pixel_from_left'] = angle_data['phi_prime_adjusted'] / ((2*3.14159)/4101)
angle_data['warped_pixel_from_bottom'] = angle_data['theta_prime_adjusted'] / (3.14159/4162)


print(angle_data.head)
# 1025, 2080
angle_data_cleaned = angle_data.dropna(axis=0)

angle_data_cleaned = angle_data_cleaned.astype({"warped_pixel_from_left": int, "warped_pixel_from_bottom": int,
                        "camera_pixel_from_left":int, "camera_pixel_from_bottom":int})
print("length:", len(angle_data))
print(angle_data_cleaned['camera_pixel_from_bottom'])

warped_picture = np.zeros((4200, 4200, 3))
print(angle_data.columns)
angle_data_cleaned = angle_data_cleaned.to_numpy().astype(int)

print("BEEP:", type(angle_data_cleaned[0][9]))


vertical_theta = 3.14159/2.0
horizontal_phi = 3.14159/2.0
counter = 0
max_value=255
horizontal_start = 1025
vertical_start = 2081
warped_picture[:][:] = [255, 255, 255]

for i in range(1, 40400, 1):
    warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][0] = \
        img_array[angle_data_cleaned[i][8]][angle_data_cleaned[i][9]][0]

    warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][1] = \
        img_array[angle_data_cleaned[i][8]][angle_data_cleaned[i][9]][1]

    warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][2] = \
        img_array[angle_data_cleaned[i][8]][angle_data_cleaned[i][9]][2]

    #warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][0] = 0
    #warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][1] = 0
    #warped_picture[angle_data_cleaned[i][6]][angle_data_cleaned[i][7]][2] = 0

    if img_array[angle_data_cleaned[i][8]][angle_data_cleaned[i][9]][0] < max_value:
        max_value = img_array[angle_data_cleaned[i][8]][angle_data_cleaned[i][9]][0]

for vertical in range(0, 4162, 1):
    for horizontal in range(0, 4101, 1):

        pass
#        if warped_picture[vertical][horizontal][:] != 0:
#            pass
#        else:
#            warped_picture[vertical][horizontal][0] = 0
print("MAX:", max_value)


warped_picture = warped_picture.astype(int)

plt.imshow(warped_picture, interpolation='nearest')
plt.show()
print("CHECK:", warped_picture[0][0])

while horizontal_phi < ((3.14159/2.0)+200*((2*3.14159)/img_width)):
    vertical_theta = 3.14159/2.0
    while vertical_theta < ((3.14159/2.0)+200*(3.14159)/img_height):






        vertical_theta += ((3.14159/img_height))

    horizontal_phi += ((2*3.14159)/img_width)

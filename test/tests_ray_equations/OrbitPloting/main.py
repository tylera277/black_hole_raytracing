

import pandas as pd
import matplotlib.pyplot as plt

orbit_positions = pd.read_csv("../output.csv")

print(orbit_positions.head())

orbit_positions = orbit_positions.to_numpy()
plt.ion()
figure, ax = plt.subplots(figsize=(10, 8))
plt.xlim(-5e11, 5e11)
plt.ylim(-5e11, 5e11)

for position_value in range(0, orbit_positions.shape[0], 100):
    plt.scatter(orbit_positions[position_value][0], orbit_positions[position_value][1])
    plt.pause(0.01)

plt.show()
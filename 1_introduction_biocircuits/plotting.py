import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the CSV file
file_path = 'transcription_translation_sundials.csv'
data = pd.read_csv(file_path)

# Plotting both mRNA concentration and Protein concentration on the same plot for direct comparison
plt.figure(figsize=(10, 6))

# Plotting mRNA concentration
plt.plot(data['Time'], data['mRNA_concentration'], label='mRNA Concentration', color='blue')

# Plotting Protein concentration
plt.plot(data['Time'], data['Protein_concentration'], label='Protein Concentration', color='green')

# Setting labels and title
plt.xlabel('Time (units)')
plt.ylabel('Concentration')
plt.title('mRNA and Protein Concentration vs Time')

# Adding legend and grid
plt.legend()
plt.grid(True)

# Display the plot
plt.show()

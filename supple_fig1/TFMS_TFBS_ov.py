import os
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

# == Plot Configuration
plt.rcParams['font.family'] = 'Arial'

# == Working Directory Paths
wdir = "motif_peak_overlap"

# == Read Peak Counts into Dictionary
p_dic = {}
with open(os.path.join(wdir, "peak_files_counting.txt"), 'r') as f_peak:
    for line in f_peak:
        count, filename = line.split()
        p_dic[filename] = int(count)

# == Read Motif Counts into Dictionary
m_dic = {}
with open(os.path.join(wdir, "motif_files_counting.txt"), 'r') as f_motif:
    for line in f_motif:
        count, filename = line.split()
        m_dic[filename] = int(count)

# == Initialize Lists for Data Collection
excel_tf = []
excel_tf_motif = []
excel_tf_peak = []
excel_m_ov_count = []
excel_m_total = []
excel_m_percent = []
excel_p_ov_count = []
excel_p_total = []
excel_p_percent = []

# == Process Overlap Files
for root, _, files in os.walk(wdir):
    for file in tqdm(files, desc="Processing files"):
        if file.endswith('.txt'):
            excel_tf.append(file.split('_')[0])
            with open(os.path.join(root, file), 'r') as ov_file:
                m_loci = set()
                p_loci = set()
                for line in ov_file:
                    ls = line.split('\t')
                    m_loci.add(f"{ls[0]}-{ls[1]}-{ls[2]}")
                    p_loci.add(f"{ls[6]}-{ls[7]}-{ls[8]}")

            # Count unique motif and peak loci
            m_count = len(m_loci)
            p_count = len(p_loci)
            excel_m_ov_count.append(m_count)
            excel_p_ov_count.append(p_count)

            # Extract file names for motif and peak
            file_name_split = file.split('.bed')
            excel_tf_motif.append(file_name_split[0])
            excel_tf_peak.append(file_name_split[1][1:])

            # Get total counts from dictionaries
            m_total = m_dic[file_name_split[0] + '.bed']
            p_total = p_dic[file_name_split[1][1:] + '.bed']
            excel_m_total.append(m_total)
            excel_p_total.append(p_total)

            # Calculate overlap percentages
            excel_m_percent.append(m_count / m_total)
            excel_p_percent.append(p_count / p_total)

# == Create DataFrame and Save to CSV
output_df = pd.DataFrame({
    'TF': excel_tf,
    'Motif_file': excel_tf_motif,
    'Peak_file': excel_tf_peak,
    'Motif_ov_Peak': excel_m_ov_count,
    'Total motif': excel_m_total,
    '% motif ov peak': excel_m_percent,
    'Peak_ov_Motif': excel_p_ov_count,
    'Total peak': excel_p_total,
    '% peak ov motif': excel_p_percent
})

output_path = os.path.join(wdir, 'TF_motif_peak_overlap_all.csv')
output_df.to_csv(output_path, index=False)
print(f"Results saved to: {output_path}")

# == Moving Average Function
def moving_average(data, window_size):
    """Calculate the moving average of a list with a specified window size."""
    smoothed_data = []
    for i in range(len(data)):
        if i < window_size:
            smoothed_data.append(data[i])
        else:
            window = data[i - window_size:i]
            smoothed_data.append(sum(window) / window_size)
    return smoothed_data

# == Plotting Overlap Percentages
# Read processed data for plotting
df = pd.read_csv(os.path.join(wdir, 'TF_motif_peak_overlap_all.csv'))
m_percent = df['% motif ov peak'].tolist()
p_percent = df['% peak ov motif'].tolist()

# Sort data based on motif overlap percentage
m_sorted, p_sorted = zip(*sorted(zip(m_percent, p_percent), key=lambda x: x[0]))

# Smooth peak overlap data with a moving average
p_sorted_smooth = moving_average(p_sorted, window_size=10)

# Create figure with dual y-axes
fig, ax1 = plt.subplots(figsize=(10, 4))
ax2 = ax1.twinx()

# Plot motif overlap on the left y-axis
ax1.plot(m_sorted, color='orange', linestyle='-', label='% of motif overlap with peak')
ax1.set_ylim(0, 0.3)
ax1.set_ylabel('% of motif overlap with peak', fontsize=15)

# Plot peak overlap on the right y-axis
ax2.plot(p_sorted_smooth, color='green', linestyle='-', label='% of peak overlap with motif (smoothed)')
ax2.set_ylim(0, 1)
ax2.set_ylabel('% of peak overlap with motif', fontsize=15)

# Add legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', bbox_to_anchor=(1.1, 1))

# Save the plot
plt.tight_layout()
plot_path = os.path.join(wdir, 'TF_motif_peak_overlap.png')
plt.savefig(plot_path, dpi=300)
print(f"Plot saved to: {plot_path}")

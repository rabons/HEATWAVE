import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

# Color palette
colors = ["#93C0DB", "#D7E6EE", "#FEFECB","#FBEB99","#F4CC68","#EBA754","#E48650","#D1624C","#B8524A","#A44642","#8A3D3D","#72372E","#5A2E38","#40204D","#352045","#2B1B3D"]

# Create figure with multiple subplots
fig = plt.figure(figsize=(14, 8))

# 1. Horizontal color bar (gradient)
ax1 = plt.subplot(3, 1, 1)
ncolors = 256
cmap = LinearSegmentedColormap.from_list("Heatwave", colors, N=ncolors)
gradient = np.linspace(0, 1, 256).reshape(1, -1)
ax1.imshow(gradient, aspect='auto', cmap=cmap, extent=[0, 1, 0, 1])
ax1.set_title('Color Gradient (Continuous)', fontsize=14, fontweight='bold', pad=20)
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 1)
ax1.axis('off')

# 2. Individual color swatches with hex codes
ax2 = plt.subplot(3, 1, 2)
n = len(colors)
for i, color in enumerate(colors):
    rect = mpatches.Rectangle((i/n, 0), 1/n, 1, facecolor=color, edgecolor='black', linewidth=0.5)
    ax2.add_patch(rect)
    # Add hex code text
    ax2.text(i/n + 1/(2*n), 0.5, color, ha='center', va='center', 
             fontsize=7, rotation=90, color='white' if i > 2 else 'black',
             fontweight='bold')
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)
ax2.set_title('Individual Color Swatches with Hex Codes', fontsize=14, fontweight='bold', pad=20)
ax2.axis('off')

# 3. Vertical color bar
ax3 = plt.subplot(1, 3, 3)
gradient_vertical = np.linspace(0, 1, 256).reshape(-1, 1)
ax3.imshow(gradient_vertical, aspect='auto', cmap=cmap, extent=[0, 1, 0, 1])
ax3.set_title('Vertical Color Bar', fontsize=14, fontweight='bold', pad=20)
ax3.set_xlim(0, 1)
ax3.set_ylim(0, 1)
ax3.axis('off')

plt.tight_layout()
plt.savefig('color_palette_visualization.png', dpi=300, bbox_inches='tight')
print("Color palette visualization saved as 'color_palette_visualization.png'")
plt.show()


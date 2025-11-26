import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# %% Constructing shapes of HCS
# X-Y coordinates of HCS, fixed
x_HCS = np.linspace(-10, 10, 2001)

# Z coordinate of HCS, varies with situations
# Situation-1: without fold, without flapping
None
# Situation-2: without fold, steady flapping
y_HCS2 = np.zeros_like(x_HCS) - 0.3
# Situation-3: without fold, kink-like flapping
x_lambda3 = 10
y_HCS23 = 3 * np.cos(2 * np.pi / x_lambda3 * x_HCS)
# Situation-4: with fold, without flapping
x_lambda4 = 10
y_HCS4 = 3 * np.cos(2 * np.pi / x_lambda4 * x_HCS)
# Situation-5: with fold, steady flapping
x_lambda5 = 10
y_HCS5 = 3 * np.cos(2 * np.pi / x_lambda5 * x_HCS)
# Situation-6: with fold, kink-like flapping
x_lambda6 = 10
y_HCS6 = 3 * np.cos(2 * np.pi / x_lambda6 * x_HCS)

# %% Constructing orbits of PSP
# orbit of PSP, fixed
x_PSP = np.linspace(-10, 10, 2001)
y_PSP = np.zeros_like(x_PSP)
orb_PSP_beg = 1
orb_PSP_quiv = -2

# %% Determining crossing locations
# situation-1
None
# situation-2
x_cross2 = np.array([-3/2, -1/2, 1/2, 3/2]) / 2 * x_lambda3
y_cross2 = np.zeros_like(x_cross2) - 0.3
# situation-3
x_cross3 = np.array([-3/2, -1/2, 1/2, 3/2]) / 2 * x_lambda3
y_cross3 = np.zeros_like(x_cross3)
# situation-4
x_cross4 = np.array([-3/2, -1/2, 1/2, 3/2]) / 2 * x_lambda4
y_cross4 = np.zeros_like(x_cross4)
# situation-5
x_cross5 = np.array([-3/2, -1/2, 1/2, 3/2]) / 2 * x_lambda5
y_cross5 = np.zeros_like(x_cross5)
# situation-6
x_cross6 = np.array([-3/2, -1/2, 1/2, 3/2]) / 2 * x_lambda6
y_cross6 = np.zeros_like(x_cross6)

# %% Setting figure properties
LineWidth = 3
FontSize = 16
QuivLength = 2
HeadWidth = 5
PointSize = 20

x_min, x_max = -10, 10
z_min, z_max = -6, 6

def get_norm(x_cross, x_lambda):
    norm = np.sqrt((6 * np.pi / x_lambda * np.sin(2 * np.pi / x_lambda * x_cross))**2 + 1)
    x_norm = (6 * np.pi / x_lambda * np.sin(2 * np.pi / x_lambda * x_cross)) / norm
    y_norm = 1 / norm
    return x_norm, y_norm

def plot_HCS(ax, x, y, alpha=1):
    cmap = plt.cm.winter
    norm = Normalize(vmin=x.min(), vmax=x.max())
    ax.scatter(x, y, c=cmap(norm(x)), s=PointSize, alpha=alpha, edgecolors='none')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(z_min, z_max)
    ax.set_aspect('equal')

def plot_HCS_normal(ax, x_cross, y_cross, x_lambda, quiv_length, head_width, line_width):
    for i in range(4):
        x_sub = x_cross[i]
        y_sub = y_cross[i]
        x_norm, y_norm = get_norm(x_sub, x_lambda)
        u = quiv_length * x_norm
        w = quiv_length * y_norm
        ax.quiver(x_sub, y_sub, u, w, color='red', linewidth=line_width,
                  headwidth=head_width, angles='xy', scale_units='xy', scale=1)

def plot_PSP_obiter(ax):
    ax.plot(x_PSP, y_PSP, 'k', linewidth=LineWidth)
    ax.quiver(orb_PSP_beg, 0, orb_PSP_quiv, 0, color='black', linewidth=LineWidth,
                 headwidth=HeadWidth, angles='xy', scale_units='xy', scale=1)
    
def plot_steady_motion(ax, x, y_up, y_dw):
    ax.quiver(x, y_up, 0,  QuivLength, edgecolor='orange', facecolor='none', 
            linewidth=LineWidth, headwidth=HeadWidth*2, angles='xy', scale_units='xy', scale=1)
    ax.quiver(x, y_dw, 0, -QuivLength, edgecolor='orange', facecolor='none', 
            linewidth=LineWidth, headwidth=HeadWidth*2, angles='xy', scale_units='xy', scale=1)

def plot_kinklike_motion(ax, x_up, y_up, x_dw, y_dw):
    ax.quiver(x_up, y_up, QuivLength, 0, edgecolor='orange', facecolor='none', 
            linewidth=LineWidth, headwidth=HeadWidth*2, angles='xy', scale_units='xy', scale=1)
    ax.quiver(x_dw, y_dw, QuivLength, 0, edgecolor='orange', facecolor='none', 
            linewidth=LineWidth, headwidth=HeadWidth*2, angles='xy', scale_units='xy', scale=1)

# %% Plotting figure
fig, axes = plt.subplots(2, 2, figsize=(10, 6))
fig.subplots_adjust(hspace=0.3, wspace=0.1)

# Situation-4
axes[0,0].set_title('(a) Static folded CS', fontsize=FontSize)
plot_HCS(axes[0,0], x_HCS, y_HCS4)
plot_HCS_normal(axes[0,0], x_cross4, y_cross4, x_lambda4, QuivLength, HeadWidth, LineWidth)
plot_PSP_obiter(axes[0,0])
axes[0,0].axis('off')

# Situation-2
axes[0,1].set_title('(b) Steady flapping of unfolded CS', fontsize=FontSize)
plot_HCS(axes[0,1], x_HCS, y_HCS2)
plot_HCS(axes[0,1], x_HCS, y_HCS2+3, alpha=1e-2)
plot_HCS(axes[0,1], x_HCS, y_HCS2-3, alpha=5e-3)
plot_PSP_obiter(axes[0,1])
plot_HCS_normal(axes[0,1], x_cross2, y_cross2, 1e10, QuivLength, HeadWidth, LineWidth)
plot_steady_motion(axes[0,1], x=-8, y_up=1, y_dw=-1)
axes[0,1].text(-9.8, 0.4, r'$t_1$', color='darkblue', fontsize=FontSize, fontweight='bold')
axes[0,1].text(-9.8, 3.3, r'$t_2$', color='darkblue', fontsize=FontSize, fontweight='bold')
axes[0,1].text(-9.8,-2.7, r'$t_3$', color='darkblue', fontsize=FontSize, fontweight='bold')
axes[0,1].axis('off')

# Situation-6
axes[1,0].set_title('(c) Kink-like flapping of folded CS', fontsize=FontSize)
plot_HCS(axes[1,0], x_HCS, y_HCS6)
plot_HCS(axes[1,0], x_HCS+3, y_HCS6, alpha=1e-2)
plot_HCS(axes[1,0], x_HCS+6, y_HCS6, alpha=5e-3)
plot_HCS_normal(axes[1,0], x_cross6, y_cross6, x_lambda6, QuivLength, HeadWidth, LineWidth)
plot_PSP_obiter(axes[1,0])
plot_kinklike_motion(axes[1,0], x_up=-7.5, y_up=2, x_dw=-9.5, y_dw=-2)
axes[1,0].text(-9.8, 3.3, r'$t_1$', color='darkblue', fontsize=FontSize, fontweight='bold')
axes[1,0].text(-6.8, 3.3, r'$t_2$', color='darkblue', fontsize=FontSize, fontweight='bold')
axes[1,0].text(-3.8, 3.3, r'$t_3$', color='darkblue', fontsize=FontSize, fontweight='bold')
axes[1,0].axis('off')

# Situation-5
axes[1,1].set_title('(d) Steady flapping of folded CS', fontsize=FontSize)
plot_HCS(axes[1,1], x_HCS, y_HCS5)
plot_HCS(axes[1,1], x_HCS, y_HCS5+1.8, alpha=1e-2)
plot_HCS(axes[1,1], x_HCS, y_HCS5-1.8, alpha=5e-3)
plot_HCS_normal(axes[1,1], x_cross5, y_cross5, x_lambda5, QuivLength, HeadWidth, LineWidth)
plot_PSP_obiter(axes[1,1])
plot_steady_motion(axes[1,1], x=-8, y_up=2.2, y_dw=0.2)
axes[1,1].text(-9.7, 3.1, r'$t_1$', color='darkblue', fontsize=FontSize, fontweight='bold')
axes[1,1].text(-9.7, 4.9, r'$t_2$', color='darkblue', fontsize=FontSize, fontweight='bold')
axes[1,1].text(-9.7, 1.3, r'$t_3$', color='darkblue', fontsize=FontSize, fontweight='bold')
axes[1,1].axis('off')

# Adding legend
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=plt.cm.winter(0.5), 
                      markersize=PointSize//2, linewidth=0),
           plt.quiver([0], [0], [1], [0], color='k', linewidth=0.5, headwidth=0.5),
           plt.quiver([0], [0], [1], [0], color='r', linewidth=0.5, headwidth=0.5)]
labels = ['Current Sheet', 'Orbit of Satellite', 'CS Normal Direction']
fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.02, 0.02),
           fontsize=FontSize-2, ncol=1)

# Adding stamp
stamp = 'plotted by plot_situations_for_successive_crossings_revised.py'
fig.text(0.99, 0.01, stamp, ha='right', fontsize=FontSize-4)

plt.show()
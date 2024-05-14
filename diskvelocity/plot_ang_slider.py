"""
show the polarization direction as a function of St for different values of alpha and beta
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
rad = np.pi / 180

# The parametrized function to be plotted
def f(St, alpha, beta):
    return np.arctan(0.5 * (-1.5*alpha + beta * St) / (-(1.5*alpha*St + beta))) / rad

St = np.geomspace(1e-3, 1e3, 200)

# Define initial parameters
init_alpha = 0.01
init_beta = 1.5

# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()
line, = ax.plot(St, f(St, init_alpha, init_beta), lw=2)
ax.set_xlabel('St')
ax.set_xscale('log')
ax.set_title('polarization angle')
ax.set_yticks(np.arange(-90, 91, 30))
ax.set_ylim(-95, 95)
for ival, txt, va in zip([-90, 0, 90], 
                    ['azimuthal', 'radial', 'azimuthal'],
                    ['bottom', 'bottom', 'top']):
    ax.axhline(ival, color='k', linestyle=':')
    ax.text(St[-1], ival, txt, color='k', va=va, ha='right')

# adjust the main plot to make room for the sliders
fig.subplots_adjust(left=0.25, bottom=0.25)

# Make a horizontal slider to control beta
axbeta = fig.add_axes([0.25, 0.1, 0.65, 0.03])
beta_slider = Slider(
    ax=axbeta,
    label='beta',
    valmin=-1,
    valmax=1,
    valinit=init_beta,
)

# Make a vertically oriented slider to control alpha
axalpha = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
alpha_slider = Slider(
    ax=axalpha,
    label="alpha",
    valmin=0,
    valmax=0.2,
    valinit=init_alpha,
    orientation="vertical"
)


# The function to be called anytime a slider's value changes
def update(val):
    line.set_ydata(f(St, alpha_slider.val, beta_slider.val))
    fig.canvas.draw_idle()


# register the update function with each slider
beta_slider.on_changed(update)
alpha_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    beta_slider.reset()
    alpha_slider.reset()
button.on_clicked(reset)

plt.show()

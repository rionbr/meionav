
#
# To include Helvetica as a font you need to place the .ttf files inside /u/<user>/.local/share/fonts
#
#import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'

#'font.sans-serif': ['DejaVu Sans', 'Bitstream Vera Sans','Computer Modern Sans Serif', 'Lucida Grande', 'Verdana', 'Geneva', 'Lucid', 'Arial', 'Helvetica', 'Avant Garde', 'sans-serif'],
#mpl.rcParams['text.usetex'] = False
#mpl.rcParams['figure.titlesize'] = 'medium'
#mpl.rcParams['axes.titlesize'] = 'small'
#mpl.rcParams['axes.labelsize'] = 'small'
#mpl.rcParams['xtick.labelsize'] = 'x-small'
#mpl.rcParams['ytick.labelsize'] = 'x-small'
#mpl.rcParams['legend.fontsize'] = 'small'
mpl.rcParams['hatch.linewidth'] = 0.5
mpl.rcParams['hatch.color'] = '#969696'

fig, ax = plt.subplots(figsize=(4, 4), nrows=1, ncols=1)
ax.plot([-1,0,1],[-1,0,1])
ax.set_title('This is Helvetica')

ax.text(-0.8,0.1, "Helvetica")
ax.text(-0.8,0.2, "Helvetica Italic", style="italic")
ax.text(-0.8,0.3, "Helvetica Bold", weight="bold")
ax.text(-0.8,0.4, "Helvetica Bold Italic", weight="bold", style="italic")
ax.text(-0.8,0.5, "Helvetica Light", weight="light")
ax.text(-0.8,0.6, "Helvetica Light Italic", weight="light", style="italic")
ax.text(-0.8,0.7, "Helvetica Light", weight="light")
ax.text(-0.2,0.-0.6, r"$ 1 = \frac{-123}{+456} \alpha \beta \Delta \sigma \sum^{\infty}_{n=1} \int_a^b$", fontsize='large')
ax.set_ylabel("Computer Modern is the default Latex Font.")
ax.set_xlabel(r"Math notation is Computer Modern")

plt.tight_layout()
plt.savefig('images/font-test.pdf')
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# Create a range of values for the x-axis (all positive)
x = np.linspace(0, 8, 100)

# Calculate the probability density function (PDF) for the shifted Gaussian distribution
mu = 4  # Shifted mean
sigma = 1  # Standard deviation
pdf = norm.pdf(x, mu, sigma)

data = np.random.normal(mu, sigma, 1000)
# Create a histogram to visualize the distribution
# plt.hist(data, bins=20, density=True, alpha=0.6, color='darkorange')

# Plot the PDF as a line
plt.plot(x, pdf, 'dimgrey', lw=2)

# Fill the area under the PDF curve between 0.8 and 0.9
x_fill = np.linspace(5, 5.5, 100)  # Shifted range
pdf_fill = norm.pdf(x_fill, mu, sigma)
plt.fill_between(x_fill, pdf_fill, alpha=0.5, color='firebrick')

plt.xlabel('Score', fontsize=14)
plt.ylabel('Frequency in population', fontsize=14)
plt.xticks([])
plt.yticks([])

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)

# plt.show()
plt.savefig('distro-part.pdf')

import os
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
from scipy.stats import kruskal, t


# noinspection PyShadowingNames
class HmDist:
	def __init__(self,
	             list_rgb_color,
	             yyyymmdd,
	             str_working_dir,
	             deg,
	             list_hm,
	             file_postfix,
	             ):
		"""
		Initialize variables
			list_rgb_color:
			yyyymmdd:
			str_working_dir:
			deg: Polynomial regression.
			list_hm: A list of heavy metals
			file_prefix: File postfix
		"""
		# Str
		self.str_working_directory = str_working_dir
		self.str_output_directory = self.str_working_directory + "3.2_HMs_Distribution\\"
		self.str_excel_path = self.str_working_directory +f"2.2_Imputation\\Imputation{yyyymmdd}\\03Revised\\Preprocessed{yyyymmdd}.xlsx"
		self.str_world_path = self.str_working_directory + "Worldmap\\Worldmap.shp"
		self.str_continent_path = self.str_working_directory+"Continent\\Continent\\continent.shp"
		# Dataframe
		self.df = None
		self.df_temp = None
		self.df_parts = None
		self.df_latitude = None
		# Num
		self.num_degree = deg
		self.num_symmetric_axis = 0
		# List
		self.list_hm = list_hm
		self.list_rgb_color = list_rgb_color
		self.list_ax_legend = []
		self.list_ax_world = []
		self.list_ax_scatter = []
		self.list_head = []
		self.list_colors = []
		# Color schemes
		self.linear_cmap = None
		self.palette = None
		# Initialize
		list_path = [self.str_working_directory,
		             self.str_output_directory,
		             self.str_excel_path,
		             self.str_world_path,
		             self.str_continent_path]
		for p in list_path:
			if not os.path.exists(p) and p[-1] == "\\":
				os.mkdir(p)
			elif not os.path.exists(p) and p[-1] != "\\":
				raise ValueError("Invalid path: " + p)

		self.color_setup()
		self.built_canvas()
		self.read_excel(delete_abnormal=True)

		for n, (HM, color) in enumerate(zip(self.list_hm, self.list_colors)):
			self.draw_scatter(self.list_ax_scatter[n], HM)
			self.draw_world(self.list_ax_world[n], HM, mode="scatter")
		self.add_serials(list_ax = self.list_ax_scatter + self.list_ax_world)

		# plt.show()
		plt.savefig(self.str_output_directory+f"HM_{file_postfix}_{yyyymmdd}.png",dpi=300)
		plt.savefig(self.str_output_directory + f"HM_{file_postfix}_{yyyymmdd}.svg")
		plt.close()

	@staticmethod
	def rgb2hex(rgb):
		"""
		Convert the color in RGB format to hexadecimal string

		___________
		Parameters:
			rgb: A set with three elements, indicating "R", "G" and "B" values

		_______
		Return:
			color: A string-like hexadecimal color
		"""

		color = "#"
		for i in rgb:  # Iterate each color
			num = int(i)  # Get integer in case of any unsupported format
			color += str(hex(num))[-2:].replace("x", "0").upper()  # Convert
		return color

	@staticmethod
	def delete_abnormal_values(df, column_drop):
		"""
		Remove the abnormal values
		Fill null values with the mean values

		___________
		Parameters:
			df: A dataframe with potential abnormal values
			column_drop: Do not delete the abnormal values in the column_drop
		_______
		Return:
			df: A dataframe without abnormal values
		"""

		# Load the Head of the DataFrame. Create a New List to Update the Head
		column_head = list(df)  # Get a list of the column names
		list_head = [i for i in column_head if i not in column_drop]
		df_copy = df.copy() # Copy the dataframe

		# Select the Heavy Metal Columns
		for num, h in enumerate(list_head): # Iterate each column name
			q1 = df_copy[h].quantile(0.25)  # Get values at the 25% quantile point
			q3 = df_copy[h].quantile(0.75)  # Get values at the 75% quantile point
			iqr = q3 - q1   # Get the data value range within 25% - 75%
			# Determine the Range to Distinguish the Abnormal Values
			num_lower_limit = q1 - 1.5 * iqr    # Upper limit
			num_upper_limit = q3 + 1.5 * iqr    # Lower limit
			num_mean = np.mean(df_copy[h])  # Calculate the Mean Value of the Column
			# Replace the Abnormal Values with the Mean Value of Corresponding Columns
			df_copy.loc[(df_copy[h] > num_upper_limit) | (df_copy[h] < num_lower_limit), h] = num_mean

		print("Abnormal data are removed!")
		return df_copy

	@staticmethod
	def add_serials(list_ax, position=None):
		"""
		Add serial marks for each subplots

		___________
		Parameters:
			list_ax: A list of axes that need serial marks
			position: A list with only two values within range [0,1], used to adjust the serial mark positions
		"""

		if position is None:
			position = [0.95, 0.05]
		# List needs to be added with serial
		list_tb = [chr(i) for i in range(97, 123)]  # A list of alphabetic serials
		for n, la in enumerate(list_ax):  # Iterate each subplot
			# Add serial marks
			la.text(x=la.get_xlim()[0] * (1 + position[1]) - la.get_xlim()[1] * position[1],  # X coordinates
			        y=la.get_ylim()[1],  # Y coordinates
			        s=list_tb[n],  # Text content
			        weight="bold",  # Bold text
			        ha="left",  # Horizontal alignment
			        va="top",  # Vertical alignment
					fontsize= 16)   # Font size
		print("Serials have been added")

	def color_setup(self):
		"""
		Setup palettes
		"""

		# Convert RGB colors to HEX ones
		color_indict = [self.rgb2hex(color) for color in self.list_rgb_color]
		# Create a linear color map according to the given color list
		self.linear_cmap = LinearSegmentedColormap.from_list("custome", color_indict)
		# Create a list of even discrete colors, used for HMs
		self.list_colors = [self.linear_cmap(i) for i in np.linspace(0, 1, len(self.list_hm))]
		# Create a palette, used when we used "seaborns" to plot
		self.palette = sns.color_palette(self.list_colors)

	def built_canvas(self, fig_width=12, h_space=0.3, w_space=0.1):
		"""
		Build up a canvas

		___________
		Parameters:
			fig_width: The width of the figure
			h_space: The vertical gap between two subplots
			w_space: The horizontal gap between two subplots
		"""
		fig_height = 3 * len(self.list_hm)  # The height of the figure
		fig = plt.figure(figsize=(fig_width, fig_height))  # Figure object
		# fig.patch.set_alpha(0.0)
		gs = GridSpec(nrows = len(self.list_hm),
		              ncols = 3,
		              height_ratios=[1]* len(self.list_hm),
		              width_ratios=[1, 1, 0.02],
					  hspace=h_space,
					  wspace=w_space)
		for i in range(len(self.list_hm)):
			self.list_ax_scatter.append(fig.add_subplot(gs[i, 0]))  # Create a list of axes to plot the scatters
			self.list_ax_world.append(fig.add_subplot(gs[i, 1]))    # Create a list of axes to plot the world map
			self.list_ax_legend.append(fig.add_subplot(gs[i, 2]))   # Create a list of axes to plot the world map legend

	def read_excel(self, delete_abnormal = True):
		"""
		Read a dataframe

		___________
		Parameters:
			delete_abnormal: Bool
		"""

		# Read a dataframe
		self.df = pd.read_excel(self.str_excel_path, header=[0])
		# Extract rows, indicating Anurans with Gosner Stage larger than 46 (that is, metamorphosed Anurans)
		self.df = self.df[self.df["Concentration-Gosner_Stage"] > 46]
		# Extract whole-body HM bio-accumulation columns
		self.list_head = [x for x in list(self.df) if "Whole_Body" in x and x.split("-")[-1] in self.list_hm]
		# Extract columns
		self.df_temp = self.df[["Longitude", "Latitude", "Country", "Continent"] + self.list_head].copy()
		# Rename the columns
		self.df_temp.columns = [x.replace("log_Concentration-Whole_Body-", "")  for x in list(self.df_temp) ]
		# Rebuild a dataframe
		self.df_parts = self.df_temp.melt(id_vars=["Continent", "Country"])
		# Drop the column "Country"
		self.df_temp.drop(["Country","Continent"], axis=1, inplace=True)
		# Delete abnormal values
		if delete_abnormal: # If needed
			self.df_temp = self.delete_abnormal_values(self.df_temp,
			                                      column_drop=["Longitude", "Latitude"])

		print("The dataframe has been read")

	def draw_world(self, ax, hm, mode="scatter"):
		"""
		Draw a world map with scatters

		___________
		Parameters:
			ax: Axes
			hm: Heavy metals. Choose from "SVL", "BM" and "CF"
			mode: The way to visualize the data. Choose between "scatter" and "country"
		"""

		world = gpd.read_file(self.str_world_path)  # Read a shape file
		world.columns = ['NAME_CHN', 'Country', 'NR_C', 'NR_C_ID', 'SOC', 'geometry']   # Rename the columns
		world.drop(world[world["Country"] == "ANTARCTICA"].index, axis=0, inplace=True) # Drop the row indicating Antarctica
		list_country = list(world["Country"])   # Get the country list from the column "Country"

		if mode == "country":   # If we need to visualize the data across countries
			list_C = [] # Create an empty list
			for i in list_country:  # Iterate all countries from the list
				list_C.append(i.upper())    # Append all country
			world["Country"] = list_C   # Update the column "Country"

			# World distribution of the studied contaminants
			self.df_parts = self.df_parts[self.df_parts['Contaminant'].isin(self.list_hm)]
			list_country = list(set(self.df_parts["Country"]))  # Create a non-duplicated list of countries

			list_continent_country = [[self.df_parts.loc[i, "Country"], self.df_parts.loc[i, "Continent"]] for i in
			                          list(self.df_parts.index)]
			dict_continent_country = {}
			previous_item = []

			for current_item in list_continent_country:
				if current_item != previous_item:
					dict_continent_country[current_item[0]] = current_item[1]
					previous_item = current_item
			world['Continent'] = world['Country'].map(dict_continent_country)   # Create a column "Continent"

			for lcc in list_country:    # Iterate all countries
				num_m = self.df_parts[(self.df_parts["Country"] == lcc) & (self.df_parts["Contaminant"] == hm)][
					"Concentration"].mean()
				world.loc[world["Country"] == lcc.upper(), "Concentration"] = num_m

			mean_min = min(world["Concentration"].dropna(how="any", axis=0))    # Get the minimum concentration
			mean_max = max(world["Concentration"].dropna(how="any", axis=0))    # Get the maximum concentration
			distance_min_max = abs(mean_min - mean_max)                         # Get the difference between max and min

			world.plot(ax=ax,                                                   # Plot a world map on "ax"
			           column="Concentration",                                  # Display values
			           cmap=self.linear_cmap,                                   # Color map
			           edgecolor='gray',                                        # Edge color
			           linewidth=0.05,                                          # Line width
			           legend=False,                                            # Do not show legend
			           missing_kwds={"color": "#ededed",                        # Set a specific face color
			                         "edgecolor": "gray",                       # Set a specific edge color
			                         "hatch": ""})                              # Set a fill pattern
			ax.axhline(y = 0,                                                   # Y coordinates
			           linestyle="--",                                          # Line style
			           linewidth=0.5,                                           # Line width
			           color="gray",                                            # Line color
			           zorder=0)                                                # Set as the bottommost layer
			ax.axhline(y = 30,                                                  # Y coordinates
			           linestyle="--",                                          # Line style
			           linewidth=0.2,                                           # Line width
			           color="gray",                                            # Line color
			           zorder=0)                                                # Set as the bottommost layer
			ax.axhline(y = -30,                                                 # Y coordinates
			           linestyle="--",                                          # Line style
			           linewidth=0.2,                                           # Line width
			           color="gray",                                            # Line color
			           zorder=0)                                                # Set as the bottommost layer

			n = self.list_hm.index(hm)                                          # Get the current order
			self.list_ax_legend[n].text(x = 1,                                  # X coordinates
			                            y = mean_min,                           # Y coordinates
			                            s = f"{mean_min:.2f}",                  # Text
			                            ha='left',                              # Horizontal alignment
			                            va='bottom')                            # Vertical alignment
			self.list_ax_legend[n].text(x = 1,                                  # X coordinates
			                            y = mean_max,                           # Y coordinates
			                            s = f"{mean_max:.2f}",                  # Text
			                            ha='left',                              # Horizontal alignment
			                            va='top')                               # Vertical alignment
			self.list_ax_legend[n].text(x = 1.5,                                # X coordinates
			                            y = (mean_min + mean_max) / 2,          # Y coordinates
			                            s = f"lg {hm}",                         # Text
			                            ha='left',                              # Horizontal alignment
			                            va='center',                            # Vertical alignment
			                            rotation=90)                            # Rotation
			# Adjust axes and legends
			ax.set_ylim(-70, 90)                                                # Set y boundary
			ax.axis("off")                                                      # Turn off all axes
			gradient = np.linspace(1, 0, 256).reshape(-1, 1)     # Build up a color list
			gradient = np.hstack((gradient, gradient))                          # Convert to a 2D array
			self.list_ax_legend[n].imshow(gradient,                             # Data
			                              cmap=self.linear_cmap,                # Color map
			                              aspect='auto',                        # Image aspect
			                              extent=[0, 1, mean_min, mean_max])    # [x_start, y_start, dx, dy]
			self.list_ax_legend[n].set_ylim(mean_min - 0.1 * distance_min_max,  # Set lower y boundary
			                                mean_max + 0.1 * distance_min_max)  # Set upper y boundary
			self.list_ax_legend[n].axis("off")                                  # Turn off all axes
			self.list_ax_legend[n].set_xticks([])                               # Clear ticks on X axis
			self.list_ax_legend[n].set_yticks([])                               # Clear ticks on Y axis

		elif mode == "scatter":
			scatter_data = self.df_temp[["Longitude", "Latitude", hm]].copy()  # Extract columns
			# Calculate Statistic parameters for scatter size on the world map
			min_val = scatter_data[hm].min()    # Get a minimum number of the column "lc"
			max_val = scatter_data[hm].max()    # Get a maximum number of the column "lc"
			scatter_data[hm + "_scaled"] = (scatter_data[hm] - min_val) / (max_val - min_val)  # Standardize

			# Build a dataframe to draw scatters on a world map
			scatter_gdf = gpd.GeoDataFrame(scatter_data,    # Original data
			                               geometry=gpd.points_from_xy(scatter_data['Longitude'],   # Set X coordinates
			                                                           scatter_data['Latitude']))   # Set Y coordinates
			# Draw world maps
			world.plot(ax=ax,                                                  # Axes
			           color="#f5f5f5",                                        # Face color
			           edgecolor='gray',                                       # Edge color
			           linewidth=0.1)                                          # Edge width
			# Draw scatters on the world map
			scatter_map = ax.scatter(x = scatter_gdf.geometry.x,               # Longitude
			                         y = scatter_gdf.geometry.y,               # Latitude
			                         c = scatter_gdf[hm],          # Color changed with values
			                         cmap = self.linear_cmap,  # Color map
			                         s = scatter_gdf[hm + "_scaled"] * 50,     # Scatter size
			                         alpha = 0.4)                              # Alpha
			# Draw color bars
			cbar = plt.colorbar(scatter_map,
			                    cax = self.list_ax_legend[self.list_hm.index(hm)])
			cbar.set_label(f"lg {hm}")                                         # Set legend title

			# Add lines
			ax.axhline(y = 0,                                                  # Y coordinates
                       linestyle="--",                                         # Line style
                       linewidth=0.5,                                          # Line width
                       color="gray",                                           # Color
                       zorder=0)                                               # Set as the bottommost layer
			ax.axhline(y = 30,                                                 # Y coordinates
	                   linestyle="--",                                         # Line style
	                   linewidth=0.2,                                          # Line width
	                   color="gray",                                           # Color
	                   zorder=0)                                               # Set as the bottommost layer
			ax.axhline(y = -30,                                                # Y coordinates
	                   linestyle="--",                                         # Line style
	                   linewidth=0.2,                                          # Line width
	                   color="gray",                                           # Color
	                   zorder=0)                                               # Set as the bottommost layer
			# Set axes
			ax.set_ylim(-70, 90)                                               # Set the limit of Y axis
			ax.axis("off")                                                     # Turn off all the axes
			ax.axhline(self.num_symmetric_axis,                                    # X coordinates
	                   linestyle="--",                                         # Line style
	                   linewidth=0.3,                                          # Line width
	                   color="red")                                            # Line color

	def draw_scatter(self, ax, hm):
		"""
		Draw quantile regression scatters

		___________
		Parameters:
			ax: Axes
			hm: Heavy metals. Choose from "SVL", "BM" and "CF"
		"""
		n = self.list_hm.index(hm)
		self.df_latitude = self.df_temp[["Latitude", hm]].copy()  # Extract columns
		list_lat = self.df_latitude["Latitude"] # Get the column "Latitude"


		color = self.list_colors[n]
		list_hm = self.df_latitude[hm]  # Get the column "hm"

		# Quadratic polynomial regression
		coefficients = np.polyfit(x = list_lat,         # Set "Latitude" as the independent values
		                          y = list_hm,          # Set "hm" as the independent values
		                          deg = self.num_degree)    # Polynomial degrees. Here we set "2"
		poly = np.poly1d(coefficients)  # Sort parameters

		# Predict y
		XX = np.linspace(min(list_lat), max(list_lat), 100) # Create a list of X evenly
		y_pred = poly(XX)   # Get predicted Y values

		# Calculate standard deviation of residuals
		num = len(list_lat) # Get the row number
		list_hat = poly(list_lat)   # Predict with "poly"
		residuals = list_hm - list_hat  # Residuals
		se = np.sqrt(np.sum(residuals ** 2) / (num - self.num_degree - 1)) # Standard deviation of residuals

		# Calculate independent matrices and a covariance matrix
		X_design = np.column_stack([list_lat ** i for i in range(self.num_degree + 1)]) # lat^2 + lat^1 + lat^0
		X_design_XX = np.column_stack([XX ** i for i in range(self.num_degree + 1)])    # XX^2 + XX^1 + XX^0
		cov_matrix = se ** 2 * np.linalg.inv(X_design.T @ X_design) # Covariance matrix

		# Standard Error
		y_pred_se = np.sqrt(np.sum(X_design_XX @ cov_matrix * X_design_XX, axis=1))

		# Calculate confidence interval
		alpha = 0.05  # 95% Confidence Interval
		t_value = t.ppf(1 - alpha / 2, num - self.num_degree - 1)
		ci = t_value * y_pred_se    # Get confidence interval

		# R^2
		y_mean = np.mean(list_hm)   # Calculate the mean concentration of "hm"
		ss_res = np.sum(residuals ** 2)
		ss_tot = np.sum((list_hm - y_mean) ** 2)
		r_squared = 1 - (ss_res / ss_tot)   # Get R^2

		# Symmetric axis
		a = coefficients[0] # The coefficient of the quadratic term
		b = coefficients[1] # The coefficient of the linear term
		self.num_symmetric_axis = -b / 2 / a                            # Get the symmetric axis

		# Annotation
		para_text = ""                                                  # Create an empty string
		list_para = [chr(i) for i in range(97, 123)]                    # Create a list of alphabet

		for i, coeff in enumerate(coefficients):                        # Iterate all coefficients
			# Decompose the scientific notation aEb into a and b
			num = str(f"{coeff:.2e}").split("e")                        # Suppose "coeff" is 3.456e-2, "num" will be ["3.46", "-2"]
			new_num = num[0] if coeff < 0 else " " + num[0]
			num[1] = int(num[1])                                        # The exponents of the scientific notation is in num[1], update in case of float format
			if num[1] != 0:                                             # If the exponent is not 0
				new_num += "×10"
				for content in str(num[1]):
					new_num += f"$^{content}$"
			para_text += f"{list_para[i]} = {new_num}\n"                # Append coefficients

		# Plot
		ax.set_ylim(-70, 99)                                            # Set Y axis range
		ax.set_xlim(min(list_hm), max(list_hm))                         # Set X axis range
		ax.axhline(y = 0,                                               # Y coordinates
		           linestyle="--",                                      # Line style
		           linewidth=0.5,                                       # Line width
		           color="gray")                                        # Line color
		ax.axhline(y = -30,                                             # Y coordinates
		           linestyle="--",                                      # Line style
		           linewidth=0.2,                                       # Line width
		           color="gray")                                        # Line color
		ax.axhline(y = 30,                                              # Y coordinates
		           linestyle="--",                                      # Line style
		           linewidth=0.2,                                       # Line width
		           color="gray")                                        # Line color
		ax.axhline(y = self.num_symmetric_axis,                         # Y coordinates
		           linestyle="--",                                      # Line style
		           linewidth=0.3,                                       # Line width
		           color="red")                                         # Line color
		ax.vlines(x = ax.get_xlim()[1] * 0.99 + ax.get_xlim()[0] * 0.01,# X coordinates
		          ymin = -50,                                           # Minimum Y coordinates
		          ymax = 50,                                            # Maximum Y coordinates
		          linestyle="-",                                        # Line style
		          linewidth=0.5,                                        # Line width
		          color="black")                                        # Line color
		ax.fill_betweenx(XX,                                            # X coordinates
		                 y_pred - ci,                                   # Lower Y boundary
		                 y_pred + ci,                                   # Upper Y boundary
		                 color='gray',                                  # Color
		                 alpha=0.05,                                    # Alpha
		                 linewidth=0)                                   # Line width
		ax.scatter(x = list_hm,                                         # X coordinates
		           y = list_lat,                                        # Y coordinates
		           label=hm,                                            # Label
		           alpha=0.2,                                           # Alpha
		           color=color)                                         # Color
		ax.plot(y_pred,                                                 # X coordinates "Latitude"
		        XX,                                                     # Y coordinates "hm"
		        color= "red",                                           # Color
		        lw = 1)                                                 # Line width
		ax.text(x = ax.get_xlim()[1] * 0.05 + ax.get_xlim()[0] * 0.95,  # X coordinates
                y = ax.get_ylim()[1],                                   # Y coordinates
                s = f"{para_text}R² = {r_squared:.3f}",                 # Text
                fontsize=fontsize,                                      # Font size
                va='top',                                               # Vertical alignment
                ha="left")                                              # Horizontal alignment
		ax.text(x = ax.get_xlim()[0] * 0.05 + ax.get_xlim()[1] * 0.95,  # X coordinates
                y = self.num_symmetric_axis,                            # Y coordinates
                s = f"{self.num_symmetric_axis:.3f}°N",                 # Text
                fontsize=fontsize,                                      # Font size
                color="red",                                            # Color
                va='top',                                               # Vertical alignment
                ha="right")                                             # Horizontal alignment
		# Adjust axes
		ax.spines["right"].set_visible(False)                           # Hide the right axis
		ax.spines["left"].set_position(("data", ax.get_xlim()[1]))      # Move the left axis to the right
		ax.spines["left"].set_visible(False)                            # Hide the left axis
		ax.spines["top"].set_visible(False)                             # Hide the top axis
		ax.yaxis.set_major_locator(MultipleLocator(50))                 # Y ticks interval
		ax.text(x = ax.get_xlim()[1],                                   # X coordinates
                y = 0,                                                  # Y coordinates
                s = "Latitude",                                         # Text
		        ha="left",                                              # Horizontal alignment
		        va="center",                                            # Vertical alignment
		        rotation=90,                                            # Y axis label
		        fontsize=fontsize)                                      # Font size
		ax.set_xlabel(f"lg {hm}",                                       # X axis label
                      fontsize=fontsize)                                # Font size

		print("Successful scatter plot")

if __name__ == "__main__":
	# Color schemes
	list_color = [[14, 91, 118],
	              [26, 134, 163],
	              [70, 172, 202],
	              [155, 207, 232],
	              [205, 205, 164],
	              [255, 202, 95],
	              [254, 168, 9],
	              [253, 152, 2],
	              [251, 132, 2]]

	degree = 2  # Quadratic regression
	fontsize = 10   # Fontsize
	fontname = "Arial"
	yyyymmdd = 20250716
	str_working_directory = "E:\\学习\\研二上\\重金属\\"

	# Main text figures
	HmDist(list_rgb_color=list_color,
	        yyyymmdd = yyyymmdd,
	        str_working_dir=str_working_directory,
	        deg=degree,
	        list_hm=["Cr", "Mn"],
	        file_postfix = "main_text")
	# Supplementary figures
	HmDist(list_rgb_color=list_color,
	        yyyymmdd=yyyymmdd,
	        str_working_dir=str_working_directory,
	        deg=degree,
	        list_hm=["Cd", "Cu", "Fe", "Pb", "Zn"],
	        file_postfix = "Supplementary_figure")


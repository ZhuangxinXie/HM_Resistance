import copy
import math
import os
from datetime import datetime
import chordplot
import geopandas as gpd
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import skewness_kurtosis as sk
import svgpath2mpl
from Bio import Phylo
from PIL import Image
from docutils.nodes import legend
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
from matplotlib.offsetbox import AnnotationBbox, DrawingArea, OffsetImage, TextArea
from matplotlib.patches import PathPatch
from mpl_chord_diagram import chord_diagram
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from openpyxl.styles.builtins import output
from scipy.stats import kruskal, linregress, t
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
from scipy import stats
from warnings import simplefilter

simplefilter(action='ignore', category=UserWarning)

class Body_size_indicators:
	def __init__(self,
	             list_color,
	             str_df_path,
	             working_directory,
	             tree_path_all,
	             tree_scale_factor,
	             width_fraction,
	             indicator_type
	             ):
		self.list_color = list_color
		self.linearcmap = LinearSegmentedColormap.from_list('custome',[self.RGB_to_Hex(color) for color in list_color])
		self.indicator_type = indicator_type
		self.working_directory = working_directory
		self.output_directory = self.working_directory + "3.1_Body_Features\\"
		self.df_data_path = str_df_path
		self.shp_path = self.working_directory + "Worldmap\\Worldmap.shp"
		self.illustration_path = self.working_directory + "3.1_Body_Features\\SVG\\"
		self.tree_path_all = tree_path_all
		self.tree_scale_factor = tree_scale_factor
		self.width_fraction = width_fraction
		self.ax_chord = None
		self.ax_scatter = None
		self.ax_scatter_itv = None
		self.ax_map = None
		self.ax_map_legend = None
		self.ax_tree = None
		self.ax_heatmap = None
		self.ax_legend = None
		self.custom_cmap = LinearSegmentedColormap.from_list('custome', self.list_color)
		self.canvas_setup()
		if indicator_type == "BS":
			self.plot_chord()
		else:
			pass
		if self.indicator_type == "BS":
			plt.close()
			self.canvas_setup()
			self.plot_world()
			self.plot_scatter()
			self.add_serial(list_axes=self.ax_scatter + self.ax_map)
			plt.savefig(f"{self.output_directory}BS_distribution_{today.strftime("%Y%m%d")}.png",dpi=300)
			plt.savefig(f"{self.output_directory}BS_distribution_{today.strftime("%Y%m%d")}.svg")
			plt.close()
		elif self.indicator_type == "ITV":
			plt.close()
			self.canvas_setup()
			self.plot_world_itv()
			self.plot_scatter_itv()
			self.add_serial(list_axes=self.ax_scatter + self.ax_map)
			plt.savefig(f"{self.output_directory}ITV_distribution_{today.strftime("%Y%m%d")}.png", dpi=300)
			plt.savefig(f"{self.output_directory}ITV_distribution_{today.strftime("%Y%m%d")}.svg")
			plt.close()

		self.plot_tree()

		plt.show()

	def RGB_to_Hex(self, rgb):
		color = '#'
		for i in rgb:
			num = int(i)
			color += str(hex(num))[-2:].replace('x', '0').upper()
		return color

	def bivariate_cmap(self, lower_left=None, lower_right=None, upper_left=None, upper_right=None):
		# if upper_right is None:
		# 	upper_right = [87, 66, 73]
		# if upper_left is None:
		# 	upper_left = [100, 172, 190]
		# if lower_right is None:
		# 	lower_right = [200, 90, 90]
		# if lower_left is None:
		# 	lower_left = [232, 232, 232]
		c00 = np.array(lower_left) / 255
		c10 = np.array(lower_right) / 255
		c01 = np.array(upper_left) / 255
		c11 = np.array(upper_right) / 255

		# Create color mesh
		n = 256
		colors = np.zeros((n, n, 3))
		for i in range(n):
			for j in range(n):
				# Bi-linear interpolation
				x = i / (n - 1)
				y = j / (n - 1)
				colors[j, i] = (1 - x) * (1 - y) * c00 + x * (1 - y) * c10 + (1 - x) * y * c01 + x * y * c11  # 双线性插值
		# Create color map
		cmap = mcolors.ListedColormap(colors.reshape(n * n, 3))
		return cmap, colors

	def add_icon(self, ax, icon_path, position='lower left', zoom=0.2):
		arr_img = plt.imread(icon_path)
		imagebox = OffsetImage(arr_img, zoom=zoom, resample=True)
		ab = AnnotationBbox(imagebox, (0.05, 0.1) if position == 'lower left' else (0.95, 0.1),
		                    xycoords='axes fraction',
		                    frameon=False)
		ax.add_artist(ab)

	def calc(self, data):
		global sd
		n = len(data)
		mean = 0.0  # Mean
		mean2 = 0.0  # Mean of x^2
		mean3 = 0.0  # Mean of x^3
		for a in data:
			mean += a
			mean2 += a ** 2
			mean3 += a ** 3
		mean /= n
		mean2 /= n
		mean3 /= n
		sd = math.sqrt(abs(mean2 - mean ** 2))
		return [mean, sd, mean3]

	def calc_stat(self, data):
		global skewness, kurtosis
		[mean, sd, mean3] = self.calc(data)
		n = len(data)
		mean4 = 0.0  # Used to set the numerator of kurtosis
		for a in data:
			a -= mean
			mean4 += a ** 4
		mean4 /= n
		if sd != 0:
			skewness = (mean3 - 3 * mean * sd ** 2 - mean ** 3) / (sd ** 3)  # Skewness
			kurtosis = mean4 / (sd ** 4)  # Kurtosis
		else:
			skewness = 0
			kurtosis = 0
		return n, mean, sd, skewness, kurtosis

	def delete_abnormal_values(self, df, column_drop):

		"""
		Remove the abnormal values
		Fill null values with the mean values
		"""

		# Load the Head of the DataFrame. Create a New List to Update the Head
		column_head = list(df)
		list_head = [i for i in column_head if i not in column_drop]
		df_copy = df.copy()

		# Select the Heavy Metal Columns
		for num, h in enumerate(list_head):
			q1 = df_copy[h].quantile(0.25)
			q3 = df_copy[h].quantile(0.75)
			iqr = q3 - q1
			# Determine the Range to Distinguish the Abnormal Values
			num_lower_limit = q1 - 1.5 * iqr
			num_upper_limit = q3 + 1.5 * iqr
			# Calculate the Mean Value of the Column
			num_mean = np.mean(df_copy[h])
			# Replace the Abnormal Values with the Mean Value of Corresponding Columns
			df_copy.loc[(df_copy[h] > num_upper_limit) | (df_copy[h] < num_lower_limit), h] = num_mean

		return df_copy

	def canvas_setup(self, fig_width=12, fig_height=7):
		# Allocate subplots
		fig = plt.figure(figsize=(fig_width, fig_height),
		                 constrained_layout=True)
		gs = GridSpec(nrows=3,
		              ncols=3,
		              height_ratios=[1, 1, 1],
		              width_ratios=[2, 2, 0.05],
		              hspace=0.3,
		              wspace=0.1)
		self.ax_scatter = [fig.add_subplot(gs[r, 0]) for r in range(0, 3)]
		self.ax_map = [fig.add_subplot(gs[r, 1]) for r in range(0, 3)]
		self.ax_map_legend = [fig.add_subplot(gs[r, 2]) for r in range(0, 3)]

	def plot_chord(self):
		"""
		Part 1. Chord Plot
		"""
		df_chordplot_path = self.output_directory + "Chordplot.xlsx"
		df_taxonomy_path = self.output_directory + "Taxonomy.xlsx"
		df_continent_path = self.output_directory + "Continent.xlsx"
		shp_continent_path = self.working_directory+"Continent\\Continent\\continent.shp"
		f_chord = plt.figure(figsize=(10, 10),
		               constrained_layout=True)
		f_chord.patch.set_alpha(0.0)
		gs_temp = GridSpec(nrows=1,
		                   ncols=1,
		                   height_ratios=[1],
		                   width_ratios=[1],
		                   hspace=0.3,
		                   wspace=0.1)

		self.ax_chord = f_chord.add_subplot(gs_temp[0, 0])
		# dict_chord_color = {k: self.RGB_to_Hex(v) for k, v in dict_color.items()}
		chordplot.chord_plot(tree_path=[self.tree_path_all],
		                     df_data_path=self.df_data_path,
		                     df_chordplot_path = df_chordplot_path,
		                     df_taxonomy_path = df_taxonomy_path,
		                     df_continent_path = df_continent_path,
		                     output_dir = self.working_directory+"2.1_Method\\",
		                     cmap=self.linearcmap,
		                     font_size=6,
		                     r=radius,
		                     ax=self.ax_chord,
		                     show_name=False,
		                     linecolor=self.RGB_to_Hex(self.list_color[0]),
		                     taxonomy_width=1,
		                     width_fraction=self.width_fraction,
		                     list_color=[self.RGB_to_Hex(color) for color in list_color],
		                     continent_path=shp_continent_path)

		f_chord.savefig(self.working_directory+f"2.1_Method\\chordplot_{today.strftime("%Y%m%d")}.png",dpi=300)
		f_chord.savefig(self.working_directory + f"2.1_Method\\chordplot_{today.strftime("%Y%m%d")}.svg")
		plt.close()
		print("The chord plot has been saved.")

	def plot_world(self):
		"""
		Part 2. Plot World Maps
		Parameter:
			mode: Choose between 'country' and 'scatter'
		"""
		# Prepare shape file dataframe
		world = gpd.read_file(self.shp_path)  # Read shape file
		world.columns = ['NAME_CHN', 'Country', 'NR_C', 'NR_C_ID', 'SOC', 'geometry']  # Rename head
		world.drop(world[world["Country"] == "ANTARCTICA"].index, axis=0, inplace=True)  # Drop Antarctica
		# list_country = list(world["Country"])  # Gather country names from shape file

		# Prepare data for chord plot
		df_chord = pd.read_excel(self.df_data_path, header=[0])[
			["Latitude",
			 "Longitude",
			 "Concentration-Gosner_Stage",
			 "log_Concentration-Physical_feature-SVL(mm)",
			 "log_Concentration-Physical_feature-Body_mass(g)",
			 "log_Concentration-Physical_feature-CF(100*BM(g)/SVL(mm)³)"]]  # Read excel
		df_temp = df_chord[
			df_chord["Concentration-Gosner_Stage"] > 46].copy()  # Extract metamorphosed (Gosner Stage > 46) data
		df_temp.drop("Concentration-Gosner_Stage", axis=1, inplace=True)  # Drop GS column
		df_temp.columns = ["Lat", "Lon", "log_SVL", "log_BM", "log_CF"]  # Rename column head
		df_temp = self.delete_abnormal_values(df_temp, column_drop=["Lat", "Lon"])
		list_BS = ["log_SVL", "log_BM", "log_CF"]

		for n, BS in enumerate(list_BS):
			scatter_data = df_temp[["Lat", "Lon", BS]].copy()

			# Calculate Statistic for scatter size on the world map
			min_val = scatter_data[BS].min()
			max_val = scatter_data[BS].max()
			if BS == "log_SVL":
				scatter_data.loc[:,BS + "_scaled"] = (scatter_data[BS] - min_val) / (max_val - min_val) * 60
			elif BS == "log_BM":
				scatter_data.loc[:,BS + "_scaled"] = (scatter_data[BS] - min_val) / (max_val - min_val) * 90
			elif BS == "log_CF":
				scatter_data.loc[:,BS + "_scaled"] = (scatter_data[BS] - min_val) / (max_val - min_val) * 80

			# Prepare data
			scatter_gdf = gpd.GeoDataFrame(scatter_data,
			                               geometry=gpd.points_from_xy(scatter_data['Lon'],
			                                                           scatter_data['Lat']))
			# Plot world map
			world.plot(ax=self.ax_map[n],
			           color="#f5f5f5",
			           edgecolor='gray',
			           linewidth=0.1)  # World map
			scatter_map = self.ax_map[n].scatter(scatter_gdf.geometry.x,  # Longtitude
			                                     scatter_gdf.geometry.y,  # Latotide
												 c=scatter_gdf[BS],  # Color changed with values
												 cmap=self.linearcmap,  # Color mapping
												 s=scatter_gdf[BS + "_scaled"] * 0.8,  # Define scatter size
												 alpha=0.5 if BS!= "log_CF" else 0.7)  # Set alpha
			self.ax_map[n].axhline(0,
			                       linestyle="--",
			                       linewidth=0.5,
			                       color="gray")
			self.ax_map[n].axhline(-30,
			                       linestyle="--",
			                       linewidth=0.2,
			                       color="gray")
			self.ax_map[n].axhline(30,
			                       linestyle="--",
			                       linewidth=0.2,
			                       color="gray")
			self.ax_map[n].axis("off")

			cbar = plt.colorbar(scatter_map, cax=self.ax_map_legend[n])
			cbar.set_label(BS)

	def plot_world_itv(self, mode="scatter"):
		"""
		Part 2. Plot World Maps
		Parameter:
			mode: Choose between 'country' and 'scatter'
		"""
		# Prepare shape file dataframe
		world = gpd.read_file(self.shp_path)  # Read shape file
		world.columns = ['NAME_CHN', 'Country', 'NR_C', 'NR_C_ID', 'SOC', 'geometry']  # Rename head
		world.drop(world[world["Country"] == "ANTARCTICA"].index, axis=0, inplace=True)  # Drop Antarctica
		# list_country = list(world["Country"])  # Gather country names from shape file

		if mode == "country":
			# Prepare data for chord plot
			df_chord = pd.read_excel(self.df_data_path, header=[0])[
				["Country",
                 "Concentration-Gosner_Stage",
                 "log_Concentration-Physical_feature-SVL(mm)",
                 "log_Concentration-Physical_feature-Body_mass(g)",
                 "log_Concentration-Physical_feature-CF(100*BM(g)/SVL(mm)³)"]]  # Read excel
			df_temp = df_chord[df_chord["Concentration-Gosner_Stage"] > 46].copy().drop("Concentration-Gosner_Stage",
			                                                                     axis=1)  # Extract metamorphosed (Gosner Stage > 46) data

			df_temp.columns = ["Country", "log_SVL", "log_BM", "log_CF"]  # Rename column head
			df_melt = df_temp.melt(id_vars=['Country'],
			                       var_name='Feature',
			                       value_name='values')  # Transfer dataframe
			list_country_documented = list(set(df_temp["Country"]))  # Get documented country names

			# Calculate interspecies trait variance (ITV)
			df_mean = df_melt.pivot_table(index=["Country"],
			                              columns="Feature",
			                              values="values",
			                              aggfunc="mean")  # Calculate mean body size indicators on the country scale
			df_itv = pd.DataFrame(data=None,
			                      columns=["Country", "ITV_SVL", "ITV_BM", "ITV_CF"])  # Create an empty dataframe to contain ITVs
			for lc in list_country_documented:  # Iterate each country
				current_row = df_itv.shape[0]  # Current row we are at
				df_itv.loc[current_row, "Country"] = lc  # Fill the current country into the current row
				for lf in list_feature:  # Iterate each body size indicators
					list_data = list(df_melt[(df_melt["Country"] == lc) & (df_melt["Feature"] == lf)][
						                 "values"])  # Extract indicator data of certain countries
					n, mean, std, skew, kurt = self.calc_stat(list_data)  # Prepare basic parameters
					CV1 = std / mean
					CV4 = CV1 - CV1 ** 3 / n + CV1 / 4 / n + CV1 ** 2 * skew / 2 / n + CV1 * kurt / 8 / n
					df_itv.loc[current_row, "ITV_" + lf] = CV4  # Fill the ITVs into dataframs
			df_itv.set_index("Country", inplace=True)  # Set column of country as index

			# Connect dataframe df_itv to world
			for lf in list_feature:  # Iterate each feature
				for lc in list_country_documented:  # Iterate each country
					world.loc[world["Country"] == lc, lf] = df_mean.loc[
						lc, lf]  # Fill indicators in shape file dataframe
					world.loc[world["Country"] == lc, "ITV_" + lf] = df_itv.loc[
						lc, "ITV_" + lf]  # Fill ITVs in shape file dataframe
					world.loc[:,"norm_" + lf] = (world[lf] - world[lf].min()) / (
							world[lf].max() - world[lf].min())  # Normalize indicators
					world.loc[:,"norm_ITV_" + lf] = (world["ITV_" + lf] - world["ITV_" + lf].min()) / (
							world["ITV_" + lf].max() - world["ITV_" + lf].min())  # Normalize ITVs
					world.loc[:,'color_' + lf] = world["norm_" + lf] + world["norm_ITV_" + lf]  # Used to get color

			# Preparation
			# list_text = ["SVL (Snout-vent Length, mm)", "BM (Body Weight, g)",
			#              "CF (Condition Factor = BM / SVL³)"]  # Titles for each feature
			list_color = [col for col in list(world) if "color" in col]  # Get color values
			# list_ITV = [col for col in list(world) if "ITV" in col]  # Get ITVs
			dict_note = {"ITV SVL": "ITV$_S$$_V$$_L$",
			             "ITV BM": "ITV$_B$$_M$",
			             "ITV CF": "ITV$_C$$_F$",
			             "Mean SVL": "Mean$_S$$_V$$_L$ ",
			             "Mean BM": "Mean$_B$$_M$ ",
			             "Mean CF": "Mean$_C$$_F$ "}

			# Plot world maps
			for n, lc in enumerate(list_color):  # Iterate colors
				# Plot world map
				world.plot(ax=self.ax_map[n],
				           column=lc,
				           cmap=self.custom_cmap,
				           edgecolor='black',
				           linewidth=0.1,
				           missing_kwds={"color": "#f5f5f5", "edgecolor": "gray", "hatch": ""})
				# self.ax_map[n].text(0, -75, list_text[n], va='top', ha='center')  # Add titles
				self.ax_map[n].axis("off")  # Turn off axes
				self.ax_map[n].set_ylim(-70, 90)

				N = 256  # Define fineness of legends
				Z = np.zeros((N, N, 3))  # Create a zero-like array
				for i in range(N):  # Iterate each row
					for j in range(N):  # Iterate each column
						x = i / (N - 1)  # Calculate percentage in the current row
						y = j / (N - 1)  # Calculate percentage in the current column
						Z[j, i] = ((1 - x) * (1 - y) * np.array(self.list_color[int(len(self.list_color)/2)]) /
						           N + x * (1 - y) * np.array(self.list_color[int(len(self.list_color)/4)]) /
						           N + (1 - x) * y * np.array(self.list_color[int(len(self.list_color)/2)]) /
						           N + x * y * np.array(upper_right) / N)  # Bi-linear interpolation

				self.ax_map_legend[n].imshow(Z,
				                             cmap=self.custom_cmap,
				                             origin='lower',
				                             extent=[0.5, 1, 0.25, 0.75])  # Plot legendS

				self.ax_map_legend[n].text(x = 0.75,
				                           y = 0.1,
				                           s = dict_note[lc.replace("color_", "ITV ")],
				                           ha="center",
				                           va="top",
				                           fontsize=10,
				                           rotation=90)  # Set label on x axis
				self.ax_map_legend[n].text(x = 0.5,
				                           y = 0.5,
				                           s = dict_note[lc.replace("color_", "Mean ")],
				                           fontsize=10,
				                           ha="right",
				                           va="center",
				                           rotation=0)  # Set label on y axis
				# # Add ticks of legends
				self.ax_map_legend[n].set_xticks([0, 1], [], fontsize=10)  # Set legend ticks on x-axis
				self.ax_map_legend[n].set_yticks([0, 1], [], fontsize=10)  # Set legend ticks on y-axis
				self.ax_map_legend[n].text(x = 0.5,
				                           y = 0.22,
				                           s = f'{world[lc.replace("color_", "ITV_")].min():.2f}',
				                           ha="left",
				                           va="top",
				                           fontsize=10,
				                           rotation=90)
				self.ax_map_legend[n].text(x = 1,
				                           y = 0.22,
				                           s = f'{world[lc.replace("color_", "ITV_")].max():.2f}',
				                           ha="center",
				                           va="top",
				                           fontsize=10,
				                           rotation=90)
				self.ax_map_legend[n].text(x = 0.5,
				                           y = 0.25,
				                           s = f'{world[lc.replace("color_", "")].min():.2f}',
				                           ha="right",
				                           va="bottom",
				                           fontsize=10)
				self.ax_map_legend[n].text(x = 0.5,
				                           y = 0.75,
				                           s = f'{world[lc.replace("color_", "")].max():.2f}',
				                           ha="right",
				                           va="top",
				                           fontsize=10)
				self.ax_map_legend[n].axis("off")  # Turn off the axes

		elif mode == "scatter":
			# Prepare data for chord plot
			df_chord = pd.read_excel(self.df_data_path, header=[0])[
				["Latitude",
	            "Longitude",
	            "Concentration-Gosner_Stage",
	            "log_Concentration-Physical_feature-SVL(mm)",
	            "log_Concentration-Physical_feature-Body_mass(g)",
	            "log_Concentration-Physical_feature-CF(100*BM(g)/SVL(mm)³)"]]  # Read excel
			df_temp = df_chord[
				df_chord["Concentration-Gosner_Stage"] > 46].copy()  # Extract metamorphosed (Gosner Stage > 46) data
			df_temp.drop("Concentration-Gosner_Stage", axis=1, inplace=True)  # Drop GS column
			df_temp.columns = ["Lat", "Lon", "lg $\\it{SVL}$", "lg $\\it{BM}$", "lg $\\it{CF}$"]  # Rename column head
			df_temp = self.delete_abnormal_values(df_temp, column_drop=["Lat", "Lon"])
			list_BS = ["lg $\\it{SVL}$", "lg $\\it{BM}$", "lg $\\it{CF}$"]

			for n, BS in enumerate(list_BS):
				scatter_data = df_temp[["Lat", "Lon", BS]].copy()

				# Calculate Statistic for scatter size on the world map
				min_val = scatter_data[BS].min()
				max_val = scatter_data[BS].max()
				if BS == "lg $\\it{SVL}$":
					scatter_data.loc[:,BS + "_scaled"] = (scatter_data[BS] - min_val) / (max_val - min_val) * 60
				elif BS == "lg $\\it{BM}$":
					scatter_data.loc[:,BS + "_scaled"] = (scatter_data[BS] - min_val) / (max_val - min_val) * 90
				elif BS == "lg $\\it{CF}$":
					scatter_data.loc[:,BS + "_scaled"] = (scatter_data[BS] - min_val) / (max_val - min_val) * 80

				# Prepare data
				scatter_gdf = gpd.GeoDataFrame(scatter_data,
				                               geometry=gpd.points_from_xy(scatter_data['Lon'],
				                                                           scatter_data['Lat']))

				# Plot world map
				world.plot(ax=self.ax_map[n],
				           color="#f5f5f5",
				           edgecolor='gray',
				           linewidth=0.1)  # World map
				scatter_map = self.ax_map[n].scatter(scatter_gdf.geometry.x,  # Longtitude
													 scatter_gdf.geometry.y,  # Latotide
													 c=scatter_gdf[BS],  # Color changed with values
													 cmap=self.linearcmap,  # Color mapping
													 s=scatter_gdf[BS + "_scaled"],  # Define scatter size
													 alpha=0.5 if BS!= "lg $\\it{CF}$" else 0.7,  # Set alpha
				                                     )
				self.ax_map[n].axhline(0,
				                       linestyle="--",
				                       linewidth=0.5,
				                       color="gray")
				self.ax_map[n].axhline(-30,
				                       linestyle="--",
				                       linewidth=0.2,
				                       color="gray")
				self.ax_map[n].axhline(30,
				                       linestyle="--",
				                       linewidth=0.2,
				                       color="gray")
				self.ax_map[n].axis("off")

				cbar = plt.colorbar(scatter_map, cax=self.ax_map_legend[n])
				cbar.set_label(BS)

	def plot_scatter(self):
		y_lower = -70
		y_upper = 90
		dict_illustration = {"lg $\\it{SVL}$": ["Shorter", "Longer"],
		                     "lg $\\it{BM}$": ["Lighter", "Fatter"],
		                     "lg $\\it{CF}$": ["Puffier", "Denser"]}
		list_color = ['color $\\it{SVL}$', 'color $\\it{BM}$', 'color $\\it{CF}$']  # Get color values
		y_range = [y_lower, y_upper]
		degree = 2
		df_scatter = pd.read_excel(self.df_data_path, header=[0])[
			["Latitude",
			 "Concentration-Gosner_Stage",
			 "log_Concentration-Physical_feature-SVL(mm)",
			 "log_Concentration-Physical_feature-Body_mass(g)",
			 "log_Concentration-Physical_feature-CF(100*BM(g)/SVL(mm)³)",
			 "Concentration-Physical_feature-ITV_SVL",
			 "Concentration-Physical_feature-ITV_BM",
			 "Concentration-Physical_feature-ITV_CF"]]  # Read excel
		df_scatter = (df_scatter[df_scatter["Concentration-Gosner_Stage"] > 46].
		              drop("Concentration-Gosner_Stage", axis=1))  # Extract metamorphosed (Gosner Stage > 46) data
		list_head = ["Lat",
		             "lg $\\it{SVL}$",
		             "lg $\\it{BM}$",
		             "lg $\\it{CF}$",
		             "$\\it{ITV}$$_{SVL}$",
		             "$\\it{ITV}$$_{BM}$",
		             "$\\it{ITV}$$_{CF}$"]
		dict_pair = {"lg $\\it{SVL}$": "$\\it{ITV}$$_{SVL}$",
		             "lg $\\it{BM}$": "$\\it{ITV}$$_{BM}$",
		             "lg $\\it{CF}$": "$\\it{ITV}$$_{CF}$"}
		df_scatter.columns = list_head
		list_alpha = [0.02, 0.02, 0.02]
		list_lat = df_scatter["Lat"]
		list_BS = ["lg $\\it{SVL}$", "lg $\\it{BM}$", "lg $\\it{CF}$"]
		fit_mode = "quadratic"
		for n, (BS, lc) in enumerate(zip(list_BS, list_color)):
			if fit_mode == "quadratic":

				list_bs = df_scatter[BS]

				coefficients = np.polyfit(list_lat, list_bs, degree)
				poly = np.poly1d(coefficients)

				# Predict y
				XX = np.linspace(min(list_lat), max(list_lat), 100)
				y_pred = poly(XX)

				# Calculate standard error
				num = len(list_lat)
				list_hat = poly(list_lat)
				residuals = list_bs - list_hat
				sigma = np.sqrt(np.sum(residuals ** 2) / np.abs(n - degree - 1))

				# Calculate Metrix
				X_design = np.column_stack([list_lat ** i for i in range(degree + 1)])
				X_design_XX = np.column_stack([XX ** i for i in range(degree + 1)])

				# Covariance Metrix
				cov_matrix = sigma ** 2 * np.linalg.inv(np.dot(X_design.T, X_design))

				# Standard Error
				y_pred_se = np.sqrt(np.sum(np.dot(X_design_XX, cov_matrix) * X_design_XX, axis=1))

				# Calculate confidence interval
				alpha = 0.05  # 95% Confidence Interval
				t_value = t.ppf(1 - alpha / 2, num - degree - 1)
				ci = t_value * y_pred_se

				# R^2
				y_mean = np.mean(list_bs)
				ss_res = np.sum(residuals ** 2)
				ss_tot = np.sum((list_bs - y_mean) ** 2)
				r_squared = 1 - (ss_res / ss_tot)

				# Symmetric axis
				a = coefficients[0]
				b = coefficients[1]
				symmetric_axis = -b / 2 / a

				para_text = ""
				list_para = [chr(i) for i in range(97, 123)]

				for i, coeff in enumerate(coefficients):
					num = str(f"{coeff:.2e}").split("e")
					new_num = num[0] if coeff < 0 else " " + num[0]
					num[1] = int(num[1])
					if num[1] != 0:
						new_num += "×10"
						for numero in str(num[1]):
							new_num += f"$^{numero}$"
					para_text += f"{list_para[i]} = {new_num}\n"

				# Plot
				self.ax_scatter[n].set_ylim(y_range[0], y_range[1])
				# self.ax_scatter[n].set_xlim(-4, 3)
				self.ax_scatter[n].axhline(0,
				                           linestyle="--",
				                           linewidth=0.5,
				                           color="gray")
				self.ax_scatter[n].axhline(-30,
				                           linestyle="--",
				                           linewidth=0.2,
				                           color="gray")
				self.ax_scatter[n].axhline(30,
				                           linestyle="--",
				                           linewidth=0.2,
				                           color="gray")
				self.ax_scatter[n].axhline(symmetric_axis,
				                           linestyle="--",
				                           linewidth=0.5,
				                           color="red")
				self.ax_scatter[n].fill_betweenx(XX,
				                                 y_pred - ci,
				                                 y_pred + ci,
				                                 color='gray',
				                                 alpha=0.2,
				                                 linewidth=0)
				self.ax_scatter[n].scatter(list_bs,
				                           list_lat,
				                           label=BS,
				                           alpha=list_alpha[n],
				                           c=df_scatter[dict_pair[BS]],
				                           cmap=self.linearcmap)
				self.ax_scatter[n].plot(y_pred,
				                        XX,
				                        label='Data',
				                        color=self.RGB_to_Hex(self.list_color[0]),
				                        linewidth=1)
				self.ax_scatter[n].text(np.mean(self.ax_scatter[n].get_xlim()),
				                        -90,
				                        BS,
				                        fontsize=10,
				                        va='top',
				                        ha="center")
				self.ax_scatter[n].text(self.ax_scatter[n].get_xlim()[0],
				                        70,
				                        f"{para_text}R² = {r_squared:.4f}",
				                        fontsize=10,
				                        va='top',
				                        ha="left")
				self.ax_scatter[n].text(self.ax_scatter[n].get_xlim()[1],
				                        symmetric_axis - 3,
				                        f"{symmetric_axis:.2f}°N",
				                        fontsize=10,
				                        color="red",
				                        va='top',
				                        ha="right")

			elif fit_mode == "segment":
				k = 1
				kmeans = KMeans(n_clusters=k)
				kmeans.fit(df_scatter[["Lat"]])
				breakpoints = kmeans.cluster_centers_.ravel()
				list_breakpoints = [min(df_scatter["Lat"])]
				for i in breakpoints:
					list_breakpoints.append(i)
				list_breakpoints.append(max(df_scatter["Lat"]))

				models = []
				annotation = ''
				x_values = np.linspace(min(df_scatter["Lat"]), max(df_scatter["Lat"]), 100)
				for i in range(k + 1):
					df_segment = df_scatter[
						(df_scatter["Lat"] > list_breakpoints[i]) & (df_scatter["Lat"] <= list_breakpoints[i + 1])].copy()

					model = LinearRegression()
					model.fit(df_segment[["Lat"]], df_segment[BS])
					models.append(model)
					lower_bound = list_breakpoints[i]
					upper_bound = list_breakpoints[i + 1]

					mask = (x_values >= lower_bound) & (x_values <= upper_bound)
					segment_x = x_values[mask]

					if len(segment_x) > 0:
						slope, intercept, r_value, p_value, std_err = linregress(df_segment["Lat"], df_segment[BS])
						equation = f'lg {BS} = {slope:.2f}x + {intercept:.2f}'
						r_squared = f'R² = {r_value ** 2:.2f}'
						p_value_text = f'p = {p_value:.2e}' if p_value > 0.01 else "p < 0.01"
						annotation += f"{equation}  {r_squared}  {p_value_text}\n"
						segment_x_reshaped = segment_x.reshape(-1, 1)
						y_pred = model.predict(segment_x_reshaped)
						self.ax_scatter[n].plot(y_pred,
						                        segment_x,
						                        label=f'Segment {i + 1} Fit',
						                        linewidth=1,
						                        color="red")
				self.ax_scatter[n].set_ylim(y_range[0], y_range[1])
				self.ax_scatter[n].scatter(df_scatter[BS],
				                           df_scatter["Lat"],
				                           label=BS,
				                           alpha=list_alpha[n],
				                           color=color)
				self.ax_scatter[n].axhline(0,
				                           linestyle="--",
				                           linewidth=0.5,
				                           color="gray")
				self.ax_scatter[n].axhline(-30,
				                           linestyle="--",
				                           linewidth=0.2,
				                           color="gray")
				self.ax_scatter[n].axhline(30,
				                           linestyle="--",
				                           linewidth=0.2,
				                           color="gray")
				self.ax_scatter[n].text(self.ax_scatter[n].get_xlim()[0],
				                        y_upper,
				                        f"{annotation}",
				                        va="top",
				                        ha="left")
				for symmetric_axis in breakpoints:
					annot = f"{symmetric_axis:.3f}° N" if symmetric_axis > 0 else f"{symmetric_axis:.3f}° S"
					self.ax_scatter[n].axhline(symmetric_axis,
					                           linestyle="--",
					                           linewidth=0.5,
					                           color="red")
					self.ax_scatter[n].text(self.ax_scatter[n].get_xlim()[1],
					                        symmetric_axis,
					                        annot,
					                        va="top",
					                        ha='right',
					                        color="red")

			self.ax_scatter[n].spines["right"].set_visible(False)
			self.ax_scatter[n].spines["left"].set_visible(False)
			self.ax_scatter[n].spines["top"].set_visible(False)
			self.ax_scatter[n].get_yaxis().set_visible(False)
			self.add_icon(self.ax_scatter[n],
			              self.illustration_path + dict_illustration[lc.replace("color", "lg")][0] + ".png",
			              position='lower left',
			              zoom=0.1)
			self.add_icon(self.ax_scatter[n],
			              self.illustration_path + dict_illustration[lc.replace("color", "lg")][1] + ".png",
			              position='lower right',
			              zoom=0.1)

	def plot_scatter_itv(self):
		y_lower = -70
		y_upper = 90
		list_xlim = [[-1, 2], [-0.87, 2], [0, 0.5]] # Axes X boundary for each subplot
		y_range = [y_lower, y_upper]
		self.list_colors = [self.linearcmap(i) for i in np.linspace(0, 1, 7)]
		index_order = [chr(i) for i in range(97, 123)]  # Serial strings
		degree = 2
		df_scatter = pd.read_excel(self.df_data_path, header=[0])[
			["Latitude",
			 "Concentration-Gosner_Stage",
			 "Concentration-Physical_feature-ITV_SVL",
			 "Concentration-Physical_feature-ITV_BM",
			 "Concentration-Physical_feature-ITV_CF",
			 "log_Concentration-Physical_feature-SVL(mm)",
			 "log_Concentration-Physical_feature-Body_mass(g)",
			 "log_Concentration-Physical_feature-CF(100*BM(g)/SVL(mm)³)"]]  # Read excel
		df_scatter = df_scatter[df_scatter["Concentration-Gosner_Stage"] > 46].\
					 drop("Concentration-Gosner_Stage", axis=1)  # Extract metamorphosed (Gosner Stage > 46) data
		list_head = ["Lat",
		             "$\\it{ITV}$$_{SVL}$",
		             "$\\it{ITV}$$_{BM}$",
		             "$\\it{ITV}$$_{CF}$",
		             "lg $\\it{SVL}$",
		             "lg $\\it{BM}$",
		             "lg $\\it{CF}$"]
		dict_pair = {"$\\it{ITV}$$_{SVL}$": "lg $\\it{SVL}$",
		             "$\\it{ITV}$$_{BM}$": "lg $\\it{BM}$",
		             "$\\it{ITV}$$_{CF}$": "lg $\\it{CF}$"}
		df_scatter.columns = list_head

		list_lat = df_scatter["Lat"]
		list_BS = ["$\\it{ITV}$$_{SVL}$", "$\\it{ITV}$$_{BM}$", "$\\it{ITV}$$_{CF}$",]
		for n, (BS, color) in enumerate(zip(list_BS, self.list_colors)):
			list_bs = df_scatter[BS]

			coefficients = np.polyfit(list_lat, list_bs, degree)
			poly = np.poly1d(coefficients)

			# Predict y
			XX = np.linspace(min(list_lat), max(list_lat), 100)
			y_pred = poly(XX)

			# Calculate standard error
			num = len(list_lat)
			list_hat = poly(list_lat)
			residuals = list_bs - list_hat
			sigma = np.sqrt(np.sum(residuals ** 2) / (num - degree - 1))

			# Calculate Metrix
			X_design = np.column_stack([list_lat ** i for i in range(degree + 1)])
			X_design_XX = np.column_stack([XX ** i for i in range(degree + 1)])

			cov_matrix = sigma ** 2 * np.linalg.inv(np.dot(X_design.T, X_design))   # Covariance Metrix

			y_pred_se = np.sqrt(np.sum(np.dot(X_design_XX, cov_matrix) * X_design_XX, axis=1))  # Standard Error

			alpha = 0.05  # 95% Confidence Interval
			t_value = t.ppf(1 - alpha / 2, num - degree - 1)
			ci = t_value * y_pred_se    # Calculate confidence interval

			# R^2
			y_mean = np.mean(list_bs)
			ss_res = np.sum(residuals ** 2)
			ss_tot = np.sum((list_bs - y_mean) ** 2)
			r_squared = 1 - (ss_res / ss_tot)

			# Symmetric axis
			a = coefficients[0]
			b = coefficients[1]
			symmetric_axis = -b / 2 / a
			axis_annot = f"{symmetric_axis:.2f}°N" if symmetric_axis > 0 else f"-{symmetric_axis:.2f}°S"

			para_text = ""
			list_para = [chr(i) for i in range(97, 123)]

			for i, coeff in enumerate(coefficients):
				num = str(f"{coeff:.2e}").split("e")
				new_num = num[0] if coeff < 0 else " " + num[0]
				num[1] = int(num[1])
				if num[1] != 0:
					new_num += "×10"
					for numero in str(num[1]):
						new_num += f"$^{numero}$"
				para_text += f"{list_para[i]} = {new_num}\n"

			# Plot
			self.ax_scatter[n].set_ylim(y_range[0],
		                                y_range[1])
			self.ax_scatter[n].set_xlim(list_xlim[n][0],
		                                list_xlim[n][1])
			self.ax_scatter[n].axhline(0,
		                               linestyle="--",
		                               linewidth=0.5,
		                               color="gray")
			self.ax_scatter[n].axhline(-30,
		                               linestyle="--",
		                               linewidth=0.2,
		                               color="gray")
			self.ax_scatter[n].axhline(30,
		                               linestyle="--",
		                               linewidth=0.2,
		                               color="gray")
			self.ax_scatter[n].axhline(symmetric_axis,
		                               linestyle="--",
		                               linewidth=0.3,
		                               color="red")
			self.ax_scatter[n].vlines(x=6.45,
		                              ymin=-50,
		                              ymax=50,
		                              linestyle="-",
		                              linewidth=0.5,
		                              color="black")
			self.ax_scatter[n].hlines(-50,
			                           self.ax_scatter[n].get_xlim()[1] * 0.97 +self.ax_scatter[n].get_xlim()[0] * 0.03,
			                           self.ax_scatter[n].get_xlim()[1],
			                          color="black",
			                          linewidth=0.5)
			self.ax_scatter[n].hlines(0,
			                          self.ax_scatter[n].get_xlim()[1] * 0.97 +self.ax_scatter[n].get_xlim()[0] * 0.03,
			                          self.ax_scatter[n].get_xlim()[1],
			                          color="black",
			                          linewidth=0.5)
			self.ax_scatter[n].hlines(50,
			                          self.ax_scatter[n].get_xlim()[1] * 0.97 +self.ax_scatter[n].get_xlim()[0] * 0.03,
			                          self.ax_scatter[n].get_xlim()[1],
			                          color="black",
			                          linewidth=0.5)
			self.ax_scatter[n].fill_betweenx(XX,            # X range
			                                 y_pred - ci,   # Lower Y
			                                 y_pred + ci,   # Upper Y
			                                 color='gray',  # Color of the confidence interval
			                                 alpha=0.05,
			                                 linewidth=0)
			self.ax_scatter[n].scatter(list_bs,                     # X coordinates
		                               list_lat,                    # Y coordinates
		                               label=BS,                    # Label
		                               alpha=0.02,         # Alpha of the scatters
		                               c=df_scatter[dict_pair[BS]], # Color
		                               cmap=self.linearcmap)        # Color map
			self.ax_scatter[n].plot(y_pred,
		                            XX,
		                            label='Data',
		                            color=self.RGB_to_Hex(self.list_color[0]),
		                            linewidth=1)
			self.ax_scatter[n].text(x=list_xlim[n][0],
		                            y=70,
		                            s=f"{para_text}R² = {r_squared:.4f}",
		                            fontsize=10,
		                            va='top',
		                            ha="left")
			self.ax_scatter[n].text(x=list_xlim[n][1],
		                            y=symmetric_axis - 3,
		                            s=axis_annot,
		                            fontsize=10,
		                            color="red",
		                            va='top',
		                            ha="right")
			self.ax_scatter[n].text(x = self.ax_scatter[n].get_xlim()[1] * 0.97 +self.ax_scatter[n].get_xlim()[0] * 0.03,
		                            y = -50,
		                            s = "-50",
		                            ha="right",
		                            va="center")
			self.ax_scatter[n].text(x = self.ax_scatter[n].get_xlim()[1] * 0.97 +self.ax_scatter[n].get_xlim()[0] * 0.03,
		                            y = 50,
		                            s = "50",
		                            ha="right",
		                            va="center")
			self.ax_scatter[n].spines["right"].set_visible(False)
			self.ax_scatter[n].spines["left"].set_visible(False)
			self.ax_scatter[n].spines["top"].set_visible(False)
			self.ax_scatter[n].set_xlabel(f"{BS}", fontsize=10)
			self.ax_scatter[n].yaxis.set_visible(False)
			x = list_xlim[n][0] - abs(list_xlim[n][1] - list_xlim[n][0]) * 0.01
			y = y_upper
			self.ax_scatter[n].text(x = x,
		                            y = y,
		                            s = f'({index_order[n]})',
		                            fontsize = 10,
		                            va="top",
		                            ha="right",
		                            weight="bold")  # Add serial strings
			self.ax_scatter[n].text(x = 6.5,
		                            y = 0,
		                            s = "Latitude",
		                            ha="left",
		                            va="center",
		                            rotation=90)

	def plot_tree(self):
		def trim_tips(tree_2trim, node_2trim):
			if node_2trim.is_terminal():
				parent = tree_2trim.prune(node_2trim)
				trim_tips(tree_2trim,parent)

		"""
		Part 3. Plot trees and H-test results
		"""

		f = plt.figure(figsize=(12, 10), constrained_layout=True)
		gs_temp = GridSpec(1, 3, height_ratios=[1], width_ratios=[1, 1, 0.04], hspace=0.1, wspace=0.1)

		self.ax_tree = f.add_subplot(gs_temp[0, 0])
		self.ax_heatmap = f.add_subplot(gs_temp[0, 1])
		self.ax_legend = f.add_subplot(gs_temp[0,2])

		# Prepare indicator data
		df = pd.read_excel(self.df_data_path, header=[0])  # Read whole file
		df = df[df["Concentration-Gosner_Stage"] > 46]  # Screen and get metamorphosed anurans data
		df = df[["Country",
		         "Latitude",
		         "Longitude",
		         'Species',
		         "log_Concentration-Physical_feature-SVL(mm)",
		         "log_Concentration-Physical_feature-Body_mass(g)",
		         "log_Concentration-Physical_feature-CF(100*BM(g)/SVL(mm)³)",
		         "Concentration-Physical_feature-ITV_SVL",
		         "Concentration-Physical_feature-ITV_BM",
		         "Concentration-Physical_feature-ITV_CF"
		         ]]
		df["Site"] = df[["Country", "Latitude","Longitude"]].astype(str).agg('-'.join, axis=1)
		df.drop(["Country", "Latitude", "Longitude"],
		        axis=1,
		        inplace=True)
		df.columns = ["S",
		              "lg $\\it{SVL}$",
		              "lg $\\it{BM}$",
		              "lg $\\it{CF}$",
		              "ITV$_{SVL}$",
		              "ITV$_{BM}$",
		              "ITV$_{CF}$",
		              "Site"]  # Rename column head

		list_species = list(set(df["S"]))  # Get a list of unique species
		list_indicator = ["lg $\\it{SVL}$",
		                  "lg $\\it{BM}$",
		                  "lg $\\it{CF}$",
		                  "ITV$_{SVL}$",
		                  "ITV$_{BM}$",
		                  "ITV$_{CF}$"]  # Features of body size

		# Prepare the tree file (Continental)
		tree = Phylo.read(self.tree_path_all, "newick")  # Read tree file
		tree.rooted = True
		list_s_tree = [str(tip.name).replace("_", " ") for tip in
		               tree.get_terminals()]  # Get terminals that are available for continental scale H-test
		# H-test to get species body size indicators differences among continents
		# If p_value < 0.05,refuse H0. We regard there is significant difference among groups.
		# Create an empty dataframe
		df_f_site = pd.DataFrame(data=None, # Empty data
		                         columns=["S"]+list_indicator)  # Column names
		for num, sp in enumerate(list_s_tree):  # Iterate each species
			if sp in list_species:  # If "sp" is the expect value
				df_species = df[df["S"] == sp]  # Get rows of certain species
				df_f_site.loc[len(df_f_site.index), "S"] = sp  # # Fill with current species
				for n, li in enumerate(list_indicator):  # Iterate each feature
					df_feature = df_species[["Site", li]]  # Get slices of dataframe
					list_site = list(set(df_feature["Site"]))  # Get a list of unique continents
					list_feature_values = list(set(df_feature[li]))    # Get a list of unique feature values to make sure H-test feasible
					if len(list_site) > 1 and len(list_feature_values)>1:  # Regard indicators in a same site as a group. If we have over two groups, then execute H-test.
						groups = [group[li].values for name, group in df_feature.groupby("Site")]
						statistic, p_value = kruskal(*groups)  # H-test
						df_f_site.loc[num, li] = np.float32(p_value)  # Fill with H-test results
						if n == 0: # When n equals to 0, "sp" has been trimmed, meaning no need to trim in the subsequent iterations (n>0).
							node = tree.find_any(name=sp.replace(" ","_"))  # Find the clade with the "sp"
							trim_tips(tree,node)    # Trim clade

			else:   # If "sp" not in studied df
				node = tree.find_any(name=sp.replace(" ", "_"))  # Find the clade with the "sp"
				trim_tips(tree, node)    # Trim clade

		df_f_site.dropna(axis=0, how="any", inplace=True)  # Drop rows with any empty cells
		list_s = list(df_f_site["S"])  # Get a list of remained species
		df_f_site.drop("S", axis=1, inplace=True)  # Drop the species column
		df_f_site = pd.DataFrame(np.float32(df_f_site)) # Convert precsion
		df_f_site.index = list_s  # Set index as list_s
		df_f_site = df_f_site.reindex(list_s_tree)  # Sort index according to tree terminals (continental)
		df_f_site.columns = list_indicator  # Rename column head
		df_f_site.dropna(how="any", axis=0, inplace=True)   # Delete rows with at least one empty cell
		df_f_site.to_excel(self.output_directory + "H-test_Site.xlsx")  # Save the dataframe
		print(f"{df_f_site.shape[0]} species were selected from {len(list_s_tree)} species.")

		# Draw phylogenetic tree
		Phylo.draw(tree,                                # Phylogenetic tree
		           axes=self.ax_tree,                   # Axes
		           do_show=False,                       # Show annotations or not
		           show_spine_bottom=False,             # Show bottom axis or not
		           show_taxon=True,                     # Show taxonomy or not
		           show_confidence=False,               # Show confidence or not
		           scale_factor=self.tree_scale_factor)
		# Draw heatmap
		sns.heatmap(df_f_site,                  # H-test results
		            annot=False,                # Show values in the cells or not
		            cmap=self.linearcmap,       # Color map
		            ax=self.ax_heatmap,         # Axes
		            cbar_ax = self.ax_legend,   # Axes of color bar legend
		            cbar=True,                  # Show color bar or not
		            fmt=".3f",                  # Decimal
		            cbar_kws={'format': '%.3f',
		                      'ticks': np.linspace(df_f_site.min().min(), df_f_site.max().max(), 3),
		                      "aspect": 6})
		for i in range(df_f_site.shape[1]):
			self.ax_heatmap.text(x = i+0.5,                     # X coordinates
			                     y = df_f_site.shape[0] + 0.5,  # Y coordinates
			                     s = list_indicator[i],         # Text
			                     ha = "center",                 # Horizontal alignment
			                     va = "top",                    # Vertical alignment
			                     color = "black")               # Text color
			for j in range(df_f_site.shape[0]):
				if df_f_site.iloc[j,i] <= 0.05:
					self.ax_heatmap.text(x= i+0.5,                          # X coordinates
										 y = j+0.5,                         # Y coordinates
					                     s = f"{df_f_site.iloc[j,i]:.3f}",  # Text
					                     ha='center',                       # Horizontal alignment
					                     va='center',                       # Vertical alignment
					                     color='white',                     # Text color
					                     fontsize = 8)                      # Font size
		self.ax_heatmap.yaxis.set_visible(False)
		self.ax_heatmap.xaxis.set_visible(False)
		plt.savefig(self.output_directory+f"heatmap_{today.strftime("%Y%m%d")}.png",dpi=300)
		plt.show()
		plt.close()

	def add_serial(self,list_axes):
		"""
		Part 4. Add serial strings
		"""
		index_order = [chr(i) for i in range(97, 123)]  # Serial strings
		# list_cor = [[-radius - 0.2, radius], [-210, 85], [-210, 85], [-210, 85], [-4.2, 85], [-4.2, 85], [-4.2, 85],
		#             [-0.8, 85], [-0.8, 85], [-0.8, 85], [0, 0]]  # Where to set serial strings
		for n, lo in enumerate(list_axes):  # Iterate each axis
			x_lim = lo.get_xlim()
			y_lim = lo.get_ylim()
			x = x_lim[0] - abs(x_lim[1] - x_lim[0]) * 0.01
			y = y_lim[1]
			lo.text(x, y, f'({index_order[n]})', fontsize=10, va="top", ha="right", weight="bold")  # Add serial strings


if __name__ == "__main__":
	radius = 1
	today = datetime.now()
	# formatted_date = today.strftime("%Y%m%d")
	formatted_date = 20250726
	plt.rcParams['legend.fontsize'] = 20
	plt.rcParams['font.sans-serif'] = "Arial"
	plt.rcParams["lines.linewidth"] = 0.5
	list_feature = ["log_SVL", "log_BM", "log_CF"]
	working_directory = "E:\\学习\\研二上\\重金属\\"
	output_directory = working_directory + "3.1_Body_Features\\"
	df_data_path = working_directory +f"2.2_Imputation\\Imputation{formatted_date}\\03Revised\\Preprocessed{formatted_date}.xlsx"
	tree_path_all = working_directory+"0.0_Fundamental_data\\Species_List\\Species_all20250512.nwk"

	width_fraction = [0.5, 0.3, 0.05, 0.15]
	tree_scale_factor = 0.7

	# list_color = [[39, 71, 83],      # Dark Green
	#               [41, 114, 112],  # Dark Turquoise
	#               [41, 157, 143],  # Turquoise
	#               [138, 176, 124],  # Green
	#               [231, 198, 107],  # Yellow
	#               [243, 163, 97],  # Orange
	# 			  [230, 109, 80],   # Red
	# 			  ]
	# list_color = [[0, 71, 105],
	# 			  [64, 128, 148],
	#               [119, 179, 178],
	#               [187, 214, 198],
	#               [255, 228, 168],
	#               [242, 217, 172],
	#               [222, 160, 145],
	#               [198, 91, 63],
	#               [105, 48, 14],
	#               ]
	# list_color = [[78,101,155],
	#               [138,140,191],
	#               [184,168,207],
	#               [231,188,198],
	#               [253,207,158],
	#               [239,164,132],
	#               [182,118,108]
	# ]
	# list_color = [[68,4,90],
	#               [65,62,133],
	#               [48,104,141],
	#               [31,146,139],
	#               [53,183,119],
	#               [145,213,66],
	#               [248,230,32]
	# ]
	list_color = [[14,91,118],
	              [26,134,163],
	              [70,172,202],
	              [155,207,232],
	              [205,205,164],
	              [255,202,95],
	              [254,168,9],
	              [253,152,2],
	              [251,132,2]
	]

	Body_size_indicators(list_color=list_color, working_directory = working_directory, str_df_path=df_data_path,
	                     tree_path_all=tree_path_all, tree_scale_factor=tree_scale_factor,width_fraction=width_fraction,
	                     indicator_type="BS")
	Body_size_indicators(list_color=list_color, working_directory = working_directory, str_df_path=df_data_path,
	                     tree_path_all=tree_path_all, tree_scale_factor=tree_scale_factor, width_fraction=width_fraction,
	                     indicator_type="ITV")

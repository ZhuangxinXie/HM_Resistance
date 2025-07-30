import os
import copy
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import geopandas as gpd
import statsmodels.formula.api as smf
import Bio.Phylo as Phylo
from tqdm import tqdm
import scipy.stats as stats
from scipy.linalg import cholesky
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, silhouette_samples
from matplotlib import transforms
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def build_phylo_vcv(tree):
	"""
	Build a phylogenetic covariance matrix from phylogenetic tree

	___________
	Parameters:
		tree: Bio.Phylo.BaseTree.Tree

	________
	Returns:
		numpy.ndarray: Covariance matrix
	"""

	list_nodes = tree.get_terminals()  # A list of tips
	n_sp = tree.count_terminals()  # Number of Species
	vcv = np.zeros((n_sp, n_sp))  # Create a zero-like matrix

	dict_node2root = {}  # A dictionary to store the distances from the nodes to the root

	def calc_dist2root(node, dict_current):
		dict_node2root[node] = dict_current
		for child in node.clades:
			calc_dist2root(child, dict_current + child.branch_length)

	calc_dist2root(tree.root, 0.0)

	for i, tip1 in enumerate(list_nodes):  # Iterate each row
		for j, tip2 in enumerate(list_nodes):  # Iterate each column
			if i == j:  # If cell [i,j] is on diagonal
				vcv[i, j] = dict_node2root[tip1]  # Get distances between tips and the tree root
			else:  # If cell [i,j] is NOT on diagonal
				mrca = tree.common_ancestor(tip1, tip2)  # Find common ancestor of "tip1" and "tip2"
				vcv[i, j] = dict_node2root[mrca]  # Get distances between mrca and the tree root

	dict_vcv = {}
	num_slice = 0
	for i in range(n_sp):
		sp1 = list_nodes[i].name
		for j in range(num_slice, n_sp):
			sp2 = list_nodes[j].name
			dict_vcv[frozenset({sp1, sp2})] = vcv[i, j]  # Make a set frozen so that it can be a hashable key
		num_slice += 1

	return dict_vcv, vcv


def compare(list_df_species, list_tree_species):
	"""
	Check if the species is compatible between the dataframe and the phylogenetic tree

	___________
	Parameters:
		list_df_species: A list of species in the dataframe
		list_tree_species: A list of species in the tree
	"""

	for sp in set(list_df_species):  # Iterate each species in the dataframe
		if sp not in list_tree_species:  # If this species is not in the tree
			raise ValueError(f"{sp} not in phylogenetic tree. It has been deleted in the dataframe.")


def build_full_phylo_cov(tree, list_df_species, list_tree_species, hm):
	"""
	Build a full phylogenetic covariance matrix according to the "list_df_species".

	___________
	Parameters:
		tree: A phylogenetic tree. All species in the tree can ce found in the "list_df_species"
		list_df_species: A list of species derived from a dataframe. The species in the list are allowed repeated.
		list_tree_species: A list of species derived from a phylogenetic tree
		hm: Heavy metals.

	_______
	Return:
		 array_cov: A n*n phylogenetic covariance matrix when the row number of "df" is n
	"""

	# Initialize
	compare(list_df_species,
	        list_tree_species)  # Compare species in both list to make sure every species (in "list_df_species") available from "list_tree_species"
	dict_vcv, vcv = build_phylo_vcv(tree)  # Get a dict that denotes the Patristic distance

	# Build a species matrix for dataframe as a covariance matrix
	num_rows = len(list_df_species)  # Get the number of rows
	df_cov = pd.DataFrame(np.zeros((num_rows, num_rows)))  # Build a zero-like dataframe
	array_df_species = np.array(
		list_df_species)  # Convert the list of species to an array, used to get col/row indexes
	for sp1 in tqdm(list_tree_species, desc=f"{hm}"):  # Iterate each column. Use "tqdm" to show the progress bar
		id_col = np.where(array_df_species == sp1)[0].tolist()  # Indexes represent "sp1"
		for sp2 in list_tree_species:  # Iterate each species to get the corresponding indexes.
			id_row = np.where(array_df_species == sp2)[0].tolist()  # Indexes represent "sp2"
			df_cov.iloc[id_row, id_col] = dict_vcv[
				frozenset({sp1, sp2})]  # Fill cells with values from "dict_vcv"

	array_cov = np.tril(np.asarray(df_cov))
	array_cov = array_cov + array_cov.T

	return array_cov


def fit_model(q, mod, x):
	res = mod.fit(q=q)  # statsmodels.formula.api
	return [q, res.params["Intercept"], res.params[x]] + res.conf_int().loc[x].tolist()


def phylo_quantreg(str_Y, str_X, tree, df, list_quantile, phylo_on=True):
	"""
	Quantile regression basing on a distance matrix revised with a phylogenetic covariance distance matrix

	___________
	Parameters:
		str_Y: Column name, used to designate the column "str_Y" as the dependent variable.
		str_X: Column name, used to designate the column "str_X" as the independent variable.
		tree: A phylogenetic tree
		df: A dataframe, including at least two columns of "str_Y" and "str_X"
		list_quantile: A list of quantile point
		phylo_on: Bool. Used to determine whether we should revise with a phylogenetic matrix

	_______
	Return:
		models: A dataframe storing regression parameters
	"""

	list_tree_species = [clade.name for clade in
	                     tree.get_terminals()]  # Get a list of unique species from phylogenetic tree
	phylo_cov = np.asarray(
		build_full_phylo_cov(tree, df["Species"], list_tree_species, str_Y))  # Get a full phylogenetic matrix
	# Make sure the matrix is positive definite matrix
	eigenvalues = np.linalg.eigvals(phylo_cov)
	if not np.allclose(phylo_cov, phylo_cov.T, 1e-10) or not np.all(eigenvalues > 0):
		phylo_cov = (phylo_cov + phylo_cov.T) / 2  # Compulsorily Convert
		print("The previous correlation matrix is not positive definite matrix. This matrix has been converted.")

	phylo_cov1 = cholesky(np.asarray(phylo_cov), lower=True)  # Convert into an array and Cholesky decompose
	phylo_cov_inv = np.linalg.inv(phylo_cov1)  # Transpose matrix

	# Make sure shape of phylogenetic covariance matches that of df
	n_samples = df.shape[0]  # Get data line number
	if phylo_cov1.shape != (n_samples, n_samples):  # If the shape is not (n_samples, n_samples)
		raise ValueError("Shape of 'phylo_cov' should match the sample size.")
	if phylo_on:
		df.loc[:, str_Y] = phylo_cov_inv @ df[str_Y]
		df.loc[:, str_X] = phylo_cov_inv @ df[str_X]
	mod = smf.quantreg(f"{str_Y} ~ {str_X}", df)
	models = [fit_model(X, mod, str_X) for X in list_quantile]
	models = pd.DataFrame(models, columns=["quantile", "intercept", "slope", "lb", "ub"])
	return models


def k_multi(tree, Y):
	"""
	Calculate K-value according to the following formulation (https://doi.org/10.1093/sysbio/syu030)

			  (Y-E(Y))^T·(Y-E(Y))
		    _______________________
			(Y-E(Y))^T·C^-1·(Y-E(Y))
	K =	_________________________________
			(tr(C)-N(1^T·C^-1·1)^-1)
			________________________
					  (N-1)

	where
		Y: A matrix with multiple features (index: species, column: features)
		E(Y): Expected value at the root of the phylogeny. It equals to (1^T·C^-1·1)^(-1)·(1^T·C^-1·1)
		tr(C): Trace of a covariance matrix
		N: Number of species
		1: A one-like matrix

	___________
	Parameters:
		tree: A phylogenetic tree
		Y: A matrix the values of the certain feature(s)

	_______
	Return:
		K: K-value
	"""

	list_taxa = [tip.name for tip in tree.get_terminals()]
	N = len(list_taxa)  # Number of species
	Y = np.array(Y.loc[list_taxa])

	# cov = build_covariance_matrix(tree)  # Phylogenetic matrix
	dict_vcv, vov = build_phylo_vcv(tree)
	cov_inv = np.linalg.inv(vov)  # Inverse matrix of phylogenetic matrix
	ones = np.ones((N, 1))  # Create a N * 1 matrix
	ones_T = ones.T  # Transpose matrix of one-like matrix
	dot_1 = np.linalg.inv(
		ones_T @ cov_inv @ ones)  # Inverse matrix of a dot product matrix("ones_T", "cov_inv" and "ones")
	dot_2 = ones_T @ cov_inv @ Y  # A dot product matrix ("ones_T", "cov_inv" and "Y")
	E = dot_1 @ dot_2  # Expected value at the root of the phylogeny
	# Numerator of K
	N_numerator_1 = Y.T - E.T  # (Y-E(Y))^T
	N_numerator_2 = N_numerator_1.T  # (Y-E(Y)
	N_Numerator = N_numerator_1 @ N_numerator_2  # (Y-E(Y))^T·(Y-E(Y)
	N_dominator_1 = N_numerator_1  # (Y-E(Y))^T
	N_dominator_2 = cov_inv  # C^-1
	N_dominator_3 = N_numerator_2  # (Y-E(Y))
	N_Dominator = N_dominator_1 @ N_dominator_2 @ N_dominator_3  # (Y-E(Y))^T·C^-1·(Y-E(Y))
	# Dominator of K
	D_numerator1 = np.trace(vov)
	D_numerator2 = N * dot_1
	D_Numerator = D_numerator1 - D_numerator2
	D_Dominator = N - 1

	K = np.divide(N_Numerator, N_Dominator) / np.divide(D_Numerator, D_Dominator)

	return list(K)


class Reg:
	def __init__(self,
	             str_working_directory,
	             str_tree_path,
	             list_heavy_metal,
	             list_quantile,
	             list_rgb_color,
	             n_clusters,
	             data_date,
	             feature):
		# Str
		self.str_feature = feature
		self.str_working_dir = str_working_directory
		self.str_output_directory = self.str_working_dir + "3.3_Correlation\\"
		self.str_dataframe_path = f"{self.str_working_dir}2.2_Imputation\\Imputation{data_date}\\03Revised\\Preprocessed{data_date}.xlsx"
		self.str_continent_path = self.str_working_dir + "\\Continent\\Continent\\continent.shp"
		self.str_world_path = self.str_working_dir + "Worldmap\\Worldmap.shp"
		self.str_tree_path = str_tree_path
		# Tree
		self.tree = Phylo.read(self.str_tree_path, "newick")
		# Dataframe
		self.df_quantile = None
		self.df = None
		self.df_pca = None
		self.df_loc = None
		self.df_radar = None
		# Int
		self.n_clusters = n_clusters
		self.data_date = data_date
		self.linear_cmap = None
		self.palette = None
		# List
		self.list_colors = None
		self.list_colors_cluster = None
		self.list_head = []
		self.list_quantile = list_quantile
		self.list_hm = list_heavy_metal
		self.list_rgb_color = list_rgb_color
		# Dict
		self.dict_replace_reverse = {}
		self.dict_group_K = {}  # Create an empty dict to store phylogenetic signal K
		self.dict_feature = {"SVL": "log_Concentration-Physical_feature-SVL(mm)",
		                     "BM": "log_Concentration-Physical_feature-Body_mass(g)",
		                     "CF": "log_Concentration-Physical_feature-CF(100*BM(g)/SVL(mm)³)"}
		# Axes
		self.ax_slope = None
		self.ax_intercept = None
		self.ax_pca = None
		self.ax_radar = None
		self.ax_heatmap = None
		self.ax_heatmap_legend = None
		self.ax_reg = None

		# Initialize
		print("__________________________________")
		print("Initializing.")
		# Check path
		list_dir = [self.str_working_dir,
		            self.str_output_directory,
		            self.str_dataframe_path,
		            self.str_continent_path,
		            self.str_world_path,
		            self.str_tree_path]
		for str_path in list_dir:
			if not os.path.exists(str_path):
				if str_path[-1] == "\\":
					os.mkdir(self.str_output_directory)
				else:
					raise ValueError("Invalid path: " + str_path)
		self.color_setup()
		self.build_canvas()
		self.read_excel()
		self.df2quantile(self.str_feature)

		# Main text figures
		print("__________________________________")
		print("Drawing the main text figures.")
		self.draw_slope(ax=self.ax_slope, feature=self.str_feature)
		self.draw_pca(delete_abnormal=True)
		self.draw_reg()
		self.draw_heatmap()
		self.add_serials([self.ax_slope, self.ax_pca, self.ax_heatmap])
		plt.savefig(self.str_output_directory+f"Main_text.png_{self.data_date}.png",dpi = 300)
		plt.savefig(self.str_output_directory + f"Main_text.png_{self.data_date}.svg")

		# Supplementary figures
		print("__________________________________")
		print("Drawing the supplementary figures.")
		self.draw_radar()

		self.draw_full_quantile()

		plt.show()

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
	def delete_abnormal_values(df2delete):
		"""
		Remove the abnormal values
		Fill null values with the mean values

		___________
		Parameters:
			df2delete: A dataframe with potential abnormal values

		_______
		Return:
			df: A dataframe without abnormal values
		"""

		df_copy = df2delete.copy()  # Copy the dataframe
		list_head = list(df_copy)  # Get a list of the column names

		# Select the Heavy Metal Columns
		for num, h in enumerate(list_head):  # Iterate each column name
			q1 = df_copy[h].quantile(0.25)  # Get values at the 25% quantile point
			q3 = df_copy[h].quantile(0.75)  # Get values at the 75% quantile point
			iqr = q3 - q1  # Get the data value range within 25% - 75%
			# Determine the Range to Distinguish the Abnormal Values
			num_lower_limit = q1 - 1.5 * iqr  # Upper limit
			num_upper_limit = q3 + 1.5 * iqr  # Lower limit
			num_mean = np.mean(df_copy[h])  # Calculate the mean value of the column
			# for i in df[h]:  # Iterate each value of the current column
			# 	if i > num_upper_limit or i < num_lower_limit:  # If the value is not in the range
			# 		df.replace({h: {i: num_mean}},
			# 		           inplace=True)  # Replace the abnormal values with the mean value of corresponding columns
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
			position = [0.9, 0.1]
		# List needs to be added with serial
		list_tb = [chr(i) for i in range(97, 123)]  # A list of alphabetic serials
		for n, la in enumerate(list_ax):  # Iterate each subplot
			# Add serial marks
			la.text(x=la.get_xlim()[0] * position[0] - la.get_xlim()[1] * position[1],  # X coordinates
			        y=la.get_ylim()[1],  # Y coordinates
			        s=f"{list_tb[n].upper()}",  # Text content
			        weight="bold",  # Bold text
			        ha="left",  # Horizontal alignment
			        va="top")  # Vertical alignment

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
		# Create a list of even discrete colors, used for clustered groups
		self.list_colors_cluster = self.linear_cmap(np.linspace(0, 1, self.n_clusters))
		# Create a palette, used when we used "seaborns" to plot
		self.palette = sns.color_palette(self.list_colors_cluster)

	def build_canvas(self, fig_width=12, fig_height=10, h_space=0.1, w_space=0.3):
		"""
		Build up a canvas

		___________
		Parameters:
			fig_width: The width of the figure
			fig_height: The height of the figure
			h_space: The vertical gap between two subplots
			w_space: The horizontal gap between two subplots
		"""

		fig_main = plt.figure(figsize=(fig_width, fig_height))  # Figure object
		gs = GridSpec(nrows=2,
		              ncols=4,
		              height_ratios=[1, 0.5],
		              width_ratios=[1, 1, 1, 0.05],
		              hspace=h_space,
		              wspace=w_space)
		self.ax_slope = fig_main.add_subplot(gs[0, 0])  # Create an axis to plot the quantile slope
		self.ax_pca = fig_main.add_subplot(gs[0, 1:])  # Create an axis to plot the PCA
		self.ax_heatmap = fig_main.add_subplot(gs[1, :3])  # Create an axis to plot the heatmap
		self.ax_heatmap_legend = fig_main.add_subplot(gs[1, 3])  # Create an axis to plot the heatmap legend

	def read_excel(self):
		"""
		Read a dataframe
		"""

		# Read a dataframe
		self.df = pd.read_excel(self.str_dataframe_path, header=[0])
		# Extract rows, indicating Anurans with Gosner Stage larger than 46 (that is, metamorphosed Anurans)
		self.df = self.df[self.df["Concentration-Gosner_Stage"] > 46]
		# Extract whole-body HM bio-accumulation columns
		self.list_head = [x for x in list(self.df) if "Whole_Body" in x and x.split("-")[-1] in self.list_hm]

	def df2quantile(self, feature="CF"):
		"""
		Extract a dataframe according to the "feature"

		___________
		Parameters:
			feature: Body size feature. Choose from "SVL", "BM" and "CF"
		"""

		self.df_quantile = self.df[["Continent",
		                            "Species",
		                            self.dict_feature[feature]] + self.list_head].copy()  # Extract columns
		self.df_quantile.columns = ["Continent",
		                            "Species",
		                            feature] + self.list_hm  # Rename column names
		self.df_quantile.loc[:, "Species"] = self.df_quantile["Species"].str.replace(" ",
		                                                                             "_")  # Add a column of "Species"

	def draw_scatter(self, ax, feature="CF"):
		"""
		Draw quantile regression scatters

		___________
		Parameters:
			ax: Axes
			feature: Body size feature. Choose from "SVL", "BM" and "CF"
		"""

		for m, (hm, color) in enumerate(
				zip(self.list_hm, self.list_colors)):  # Iterate each HM and its corresponding color
			df_temp = self.df_quantile[[hm, "Species", feature]]  # Select required columns

			if df_temp.shape[0] > 0:  # If there are at least one row in "df_temp"
				ax.scatter(df_temp[hm],  # X coordinates
				           df_temp[feature],  # Y coordinates
				           label=hm,  # Mark with the name "hm"
				           alpha=0.05,  # Alpha
				           color=color,  # Scatter color
				           zorder=0)  # Set as the bottommost layer
				# Get phylogenetic quantile regression results
				models = phylo_quantreg(str_Y=feature,  # The column name of "df_temp", set as the X axis
				                        str_X=hm,  # The column name of "df_temp", set as the Y axis
				                        tree=self.tree,  # A phylogenetic tree
				                        df=df_temp,  # Data
				                        list_quantile=self.list_quantile)  # A list of quantiles
				models = pd.DataFrame(models, columns=["quantile", "intercept", "slope", "lb",
				                                       "ub"])  # Build a dataframe of model results
				x = np.arange(df_temp[hm].min(), df_temp[hm].max(),
				              0.1)  # Create a list of x coordinates within the range of HM values
				get_y = lambda a, b: a + b * x  # Create a formula to get y
				for ii in range(models.shape[0]):  # Iterate each quantile point
					y = get_y(models.intercept[ii], models.slope[ii])  # Get a list of y coordinates
					ax.plot(x,  # X coordinates
					        y,  # Y coordinates
					        linestyle="-",  # Line style
					        color=color,  # Line color
					        alpha=0.4)  # Alpha

		# Adjust axes and tick labels
		ax.set_xlabel("lg $\\it{HM}$")  # Set labels on the X axes
		ax.set_ylabel(feature)  # Set labels on the Y axes
		ax.yaxis.set_major_formatter(FormatStrFormatter("%1.3f"))  # Decimals

		print("Successfully plot scatters")

	def draw_slope(self, ax, feature = "CF", phylo_on = True):
		"""
		Draw quantile regression slopes

		___________
		Parameters:
			ax: Axes
			feature: Body size feature. Choose from "SVL", "BM" and "CF"
		"""

		for m, (hm, color) in enumerate(
				zip(self.list_hm, self.list_colors)):  # Iterate each HM and its corresponding color
			df_temp = self.df_quantile[[hm, "Species", feature]]  # Select required columns

			if df_temp.shape[0] > 0:  # If there are at least one row in "df_temp"
				# Get phylogenetic quantile regression results
				models = phylo_quantreg(str_X=feature,  # The column name of "df_temp", set as the X axis
				                        str_Y=hm,  # The column name of "df_temp", set as the Y axis
				                        tree=self.tree,  # A phylogenetic tree
				                        df=df_temp,  # Data
				                        list_quantile=self.list_quantile,  # A list of quantiles
				                        phylo_on=phylo_on)
				models = pd.DataFrame(models, columns=["quantile", "intercept", "slope", "lb",
				                                       "ub"])  # Build a dataframe of model results
				ax.plot(models["quantile"],  # X coordinates
				        models["slope"],  # Y coordinates
				        label=hm,  # Mark with the name "hm"
				        color=color)  # Line color
				ax.fill_between(models["quantile"],  # X coordinates
				                models["lb"],  # Y coordinates of lower boundary
				                models["ub"],  # Y coordinates of upper boundary
				                color=color,  # Line color
				                alpha=0.4)  # Alpha

		# Adjust axes and tick labels
		ax.legend()  # Show legends
		ax.set_xlabel("Quantile")  # Set labels on X axis
		ax.set_ylabel("Slope")  # Set labels on Y axis
		ax.xaxis.set_major_locator(
			MultipleLocator(0.25))  # Set the major tick locator for the x-axis with 0.25 interval
		ax.yaxis.set_major_formatter(FormatStrFormatter("%1.3f"))  # Decimals

		# Add vertical lines
		list_percentile = [0.25, 0.50, 0.75]  # X coordinates of the vertical lines
		for p in list_percentile:  # Iterate each position
			ax.axvline(p,  # X coordinates
			           color="red",  # Line colors
			           linewidth=0.7,  # Line width
			           linestyle="--",  # Line style
			           alpha=0.7)  # Line alpha

		print("Successfully plot slope")

	def draw_intercept(self, feature="CF"):
		"""
		Draw quantile intercept

		___________
		Parameters:
			feature: Body size feature. Choose from "SVL", "BM" and "CF"
		"""

		for m, (hm, color) in enumerate(
				zip(self.list_hm, self.list_colors)):  # Iterate each HM and its corresponding color
			df_temp = self.df_quantile[hm].copy()  # Select required columns
			df_temp = df_temp[["Concentration", feature]]
			df_temp.dropna(how="any", axis=0, inplace=True)

			if df_temp.shape[0] > 0:
				mod = smf.quantreg(f"{feature} ~ Concentration", df_temp)
				models = [self.fit_model(X, mod, "Concentration") for X in self.list_quantile]
				models = pd.DataFrame(models, columns=["quantile", "intercept", "slope", "lb", "ub"])
				self.ax_intercept.plot(models["quantile"], models["intercept"], label=hm, color=color)

		self.ax_intercept.set_ylabel("Intercept")
		self.ax_slope.xaxis.set_major_locator(MultipleLocator(0.25))
		self.ax_intercept.yaxis.set_major_formatter(FormatStrFormatter("%1.3f"))

		print("Successfully plot intercept")

	def draw_pca(self, list_head_original=None, list_head_numeric=None, delete_abnormal=True):
		"""
		Draw PCA

		___________
		Parameters:
			list_head_original: A list of column names, used to extract from a dataframe
			list_head_numeric:  A list of numeric columns, used to rename the column names
			delete_abnormal:    Delete the abonormal values or not
		"""

		# Used to decide the cluster number
		# def WGSS(df):
		# 	wgss = []
		# 	for i in range(12):
		# 		clustering = KMeans(n_clusters=i + 1, random_state=0).fit(df)
		# 		wgss.append(clustering.inertia_)
		#
		# 	plt.plot([i + 1 for i in range(12)], wgss, marker="o")
		# 	plt.show()
		#
		# def silhouette(df):
		# 	silhouette_scores = []
		# 	for i in range(1, 24):
		# 		cluster = KMeans(n_clusters=i + 1, random_state=0).fit(df)
		# 		y_pred = cluster.labels_    # Get clustering results
		# 		silhouette_avg = silhouette_score(df, y_pred)   # Calculate silhouette
		# 		silhouette_scores.append(silhouette_avg)
		# 	plt.plot([i + 1 for i in range(1, 24)], silhouette_scores, marker="d")  # Draw "silhouette_scores"
		# 	plt.show()

		def confidence_ellipse(x,
		                       y,
		                       ax,
		                       n_std=2.0,
		                       facecolor="none",
		                       **kwargs):
			cov = np.cov(x, y)
			pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
			ell_radius_x = np.sqrt(1 + pearson)
			ell_radius_y = np.sqrt(1 - pearson)
			ellipse = Ellipse(xy=(0, 0),
			                  width=ell_radius_x * 2,
			                  height=ell_radius_y * 2,
			                  facecolor=facecolor,
			                  **kwargs)

			scale_x = np.sqrt(cov[0, 0]) * n_std
			scale_y = np.sqrt(cov[1, 1]) * n_std

			transform = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(np.mean(x), np.mean(y))

			ellipse.set_transform(transform + ax.transData)
			return ax.add_patch(ellipse)

		# Initialize
		if list_head_original is None:
			list_head_original = ["log_Concentration-Physical_feature-SVL(mm)",
			                      "log_Concentration-Physical_feature-Body_mass(g)",
			                      "log_Concentration-Physical_feature-CF(100*BM(g)/SVL(mm)³)",
			                      "Concentration-Physical_feature-ITV_SVL",
			                      "Concentration-Physical_feature-ITV_BM",
			                      "Concentration-Physical_feature-ITV_CF",
			                      "log_Concentration-Whole_Body-Cd",
			                      "log_Concentration-Whole_Body-Cr",
			                      "log_Concentration-Whole_Body-Cu",
			                      "log_Concentration-Whole_Body-Fe",
			                      "log_Concentration-Whole_Body-Mn",
			                      "log_Concentration-Whole_Body-Pb",
			                      "log_Concentration-Whole_Body-Zn"]

		if list_head_numeric is None:
			list_head_numeric = ["SVL",
			                     "BM",
			                     "CF",
			                     "$\\it{ITV}$$_{SVL}$",
			                     "$\\it{ITV}$$_{BM}$",
			                     "$\\it{ITV}$$_{CF}$",
			                     "lg $\\it{Cd}$",
			                     "lg $\\it{Cr}$",
			                     "lg $\\it{Cu}$",
			                     "lg $\\it{Fe}$",
			                     "lg $\\it{Mn}$",
			                     "lg $\\it{Pb}$",
			                     "lg $\\it{Zn}$"]

		self.df_pca = self.df[list_head_original].copy()  # Select required columns
		self.df_pca.columns = list_head_numeric  # Rename the column names

		if delete_abnormal:  # If we need to delete the abnormal values
			self.df_pca = self.delete_abnormal_values(self.df_pca)  # Delete the abnormal values

		# WGSS(self.df_pca)
		# silhouette(self.df_pca)

		# Standardize the data
		scaler = StandardScaler()
		X_scaled = scaler.fit_transform(self.df_pca)

		# Cluster and remark the cluster name
		kmeans = KMeans(n_clusters=self.n_clusters, random_state=42)  # KMeans object
		self.df_pca.loc[:, "Cluster_KMeans"] = kmeans.fit_predict(
			X_scaled)  # Cluster according to the features of "df_pca". Store the clustered results in column "Cluster_KMeans"
		self.df_pca.sort_values(by=["Cluster_KMeans"], ascending=True,
		                        inplace=True)  # Sort according to the cluster values
		list_tb = [chr(i) for i in range(97, 123)]  # A list of alphabet
		dict_group = {i: f"Group {list_tb[i].upper()}" for i in
		              range(self.n_clusters)}  # Create a dictionary used to rename the clustered groups
		self.df_pca.replace({"Cluster_KMeans": dict_group},
		                    inplace=True)  # Update the cluster name according to the "dict_group"
		self.df_pca.replace(self.dict_replace_reverse).to_excel(self.str_working_dir + "Kmeans.xlsx",
		                                                        index=True)  # Save "df_pca"
		self.dict_replace_reverse[
			"Cluster_KMeans"] = dict_group  # Update the cluster name according to the "dict_group"

		# Dimension reduction to 2D with PCA
		pca = PCA(n_components=4)
		X_pca = pca.fit_transform(X_scaled)

		# Add scatters in PCA
		sns.scatterplot(x=X_pca[:, 0],  # X coordinates
		                y=X_pca[:, 1],  # Y coordinates
		                hue=self.df_pca["Cluster_KMeans"],  # Groups
		                linewidth=0.02,  # Scatter edge width
		                s=40,  # Scatter size
		                palette=self.palette,  # Color scheme
		                alpha=0.3,  # Alpha
		                ax=self.ax_pca)  # Axes

		# Add confidence interval ellipses for each cluster
		for cluster, color in zip(self.df_pca["Cluster_KMeans"].unique(),
		                          self.list_colors_cluster):  # Iterate each cluster
			mask = self.df_pca[
				       "Cluster_KMeans"] == cluster  # Get the position of the corresponding group "Cluster_KMeans"
			# Add confidence ellipses
			confidence_ellipse(X_pca[mask, 0],
			                   X_pca[mask, 1],
			                   ax=self.ax_pca,  # Axes
			                   n_std=2,
			                   linewidth=2,  # Edge width
			                   edgecolor=color,  # Edge color of the ellipse
			                   facecolor=color,  # Face color of the ellipse
			                   alpha=0.15)  # Alpha

		# Add arrows, starting from (0, 0)
		arrow_scale = 3.0  # Scale factor to adjust the size of the arrows
		feature_vectors = pca.components_.T * np.sqrt(pca.explained_variance_)  # The ends the arrows point to
		for i, col in enumerate(list_head_numeric):  # Iterate each feature
			# Add arrows
			self.ax_pca.arrow(0,  # X coordinate of the start of the arrow
			                  0,  # Y coordinate of the start of the arrow
			                  feature_vectors[i, 0] * arrow_scale,  # X deviation
			                  feature_vectors[i, 1] * arrow_scale,  # Y deviation
			                  width=0.01,  # Width of the arrow tail
			                  head_width=0.1,  # Width of the head tail
			                  alpha=0.8,  # Arrow alpha
			                  color="red")  # Arrow color

			# Add annotation for the current arrow
			self.ax_pca.text(x=feature_vectors[i, 0] * arrow_scale * 1.15,  # X coordinates
			                 y=feature_vectors[i, 1] * arrow_scale * 1.15,  # Y coordinates
			                 s=col,  # Text
			                 fontsize=10,  # Font size
			                 ha="center",  # Horizontal alignment
			                 va="center",  # Vertical alignment
			                 color="red")  # Text color

		# Adjust the axes
		self.ax_pca.spines["left"].set_position("zero")  # Move the left axis as the Y axis
		self.ax_pca.spines["bottom"].set_position("zero")  # Move the bottom axis as the X axis
		self.ax_pca.spines["right"].set_visible(False)  # Do not show the right axis
		self.ax_pca.spines["top"].set_visible(False)  # Do not show the top axis
		self.ax_pca.grid(True, linestyle="--", alpha=0.6)  # Add a background grid

		# Add legend and text
		# Add labels for X axis
		self.ax_pca.text(x=self.ax_pca.get_xlim()[1],  # X coordinates
		                 y=abs(self.ax_pca.get_ylim()[1] - self.ax_pca.get_ylim()[0]) / 2 * 0.1,  # Y coordinates
		                 s="PC 1 (%.1f%%)" % (pca.explained_variance_ratio_[0] * 100),  # Text
		                 ha="right",  # Horizontal alignment
		                 va="top")  # Vertical alignment
		# Add labels for Y axis
		self.ax_pca.text(x=abs(self.ax_pca.get_xlim()[1] - self.ax_pca.get_xlim()[0]) / 2 * 0.05,  # X coordinates
		                 y=self.ax_pca.get_ylim()[0],  # Y coordinates
		                 s="PC 2 (%.1f%%)" % (pca.explained_variance_ratio_[1] * 100),  # Text
		                 rotation=90,  # Rotate the text
		                 ha="left",  # Horizontal alignment
		                 va="bottom")  # Vertical alignment
		self.ax_pca.legend_.set_title("")  # Do not show the title of the legend
		self.ax_pca.legend_.set_frame_on(False)  # Do not show the frame of the legend

		print("Successful PCA plot!")

	def draw_heatmap(self):
		"""
		Draw a heatmap, showing the correlation between the HMs and body size
		"""

		df_heatmap = self.df_pca.copy()  # Copy the dataframe
		df_heatmap.drop("Cluster_KMeans", axis=1, inplace=True)  # Remove the column "Cluster_KMeans"
		list_head = list(df_heatmap)  # Get a list of the head

		df_cor = df_heatmap.corr()  # Get a correlation matrix
		list_hm_formatted = ["lg $\\it{Cd}$",
		                     "lg $\\it{Cr}$",
		                     "lg $\\it{Cu}$",
		                     "lg $\\it{Fe}$",
		                     "lg $\\it{Mn}$",
		                     "lg $\\it{Pb}$",
		                     "lg $\\it{Zn}$"]  # A list of formatted HMs

		# Screen to get a 7 * 9 dataframe
		# Seven rows for studied HMs
		# Nine columns for body size features
		for head in list_head:  # Iterate "list_head"
			if head in list_hm_formatted:  # If the head is one of the HMs
				df_cor.drop(head, axis=1, inplace=True)  # Drop the current column
		df_cor = df_cor.loc[list_hm_formatted, :]  # Extract rows involving the HMs

		# Draw heatmap of correlation matrix
		sns.heatmap(df_cor,  # The correlation matrix
		            fmt=".3f",  # Decimals
		            annot=False,  # Do not annotate
		            ax=self.ax_heatmap,  # Axes
		            cbar_ax=self.ax_heatmap_legend,  # Color bar
		            cmap=self.linear_cmap,  # Color scheme
		            cbar_kws={"aspect": 30})  # Width / Height of the color bar

		# Annotate the values in cells
		for i in range(df_cor.shape[1]):  # Iterate each column
			for j in range(df_cor.shape[0]):  # Iterate each row
				self.ax_heatmap.text(x=i + 0.5,  # X coordinates
				                     y=j + 0.5,  # Y coordinates
				                     s=f"{df_cor.iloc[j, i]:.3f}",  # Annotation content
				                     ha="center",  # Horizontal alignment
				                     va="center",  # Vertical alignment
				                     color="white")  # Text color
		self.ax_heatmap.set_yticklabels(list_hm_formatted, rotation=0, ha="right")  # Set tick labels on y axis

		print("Successful heatmap plot!")

	def draw_radar(self, norm_method="minmax"):
		"""
		Draw radar plot

		___________
		Parameters:
			norm_method: Method of standardizing the data
		"""

		f = plt.figure(figsize=(12, 6),
		               constrained_layout=True)
		gs_temp = GridSpec(nrows=1,
		                   ncols=1,
		                   height_ratios=[1],
		                   width_ratios=[1],
		                   hspace=0.3,
		                   wspace=0.1)
		self.ax_radar = f.add_subplot(gs_temp[0, 0], polar=True)

		r = 1.1  # Radius of the plot
		self.df_radar = self.df_pca.copy()  # Copy the dataframe
		self.df_radar = self.df_radar.join(self.df.loc[self.df_pca.index, "Latitude"],  # Concat corresponding latitudes
		                                   lsuffix="_caller",
		                                   rsuffix="_other",
		                                   how="inner")
		self.df_radar.rename(columns={"Latitude": "Lat"}, inplace=True)  # Rename the index "Latitude"
		list_head = list(self.df_radar)  # Get a list of the heads of "self.df_radar"
		list_head.remove("Cluster_KMeans")  # Remove "Cluster_KMeans" from the list
		# Standardize
		if norm_method == "minmax":  # If use "min-max" to standardize
			for head in list_head:  # Iterate each feature
				self.df_radar[head] = (self.df_radar[head] - self.df_radar[head].min()) / (
						self.df_radar[head].max() - self.df_radar[head].min())  # Standardize
		elif norm_method == "zscore":  # If use "z-score" to standardize
			self.df_radar = (self.df_radar - self.df_radar.mean()) / (self.df_radar.std())

		df_pivot = self.df_radar.pivot_table(index=["Cluster_KMeans"], aggfunc="mean")  # Pivot the dataframe
		df_pivot.rename(index=self.dict_replace_reverse["Cluster_KMeans"],
		                inplace=True)  # Update the name of each clustered group
		df_pivot = df_pivot.reindex(columns=list_head)  # Replace the index with "list_head"
		df_pivot.to_excel(self.str_output_directory + f"Radar_{self.data_date}.xlsx",
		                  index=False)  # Save dataframe of radar plot

		# Build theta list
		angles = np.linspace(0, 2 * np.pi, df_pivot.shape[1], endpoint=False)
		angles = np.concatenate((angles, [angles[0]]))

		# Draw background grid
		self.ax_radar.plot(angles,
		                   [r] * len(angles),
		                   linewidth=0.5,
		                   color="gray",
		                   alpha=0.4)

		# Draw
		for cluster, color in zip(df_pivot.index,
		                          self.list_colors_cluster):  # Iterate each group and its corresponding color
			list_cluster = list(df_pivot.loc[cluster, :])
			values = np.concatenate((list_cluster, [list_cluster[0]]))

			self.ax_radar.plot(angles,  # A theta list
			                   values,  # A radius list
			                   "o-",  # Vertice style
			                   linewidth=2,  # Line width
			                   color=color)  # Color
			self.ax_radar.fill(angles,  # A theta list
			                   values,  # A radius list
			                   alpha=0.25,  # Alpha
			                   color=color)  # Filling color

		# Add grids and labels
		for angle, label in zip(angles[:-1], list(df_pivot)):  # Remove the last angle to avoid duplication
			# Plot radial lines
			self.ax_radar.plot([angle, angle],  # A theta list
			                   [0, r],  # A radius list
			                   "gray",  # Color
			                   linewidth=0.5,  # Line width
			                   zorder=0,  # Put the layer order to the bottommost
			                   alpha=0.4)  # Alpha
			# Add text
			self.ax_radar.text(x=angle,
			                   y=r * 1.1,
			                   s=label,
			                   ha="center",
			                   va="center",
			                   fontsize=10)
		self.ax_radar.legend(loc="right")
		self.ax_radar.axis("off")
		self.ax_radar.grid(True)
		f.savefig(self.str_output_directory + "Radar.png", dpi=300)
		f.savefig(self.str_output_directory + "Radar.svg")
		plt.close()

		print("Successful radar plot!")

	def draw_dist(self):
		"""
		Draw feature distribution of different groups
		"""

		df = pd.read_excel(self.str_working_dir + "Kmeans_Location.xlsx", sheet_name="Sheet1")[["Cluster_KMeans", "CF"]]
		sns.histplot(data=df,
		             x="CF",
		             hue="Cluster_KMeans",
		             kde=True,
		             bins=200,
		             ax=self.ax_dist,
		             palette=self.palette,
		             stat="density",
		             line_kws={"color": "white"})

		list_percentile = [25, 50, 75]
		quantiles_temp = np.percentile(df["CF"], list_percentile)
		# 绘制分位数竖线
		for (q, p) in zip(quantiles_temp, list_percentile):
			self.ax_dist.axvline(q,
			                     color="red",
			                     linewidth=0.7,
			                     linestyle="--",
			                     alpha=0.7)
			self.ax_dist.text(q,
			                  1.9,
			                  f"{p}" + "$^{th}$",
			                  rotation=90,
			                  va="bottom",
			                  ha="right",
			                  fontsize=9,
			                  color="red")
		self.ax_dist.legend_.set_title("")
		self.ax_dist.legend_.set_frame_on(False)
		self.ax_dist.set_xlabel("lg CF")
		self.ax_dist.spines["top"].set_visible(False)
		self.ax_dist.spines["right"].set_visible(False)

	def draw_map(self):
		self.df_loc = self.df[["Latitude", "Longitude"]]
		self.df_loc = self.df_pca.join(self.df_loc, lsuffix="_caller", rsuffix="_other", how="inner")
		self.df_loc.replace(self.dict_replace_reverse, inplace=True)
		self.df_loc.to_excel(self.working_dir + "Kmeans_Location.xlsx", index=False)
		world = gpd.read_file(self.world_path)
		world.plot(ax=self.ax_world, color="gray", edgecolor="white", linewidth=0.05, legend=False, alpha=0.2)
		list_group = list(set(self.df_loc["Cluster_KMeans"]))
		for n, lg in enumerate(list_group):
			list_lat = self.df_loc[self.df_loc["Cluster_KMeans"] == lg]["Latitude"]
			list_lon = self.df_loc[self.df_loc["Cluster_KMeans"] == lg]["Longitude"]
			self.ax_world.scatter(list_lon, list_lat, label=lg, color=self.list_colors_cluster[n], alpha=0.1, s=12)
		self.ax_world.axhline(30, linestyle="--", linewidth=0.4, color="gray")
		self.ax_world.axhline(0, linestyle="--", linewidth=0.5, color="gray")
		self.ax_world.axhline(-30, linestyle="--", linewidth=0.4, color="gray")
		# self.ax_world.set_ylim(-70,90)
		self.ax_world.axis("off")
		self.ax_world.legend()

		print("World map plotted")

	def draw_reg(self):
		"""
		Draw regression curves of different groups
		"""

		def trim_tips(tree, node):
			"""
			Trim a node in a phylogenetic tree

			___________
			Parameters:
				tree: Phylogenetic tree
				node: The node needs to be removed
			"""
			if node.is_terminal():
				parent = tree.prune(node)
				trim_tips(tree, parent)

		f = plt.figure(figsize=(12, 12),
		               constrained_layout=True)
		gs_temp = GridSpec(nrows=2,
		                   ncols=2,
		                   height_ratios=[1, 1],
		                   width_ratios=[1, 1],
		                   hspace=0.3,
		                   wspace=0.3)
		self.ax_reg = []
		for r in range(2):
			for c in range(2):
				self.ax_reg.append(f.add_subplot(gs_temp[r, c]))

		list_group = list(set(self.df_pca["Cluster_KMeans"]))  # Get a list of group marks
		list_group.sort()
		# A list of features for phylogenetic signal calculation
		list_feature = ["SVL",
		                "BM",
		                "CF",
		                "$\\it{ITV}$$_{SVL}$",
		                "$\\it{ITV}$$_{BM}$",
		                "$\\it{ITV}$$_{CF}$"]

		df_cluster = pd.concat([self.df_pca[list_feature + ["Cluster_KMeans"]], self.df["Species"]],
		                       axis=1)  # Build a dataframe

		for n, (group, color) in enumerate(zip(list_group, self.list_colors_cluster)):  # Iterate each group
			anno_text = ""  # Create an empty text to store
			df_group = df_cluster[df_cluster["Cluster_KMeans"] == group].copy()  # Select data in "group"
			df_group.drop("Cluster_KMeans", axis=1, inplace=True)  # Drop the column marking cluster groups
			df_group_pivot = df_group.pivot_table(values=list_feature,
			                                      index=["Species"],
			                                      aggfunc="mean")  # Reclassify the dataframe
			list_group_species = [i.replace(" ", "_") for i in df_group_pivot.index]  # Get a list of species in "group"
			df_group_pivot.index = list_group_species  # Rename the index of "df_group_pivot", making the Latin Names compatible between tree and dataframe
			# Build a phylogenetic tree of "group"
			group_tree = copy.deepcopy(self.tree)  # Copy a complete tree
			for sp in group_tree.get_terminals():  # Iterate each tip
				if sp.name not in list_group_species:  # If the tip "sp" is not in "group"
					trim_tips(group_tree, sp)  # Trim the irrelevant tips

			# Calculate the phylogenetic signals to show the evolutionary difference of each "group"
			self.dict_group_K[group] = {}  # Create an empty dict to store the parameters
			for feature in list_feature[:3]:  # Iterate each feature
				self.dict_group_K[group][feature] = k_multi(group_tree, df_group_pivot[
					feature])  # Calculate the phylogenetic signals K of "feature"

			results = smf.ols("BM~SVL", data=df_group_pivot).fit()  # Fit the data

			# Store regression results into "dict_para"
			dict_para = {"Slope": results.params["SVL"],  # Slope of the regression curves
			             "Intercept": results.params["Intercept"]}  # Intercept of the regression curves

			# Organize an annotation text
			anno_text += f"a = {dict_para["Slope"]:.3f}\nb = {dict_para["Intercept"]:.3f}\n"
			itv_text = ""
			for feature, K in self.dict_group_K[group].items():  # Iterate each feature to get phylogenetic signals K
				self.df_pca.loc[df_group.index, "$\\it{K}$$_{" + feature + "}$"] = np.sum(K)
				anno_text += "$\\it{K}$$_{" + feature + "}" + f"$: {np.sum(K):.3f}\n"  # Append parameters of phylogenetic signals K
				itv_text += ("$\\overline{ITV}$" + "$_{" + feature + "}$:"
				             + f" {np.mean(df_group_pivot["$\\it{ITV}$$_{" + feature + "}$"]):.3f} ± {np.std(df_group_pivot["$\\it{ITV}$$_{" + feature + "}$"]):.3f}\n")
			anno_text += itv_text

			# Construct x and y values to plot the regression results
			X = np.linspace(min(df_group["SVL"]), max(df_group["BM"]), 50)  # X range
			Y = dict_para["Intercept"] + dict_para["Slope"] * X  # Y values of regression curves

			sns.scatterplot(x=df_group["SVL"],  # X coordinates
			                y=df_group["BM"],  # Y coordinates
			                linewidth=0.02,  # Edge width of scatters
			                s=40,  # Size of the scatters
			                color=self.linear_cmap((n + 1) / (self.n_clusters - 1)),  # Color schemes
			                alpha=0.15,  # Alpha of scatters
			                ax=self.ax_reg[n])  # Axes
			self.ax_reg[n].set_xlim(1.4, 2.3)
			self.ax_reg[n].set_ylim(0, 3)
			# Draw regression results
			self.ax_reg[n].scatter(df_group["SVL"],  # X coordinates
			                       df_group["BM"],  # Y coordinates
			                       color=color,  # Color of the scatters
			                       s=30,  # Size of the scatters
			                       alpha=0.1)  # Alpha
			# Plot regression curves
			self.ax_reg[n].plot(X,
			                    Y,  # Y values
			                    color=color,  # Regression curve colors
			                    label=group)  # Cluster labels
			# Add regression parameters
			self.ax_reg[n].text(x=self.ax_reg[n].get_xlim()[0] * 0.97 + self.ax_reg[n].get_xlim()[1] * 0.03,
			                    # X coordinates of annotations
			                    y=self.ax_reg[n].get_ylim()[1] * 0.97 + self.ax_reg[n].get_ylim()[0] * 0.03,
			                    # Y coordinates of annotations
			                    s=anno_text,  # Regression parameters
			                    ha="left",  # Horizontal alignment
			                    va="top",  # Vertical alignment
			                    fontsize=10)  # Font size
			self.ax_reg[n].set_xlabel("lg $\\it{SVL}$")
			self.ax_reg[n].set_ylabel("lg $\\it{BM}$")
			self.ax_reg[n].legend(loc="upper right")
		self.add_serials(self.ax_reg, position=[0.97, 0.03])
		f.savefig(self.str_output_directory + f"OLS_{self.data_date}.png", dpi=300)
		f.savefig(self.str_output_directory + f"OLS_{self.data_date}.svg")
		plt.close()
		print("Successful regression plot!")

	def draw_full_quantile(self):
		"""
		Draw scatter and quantile slope of "SVL", "BM" and "CF"
		"""

		list_feature = ["SVL", "BM", "CF"]
		dict_quantile = {"Scatter": lambda ax, f: self.draw_scatter(ax=ax, feature=f),
	                     "Slope_phylo": lambda ax, f: self.draw_slope(ax=ax, feature=f),
			             "Slope": lambda ax, f: self.draw_slope(ax=ax, feature=f,phylo_on=False)}

		fig_temp = plt.figure(figsize=(12, 12), constrained_layout=True)
		gs_temp = GridSpec(nrows=3,
		                   ncols=3,
		                   height_ratios=[1, 1, 1],
		                   width_ratios=[1, 1, 1],
		                   hspace=0.3,
		                   wspace=0.3)
		ax_quantile = [fig_temp.add_subplot(gs_temp[r, c]) for c in range(3) for r in range(3)]

		n = 0
		for fea in list_feature:
			self.df2quantile(fea)
			for quan in dict_quantile.keys():
				dict_quantile[quan](ax_quantile[n], fea)
				n += 1
		self.add_serials(ax_quantile)

		fig_temp.savefig(self.str_output_directory + f"Quantile_{self.data_date}.png", dpi=300)
		fig_temp.savefig(self.str_output_directory + f"Quantile_{self.data_date}.svg")
		plt.close()

		print("Successful quantile plot!")

	def draw_anova(self):
		fig_temp = plt.figure(figsize=(12, 12), constrained_layout=True)
		gs_temp = GridSpec(nrows=2,
		                   ncols=3,
		                   height_ratios=[1, 1],
		                   width_ratios=[1, 1, 1],
		                   hspace=0.3,
		                   wspace=0.3)
		ax_anova = [fig_temp.add_subplot(gs_temp[r, c]) for c in range(3) for r in range(2)]



		# 进行单因素方差分析（One-way ANOVA）
		f_statistic, p_value = stats.f_oneway(fertilizer_A, fertilizer_B, fertilizer_C)

		# 打印结果
		print("F-statistic:", f_statistic)
		print("P-value:", p_value)

		# 判断是否拒绝零假设
		alpha = 0.05
		if  p_value < alpha:
			    print("拒绝零假设：不同肥料对植物生长的影响有显著差异")
		else:
			    print("无法拒绝零假设：不同肥料对植物生长的影响没有显著差异")

		# 绘制箱线图
		plt.boxplot([fertilizer_A, fertilizer_B, fertilizer_C], labels=["肥料 A", "肥料 B", "肥料 C"])
		plt.title("不同肥料对植物生长的影响")
		plt.ylabel("植物生长高度 (cm)")
		plt.show()


if __name__ == "__main__":
	plt.rcParams["legend.fontsize"] = 10
	plt.rcParams["font.sans-serif"] = "Arial"
	warnings.filterwarnings("ignore", category=FutureWarning)

	date = 20250716
	str_working_dir = "E:\\学习\\研二上\\重金属\\"
	tree_path_all = str_working_dir + "0.0_Fundamental_data\\Species_List\\Species_all20250512.nwk"

	quantiles = np.arange(0.05, 0.95, 0.01)
	list_hm = ["Cd", "Cr", "Cu", "Fe", "Mn", "Pb", "Zn"]

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


	Reg(str_working_directory=str_working_dir, # Working directory
	         str_tree_path=tree_path_all,           # Tree path
	         list_heavy_metal=list_hm,              # A list of heavy metals
	         list_quantile=quantiles,               # A list of quantile points
	         list_rgb_color=list_color,             # A list of RGB colors
	         n_clusters=4,                          # The number of clustered groups
	         data_date=date,                        # The date that the preprocessed data was generated
	         feature="CF")                          # The topic body size feature

import math
import os
import random
import warnings
import numpy as np
import pandas as pd
from Bio import Phylo
from tqdm import tqdm
from datetime import datetime
from matplotlib import pyplot as plt
from fancyimpute import KNN, SoftImpute
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression

warnings.simplefilter(action='ignore', category=FutureWarning)

class Preprocessing:
	def __init__(self,
	             str_dir,
	             str_df_path,
	             str_sheetname,
	             fill_method,
	             num_iter,
	             delete_abnormal_data=True,
	             num_date=None):
		# Num
		self.num_moisture_adult = 77.80
		self.num_moisture_metamorphosis = 87.03
		self.num_moisture_larvae = 90.12
		self.num_moisture_hatchling = 93.21
		self.num_moisture_embryos = 96.30
		self.num_iter = num_iter
		self.num_date = num_date
		# Dict
		self.dict_moisture = {"Embryos": self.num_moisture_embryos, "Hatchlings": self.num_moisture_hatchling,
		                      "Larvae": self.num_moisture_larvae, "Metamorphosis": self.num_moisture_metamorphosis,
		                      "Adult": self.num_moisture_adult}
		self.dict_parts = {"Whole_Body": [], "Skin": [], "Muscle": [], "Liver": [], "Kidney": [], "Lung": [],
		                   "Gonad": [], "Blood": [], "Toe_Clip": [], "Intestine": [], "Bone": []}
		self.dict_Stage_divide = {"Embryos": [0, 19], "Hatchlings": [19, 25], "Larvae": [25, 41],
		                          "Metamorphosis": [41, 46], "Adult": [47, np.inf]}
		self.dict_pseudo = {}
		self.dict_pseudo_backward = {}
		# Str
		self.str_fill_method = fill_method
		self.str_dir = str_dir
		self.str_working_directory = self.str_dir +f"2.2_Imputation\\Imputation{formatted_date}\\"
		self.str_df_path = str_df_path
		self.str_tree_path = self.str_dir+"0.0_Fundamental_data\\Species_List\\Species_all20250420.nwk"
		self.str_iteration_path = self.str_working_directory + f"01{self.str_fill_method}_impute\\"
		self.str_sheetname = str_sheetname
		# List
		self.list_2pseudo = ['Country', 'Continent', 'Species', 'Concentration-Sex']
		# Dataframe
		self.df = pd.DataFrame(data=None)
		self.df_impute_linear = pd.DataFrame(data=None)
		self.df_impute_knn = None
		self.df_impute_soft = pd.DataFrame(data=None)
		self.df_impute_tree = pd.DataFrame(data=None)
		self.df_weight = pd.DataFrame(data=None)
		# Initialize
		list_path = [self.str_dir,
					 self.str_working_directory,
					 self.str_df_path,
					 self.str_tree_path,
		             self.str_iteration_path]
		for lp in list_path:
			if not os.path.exists(lp):
				if lp[-1] == "\\":
					os.mkdir(lp)
				else:
					raise ValueError(f"Invalid path: {lp}")

		self.df = self.read_df()
		self.df_head = list(self.df)
		self.rename_columns()
		self.get_slices()
		self.execute_ww2dw()
		if delete_abnormal_data:
			self.delete_abnormal_values()
		self.log_convert(10)
		self.pseudo_data()
		self.fill_blanks()
		self.mean_interpolation()
		self.linear_modeling_revise(self.df_impute_tree,
		                            independent_value=["Latitude",
		                                               "Climate-Actual_Evapotranspiration(mm)_mean",
		                                               "Climate-Actual_Evapotranspiration(mm)_min",
		                                               "Climate-Actual_Evapotranspiration(mm)_max",
		                                               "Climate-Actual_Evapotranspiration(mm)_median",
		                                               "Climate-Actual_Evapotranspiration(mm)_stdDev",
		                                               "Climate-Climate_Water_Deficit(mm)_mean",
		                                               "Climate-Climate_Water_Deficit(mm)_min",
		                                               "Climate-Climate_Water_Deficit(mm)_max",
		                                               "Climate-Climate_Water_Deficit(mm)_median",
		                                               "Climate-Climate_Water_Deficit(mm)_stdDev",
		                                               "Climate-Palmer_Drought_Severity_Index_mean",
		                                               "Climate-Palmer_Drought_Severity_Index_min",
		                                               "Climate-Palmer_Drought_Severity_Index_max",
		                                               "Climate-Palmer_Drought_Severity_Index_median",
		                                               "Climate-Palmer_Drought_Severity_Index_stdDev",
		                                               "Climate-Reference_Evapotranspiration(mm)_mean",
		                                               "Climate-Reference_Evapotranspiration(mm)_min",
		                                               "Climate-Reference_Evapotranspiration(mm)_max",
		                                               "Climate-Reference_Evapotranspiration(mm)_median",
		                                               "Climate-Reference_Evapotranspiration(mm)_stdDev",
		                                               "Climate-Precipitation_Accumulation(mm)_mean",
		                                               "Climate-Precipitation_Accumulation(mm)_min",
		                                               "Climate-Precipitation_Accumulation(mm)_max",
		                                               "Climate-Precipitation_Accumulation(mm)_median",
		                                               "Climate-Precipitation_Accumulation(mm)_stdDev",
		                                               "Climate-Runoff(mm)_mean",
		                                               "Climate-Runoff(mm)_min",
		                                               "Climate-Runoff(mm)_max",
		                                               "Climate-Runoff(mm)_median",
		                                               "Climate-Runoff(mm)_stdDev",
		                                               "Climate-Soil_Moisture(mm)_mean",
		                                               "Climate-Soil_Moisture(mm)_min",
		                                               "Climate-Soil_Moisture(mm)_max",
		                                               "Climate-Soil_Moisture(mm)_median",
		                                               "Climate-Soil_Moisture(mm)_stdDev",
		                                               "Climate-Downward_Surface_Shortwave_Radiation(W/m²)_mean",
							       "Climate-Downward_Surface_Shortwave_Radiation(W/m²)_min",
							       "Climate-Downward_Surface_Shortwave_Radiation(W/m²)_max",
							       "Climate-Downward_Surface_Shortwave_Radiation(W/m²)_median",
							       "Climate-Downward_Surface_Shortwave_Radiation(W/m²)_stdDev",
							       "Climate-Snow_Water_Equivalent(mm)_mean",
							       "Climate-Snow_Water_Equivalent(mm)_min",
							       "Climate-Snow_Water_Equivalent(mm)_max",
							       "Climate-Snow_Water_Equivalent(mm)_median",
							       "Climate-Snow_Water_Equivalent(mm)_stdDev",
							       "Climate-Minimum_Temperature(°C)_mean",
							       "Climate-Minimum_Temperature(°C)_min",
							       "Climate-Minimum_Temperature(°C)_max",
							       "Climate-Minimum_Temperature(°C)_median",
							       "Climate-Minimum_Temperature(°C)_stdDev",
							       "Climate-Maximum_Temperature(°C)_mean",
							       "Climate-Maximum_Temperature(°C)_min",
							       "Climate-Maximum_Temperature(°C)_max",
							       "Climate-Maximum_Temperature(°C)_median",
							       "Climate-Maximum_Temperature(°C)_stdDev",
							       "Climate-Vapor_Pressure(kPa)_mean",
							       "Climate-Vapor_Pressure(kPa)_min",
							       "Climate-Vapor_Pressure(kPa)_max",
							       "Climate-Vapor_Pressure(kPa)_median",
							       "Climate-Vapor_Pressure(kPa)_stdDev",
							       "Climate-Vapor_Pressure_Deficit(kPa)_mean",
							       "Climate-Vapor_Pressure_Deficit(kPa)_min",
							       "Climate-Vapor_Pressure_Deficit(kPa)_max",
							       "Climate-Vapor_Pressure_Deficit(kPa)_median",
							       "Climate-Vapor_Pressure_Deficit(kPa)_stdDev",
							       "Climate-Wind_Speed_At_10m(m/s)_mean",
							       "Climate-Wind_Speed_At_10m(m/s)_min",
							       "Climate-Wind_Speed_At_10m(m/s)_max",
							       "Climate-Wind_Speed_At_10m(m/s)_median",
							       "Climate-Wind_Speed_At_10m(m/s)_stdDev",
							       "Water-Altitude(m)", "Concentration-Sex"],
		                            dependent_value=["log_Concentration-Physical_feature-SVL(mm)",
		                                             "log_Concentration-Physical_feature-Body_mass(g)",
		                                             "log_Concentration-Physical_feature-CF(100*BM(g)/SVL(mm)³)"],
		                                             log_transform=False, log_base=10)
		print("___________________________________\nAll preprocessing done!\n___________________________________")

	def read_df(self):
		df = pd.read_excel(self.str_df_path, sheet_name=self.str_sheetname, header=[0, 1, 2])  # Load the data
		print("DataFrame is successfully read!")
		return df

	'''
	Preprocessing the data, 
	including: (1) Rename the head;
			   (2) Digitalize the Gosner Stage;
			   (3) Convert data in wet weight into dry weight;
			   (4) Remove abnormal values
	'''

	def rename_columns(self):
		"""
		Rename the column heads
		"""

		list_new_head = []
		# Iterate each column
		for n in range(self.df.shape[1]):
			# Read the head of the current column
			lh = self.df_head[n]
			# Define a new head of the current column
			temp_head = ""
			i = 0
			for h in lh:
				if i == 0:
					temp_head += h.replace(" ", "_")
				elif i != 0 and "Unnamed" not in h:
					temp_head += "-" + h.replace(" ", "_")
				i += 1
			list_new_head.append(
				temp_head.replace("Physicochemical_property-", "").replace("(mg/kg)", "").replace("(mg/L)", ""))
		self.df.columns = list_new_head  # Update the head
		self.df_head = list_new_head
		print("Rename done!")

	def get_slices(self):
		"""
		Get the column ranges of different body parts
		"""

		for num, part in enumerate(self.dict_parts):
			i = 0
			for lnh in list(self.df.columns):
				if part in lnh:
					self.dict_parts[part].append(i)
				i += 1
			self.dict_parts[part] = [np.min(self.dict_parts[part]), np.max(self.dict_parts[part]) + 1]
		print("DataFrame is successfully divided!")

	def ww2dw(self, df, stage):
		"""
		Convert concentrations in wet weight to dry weight
		___________
		Parameters:
			df: A dataframe with contaminant concentrations in wet weight
			stage: Discrete stages of Anurans

		________
		Returns:
			df_dw: A dataframe with all dry weigt data
		"""
		df_dw = df
		if isinstance(stage, int):
			if stage >= 47:
				df_dw = df / (100 - self.dict_moisture["Adult"]) * 100
			elif 42 <= stage < 47:
				df_dw = df / (100 - self.dict_moisture["Metamorphosis"]) * 100
			elif 26 <= stage < 42:
				df_dw = df / (100 - self.dict_moisture["Larvae"]) * 100
			elif 20 <= stage < 26:
				df_dw = df / (100 - self.dict_moisture["Hatchlings"]) * 100
			elif stage < 20:
				df_dw = df / (100 - self.dict_moisture["Embryos"]) * 100
		return df_dw

	def execute_ww2dw(self):
		df_temp = self.df.copy()
		# Initialize the average moisture of each stage with two columns ["GS", "Moisture"]
		df_moisture = pd.DataFrame(data=None, columns=["GS", "Moisture"])
		for lnh in df_temp.columns:
			if "Gosner_Stage" in lnh:
				df_temp.replace(
					{lnh: {"Adult": 47, "null": np.nan, "Immatured": np.nan, "Tadpole": np.nan, "Juvenile": 30.5}},
					inplace=True)
				df_moisture["GS"] = df_temp[lnh]
			elif "Moisture" in lnh:
				df_moisture["Moisture"] = list(self.df[lnh])
		self.df = df_temp.copy()
		# Drop null rows
		df_moisture.dropna(axis=0, how="any", inplace=True)

		# Iterate each stage to get its corresponding default moisture
		for num, dsd in enumerate(self.dict_Stage_divide):
			if df_moisture[
				(self.dict_Stage_divide[dsd][0] <= df_moisture["GS"]) & (
						df_moisture["GS"] < self.dict_Stage_divide[dsd][1])].shape[
				0] != 0:
				# Update by replacing the default moisture values with the mean values
				self.dict_moisture[dsd] = \
					np.mean(df_moisture[(self.dict_Stage_divide[dsd][0] <= df_moisture["GS"]) & (
							df_moisture["GS"] < self.dict_Stage_divide[dsd][1])]["Moisture"])

		# Get the indexes of the lines in wet weight
		list_dwww = list(self.df[self.df["dry/wet/lipid_weight"] == "ww"].index)
		# Get concentration slice column range
		list_range = [i for k, v in self.dict_parts.items() for i in v]
		# Iterate each line in wet weight
		for ld in list_dwww:
			# Get the Gosner Stage in the current line
			num_Gosner_stage = self.df.at[ld, "Concentration-Gosner_Stage"]
			# Replace the date in wet weight with one in dry weight
			self.df.iloc[ld, min(list_range):max(list_range)] = self.ww2dw(
				df=self.df.iloc[ld, min(list_range):max(list_range)], stage=num_Gosner_stage)

		del df_temp, df_moisture, list_dwww, list_range
		print("All data are in dry weight unit!")

	def delete_abnormal_values(self):
		"""
		Remove the abnormal values
		Fill null values with the mean values
		"""
		# Get the Slice from the Original Data. Only with HM columns
		df_slice = self.df.iloc[:, self.dict_parts["Whole_Body"][0]:self.dict_parts["Bone"][-1]]
		# Load the Head of the DataFrame. Create a New List to Update the Head
		list_head = list(df_slice)
		list_new_head = []
		list_sort = {}

		# Select the Heavy Metal Columns
		for num, h in enumerate(list_head):
			str_element_name = h.split("(")[0].split("-")[-1]
			list_new_head.append(str_element_name)
			q1 = df_slice[h].quantile(0.25)
			q3 = df_slice[h].quantile(0.75)
			iqr = q3 - q1
			# Determine the Range to Distinguish the Abnormal Values
			num_lower_limit = q1 - 1.5 * iqr
			num_upper_limit = q3 + 1.5 * iqr
			# Calculate the Mean Value of the Column
			num_mean = np.mean(df_slice[h])
			# df_slice.fillna({h: num_mean}, inplace=True)
			list_sort[num_mean] = str_element_name
			# Iterate Each Value of the Current Column
			for i in df_slice[h]:
				if i > num_upper_limit or i < num_lower_limit:
					# Replace the Abnormal Values with the Mean Value of Corresponding Columns
					df_slice.replace({h: {i: num_mean}}, inplace=True)

		for lh in list_head:
			self.df[lh] = df_slice[lh]

		del df_slice, list_head, list_new_head, list_sort
		print("Abnormal data are removed!")

	def log_convert(self, log_base=10):
		"""
		The body size features and the bioaccumulated HMs are right-skewed.
		It is necessary to do log-transformation so that the distributions are normal
		___________
		Parameters:
			log_base: Logarithm base
		"""

		df_trans = pd.DataFrame(data=None, columns=self.df_head)
		list_log_head = []  # New head for log-transformed dataframe
		# Convert the parameters into new columns
		list_parameters = ["Concentration-Physical_feature-SVL(mm)",
		                   "Concentration-Physical_feature-Body_mass(g)",
		                   "Concentration-Physical_feature-CF(100*BM(g)/SVL(mm)³)"
		                   ] + self.df_head[self.dict_parts["Whole_Body"][0]:self.dict_parts["Bone"][
			-1]]  # A list of Right-skewed features that needed log-transformed
		for head in self.df_head:   # Iterate each column
			if head not in list_parameters: # If the column needs no log-transformation
				list_log_head.append(head)
				df_trans[head] = self.df[head]  # Reserve the original column for "df_trans
			else:   # If the column needs to be log-transformed
				list_log_head.append(f"log_{head}") # Append a new head
				df_trans[head] = log_transform_method[str(log_base)](self.df[head])  # Only apply log transformation

		self.df = df_trans
		self.df.columns = list_log_head
		del df_trans
		print("Data are log-transformed!")

	def pseudo_data(self):
		"""
		Convert string-like data to pseudo code
		"""

		list_head = list(self.df)
		i = 0
		list_pseudo_col = []
		for c in list_head:
			for ls in self.list_2pseudo:
				if ls in c:
					list_pseudo_col.append(c)
			i += 1
		df_pseudo = self.df.loc[:, list_pseudo_col]

		for num_col, str_col in enumerate(list_pseudo_col):
			list_col = list(set(df_pseudo[str_col]))
			dict_col = {}
			for i, items in enumerate(list_col):
				dict_col[items] = i
			df_pseudo.replace({str_col: dict_col}, inplace=True)
			self.dict_pseudo[str_col] = dict_col
		self.df[list_pseudo_col] = df_pseudo

		for key, value in self.dict_pseudo.items():
			self.dict_pseudo_backward[key] = {sub_value: sub_key for sub_key, sub_value in value.items()}

		with open(self.str_working_directory + "\\pseudo_code.txt", "w", encoding="utf-8") as pseudo:
			pseudo.write(str(self.dict_pseudo_backward).replace("nan","np.nan"))    # Replace to make sure "nan" can be identified as np.nan (that is, a null value)

		print("String data are converted into pseudo code!")

	def pseudo_data_backward(self, df):
		"""
		Convert the digitized string data backward to their corresponding strings

		___________
		Parameters:
			df: A dataframe with pseudo codes indicating the string-like data
		_______
		Return:
			 df: A dataframe with string data
		"""

		for key, value in self.dict_pseudo.items():
			dict_temp = {sub_value: sub_key for sub_key, sub_value in value.items()}
			df.replace({key: dict_temp}, inplace=True)
		print("All pseudo codes are converted backward to string type!")
		return df

	def impute_linear(self, df):
		"""
		Impute with mean values
		"""
		mean_impute = SimpleImputer(missing_values=np.nan, strategy='mean')
		df_impute_linear = pd.DataFrame(mean_impute.fit_transform(df), columns=df.columns)
		return df_impute_linear

	def impute_knn(self, df, k_values=5, phylo=None):
		"""
		Impute by KNN method
		k_values: a list of possible k values
		"""

		X_filled_knn = KNN(k=k_values, verbose=False, phylo = phylo, use_argpartition=True).fit_transform(df)
		df_impute_knn = pd.DataFrame(X_filled_knn, columns=list(df))
		return df_impute_knn

	def impute_phylo_knn(self, df, phylo_cov=None):
		return self.impute_knn(df, phylo = phylo_cov)

	def impute_soft(self, df):
		df_impute_soft = pd.DataFrame(SoftImpute().fit_transform(df), columns=df.columns)
		return df_impute_soft

	def fill_blanks(self):
		"""
		Fill empty cells
		"""

		df_species = self.df["Species"]
		list_drop = ['ID', 'Time', 'Year', 'Season', 'Location', 'Sample_code', 'Code_of_Source', 'Source',
		             'Water-Temperature(℃)', "dry/wet/lipid_weight", "Species"]
		df_2fill = self.df.iloc[:, 0:self.dict_parts["Bone"][-1]].drop(list_drop, axis=1).copy()    # Drop columns

		if self.str_fill_method.upper() == "KNN":
			for i in tqdm(range(1, KNN_iter_times + 1, 1), desc="Executing KNN"):
				X_filled_knn = KNN(k=5, verbose=False).fit_transform(df_2fill)
				df_impute = pd.DataFrame(X_filled_knn, columns=list(df_2fill))
				df_impute["k_value"] = i
				df_impute = pd.concat([df_impute, df_species], axis=1)
				df_impute.to_excel(f"{self.str_iteration_path}\\KNN_imputation_step{i:03}.xlsx", index=False)

		elif self.str_fill_method == "PhyloKNN":
			def build_phylo_vcv(phylo_tree):
				"""
				Build a phylogenetic covariance matrix from phylogenetic tree

				___________
				Parameters:
					phylo_tree: Bio.Phylo.BaseTree.Tree

				________
				Returns:
					numpy.ndarray: Covariance matrix
				"""

				list_nodes = phylo_tree.get_terminals()  # A list of tips
				list_species = [node.name for node in list_nodes]
				n_sp = phylo_tree.count_terminals()  # Number of Species
				vcv_matrix = np.zeros((n_sp, n_sp))  # Create a zero-like matrix

				dict_node2root = {}  # A dictionary to store the distances from the nodes to the root

				def calc_dist2root(node, dict_current):
					dict_node2root[node] = dict_current
					for child in node.clades:
						calc_dist2root(child, dict_current + child.branch_length)

				calc_dist2root(phylo_tree.root, 0.0)

				for i, tip1 in enumerate(list_nodes):  # Iterate each row
					for j, tip2 in enumerate(list_nodes):  # Iterate each column
						if i == j:  # If cell [i,j] is on diagonal
							vcv_matrix[i, j] = dict_node2root[tip1]  # Get distances between tips and the tree root
						else:  # If cell [i,j] is NOT on diagonal
							mrca = phylo_tree.common_ancestor(tip1, tip2)  # Find common ancestor of "tip1" and "tip2"
							vcv_matrix[i, j] = dict_node2root[mrca]  # Get distances between mrca and the tree root

				dict_vcv = {}
				num_slice = 0
				for i in range(n_sp):
					sp1 = list_nodes[i].name
					for j in range(num_slice, n_sp):
						sp2 = list_nodes[j].name
						dict_vcv[frozenset({sp1, sp2})] = vcv_matrix[i, j]  # Make a set frozen so that it can be a hashable key
					num_slice += 1
				vcv_matrix = pd.DataFrame(vcv_matrix, columns=list_species, index = list_species)
				return vcv_matrix, dict_vcv

			def compare(list_species_df, list_species_tree):
				# Check if the species is compatible between the dataframe and the phylogenetic tree
				for sp in set(list_species_df):  # Iterate each species in the dataframe
					if sp not in list_species_tree:  # If this species is not in the tree
						raise ValueError(f"{sp} not in phylogenetic tree. It has been deleted in the dataframe.")

			def build_full_phylo_cov(phylo_tree, list_species_df, list_species_tree):
				# Initialize
				compare(list_species_df,
				        list_species_tree)  # Compare species in both list to make sure every species (in "list_df_species") available from "list_tree_species"
				vcv_matrix, dict_vcv = build_phylo_vcv(phylo_tree)  # Get a dict that denotes the Patristic distance

				# Build a species matrix for dataframe as a covariance matrix
				num_rows = len(list_species_df)  # Get the number of rows
				df_cov = pd.DataFrame(np.zeros((num_rows, num_rows)))  # Build a zero-like dataframe
				array_df_species = np.array(
					list_species_df)  # Convert the list of species to an array, used to get col/row indexes
				for sp1 in tqdm(list_species_tree,desc="Building the phylogenetic matrix"):  # Iterate each column. Use "tqdm" to show the progress bar
					id_col = np.where(array_df_species == sp1)[0].tolist()  # Indexes represent "sp1"
					for sp2 in list_species_tree:  # Iterate each species to get the corresponding indexes.
						id_row = np.where(array_df_species == sp2)[0].tolist()  # Indexes represent "sp2"
						df_cov.iloc[id_row, id_col] = dict_vcv[
							frozenset({sp1, sp2})]  # Fill cells with values from "dict_vcv"

				array_cov = np.tril(np.asarray(df_cov))
				array_cov = array_cov + array_cov.T
				return vcv_matrix, array_cov

			def partial_impute(df, df_ref, iteration=0, method="weight", weight_matrix=None):
				"""
				Impute a dataframe column by column. A imputed datum is acquired from the weighted observed data.

				__________
				Parameters:
					df: A n×m dataframe with missing data with digitalize data, not the string-like data.
					df_ref: A n×1 dataframe-like column. Here we set "Species" as the reference column
					method: Partial imputation method is chosen among "mean", "median" and "weight". If "weight" is chosen, the next parameter "weight_matrix" should be given.
					weight_matrix: Weight matrix is based on the correlation of the reference column. Given that we set "Species" as our reference column, this matrix represents the MRCA (distance of between the phylogenetic tree and the common ancestor of a pair of species)

				_______
				Return:
					 df: Imputed matrix

				__________
				For example:
				Suppose we have a "df" (5×2)
					Feature1    Feature2
				0   a   NA
				1   NA  d
				2   e   f
				3   g   NA
				4   i   j

				"df_ref" (5×1)
					Species
				0   sp1
				1   sp2
				2   sp2
				3   sp3
				4   sp1

				"weight_matrix"(3×3, its shape is determined by the "df_ref" content)
						sp1 sp2 sp3
				sp1     w1  w2  w3
				sp2     w2  w1  w4
				sp3     w3  w4  w1

				Here we have null cells in df[1,0], df[0,1] and df[3,1]
				df[1,0] = w2*a + w1*e + w4*g + w2*i
				df[0,1] = w2*d + w2*f + w1*j
				df[3,1] = w4*d + w4*f + w3*j
				"""
				num_row, num_col = df.shape  # Get the shape of "df"
				list_head = list(df)  # Get the head of "df"
				head_ref = list(df_ref)[0]  # Get the head of "df_ref"
				num_pre_fill = int(num_row / self.num_iter)  # The number of row to fill in advance in every iteration
				num_pre_fill = num_pre_fill if num_pre_fill > 1 else int(num_row * 0.2)  # Re-define the number to fill
				mask_non_missing_dict = {head: np.where(~np.isnan(df[head]))[0] for head in
				                         list_head}  # Find null position in each column

				# Partial Imputation
				# Fill a slice of rows according to the phylogenetic distance
				list_row_range = np.arange(iteration*num_pre_fill, (iteration+1)*num_pre_fill,num_pre_fill).astype(int)  # Partial imputation rows
				df_slice = df.iloc[list_row_range, :]  # Extract dataframe slice
				for head in list_head:  # Iterate each column
					mask_missing = np.where(np.isnan(df_slice[head]))[0]
					if len(mask_missing) > 0:  # If there is any missing data in column "head" within "df_slice"
						df_observed = df.loc[mask_non_missing_dict[head], head]  # Observed data position in column "head"
						# list_ref_observed = [self.dict_pseudo_backward[head_ref][i] for i in df_ref[df_observed.index]]
						# list_ref_missing = set(self.dict_pseudo_backward[head_ref][i] for i in df_ref[df_slice.index])
						if method == "weight":
							for id_missing in df_slice.index:  # Iterate each null cell in the column "head" within sliced rows
								list_id_observed = df_observed.index  # Get the value index of "df_ref"
								weight=weight_matrix[np.ix_([id_missing], list_id_observed)]    # Get weight matrix of observed data
								# Fill the null position
								if np.all(np.isnan(weight)):    # If all cells are null
									weight[:] = 1    # Fill the weight matrix with 1
								elif np.any(np.isnan(weight)):  # If there is at least one cell but not fully null.
									mask = np.isnan(weight) # Get the mask of missing data
									weight[mask] = np.nanmean(weight)  # Fill the weight matrix with the mean
								sum_weight = np.sum(weight) # Get the sum of the weight
								if sum_weight != 0:  # Divide in case of a zero dominator
									weight = weight / sum_weight    # Scale to make the sum of the weight matrix one
								# Add random fluctuation to the weight matrix
								np.random.seed(iteration)
								array_random = np.random.normal(loc=np.mean(weight), scale=np.std(weight), size=weight.shape)
								mask_negative = np.where(array_random<0)
								array_random[mask_negative] = 0
								weight += array_random * 0.01
								df.loc[id_missing, head] = weight @ df_observed  # Impute with weighted observed data
						elif method == "mean":
							for id_missing in df_slice.index:  # Iterate each null cell in the column "head" within sliced rows
								df.loc[id_missing, head] = np.mean(
									df_observed)  # Impute with the mean of observed data
						elif method == "median":
							for id_missing in df_slice.index:  # Iterate each null cell in the column "head" within sliced rows
								df.loc[id_missing, head] = np.median(
									df_observed)  # Impute with the median of observed data
				return df

			# Get a phylogenetic distance matrix "phylo_cov"
			tree = Phylo.read(self.str_tree_path, "newick")
			list_tree_species = [t.name for t in tree.get_terminals()]
			list_df_species = [self.dict_pseudo_backward["Species"][sp].replace(" ", "_") for sp in list(df_species)]
			vcv, phylo_cov = build_full_phylo_cov(tree, list_df_species, list_tree_species)

			# Iterate and impute
			# Including "Partial Imputation" and "Full Imputation"
			# If we only impute directly with KNN, the imputed results are highly homogeneous or almost the same
			# Partial imputation: Add potential fluctuation to the data
			# Full Imputation: Impute all null cells basing on the partially imputed data
			for n in tqdm(range(0, self.num_iter), desc="Executing Phylo-KNN"):  # Iterate each slice
				df_2fill = partial_impute(df_2fill, df_species, n, method="weight",weight_matrix=phylo_cov)   # Partial imputation
				df_impute = self.impute_phylo_knn(df_2fill,phylo_cov)   # Fully impute basing on the previously imputed
				df_impute = pd.concat([df_impute,df_species],axis=1)
				df_impute.to_excel(f"{self.str_iteration_path}\\PhyloKNN_imputation_step{n:03}.xlsx", index=False)

		elif self.fill_method == "Soft":
			for i in tqdm(range(1, self.num_iter), desc="Executing Phylo-KNN"):
				df_impute = self.impute_soft(df_2fill)
				df_impute = pd.concat([df_impute, df_species], axis=1)
				df_impute.to_excel(f"{self.str_iteration_path}\\Soft_imputation_step{i:03}.xlsx", index=False)

	def mean_interpolation(self):
		df_whole_knn = []

		for dp in tqdm(os.listdir(self.str_iteration_path)):
			# print(dp)
			file_path = self.str_iteration_path + dp
			df_knn = pd.read_excel(file_path)
			if self.str_fill_method == "KNN":
				df_knn.drop("k_value", axis=1, inplace=True)
			df_whole_knn.append(df_knn)
		stacked_array = np.dstack([d.values for d in df_whole_knn])
		# if mean_method == "arithmetic":
		mean_array = np.mean(stacked_array, axis=2)
		# elif mean_method == "geometric":
		# 	mean_array ==

		self.df_impute_tree = pd.DataFrame(mean_array, index=df_whole_knn[0].index, columns=df_whole_knn[0].columns)
		str_savemean = self.str_working_directory+"\\02Mean_imputation\\"
		if not os.path.exists(str_savemean):
			os.mkdir(str_savemean)
		self.df_impute_tree.to_excel(str_savemean + "\\mean.xlsx")

	def linear_modeling_revise(self, df, independent_value=None, dependent_value=None, log_transform=True, log_base=10):
		if independent_value is None:
			independent_value = []
		if dependent_value is None:
			dependent_value = []

		X = df[independent_value]
		if log_transform:
			print("It is recommended not to execute logarithm transformation. The subsequent outcome is transformed")
			for dv in dependent_value:
				y = log_transform_method[str(log_base)](df[dv])
				lr = LinearRegression().fit(X, y)
				list_d = [lr.coef_[n] * df[independent_value[n]] for n in range(len(independent_value))]
				df[dv] = df[dv] - sum(list_d)
		else:
			for dv in dependent_value:
				y = df[dv]
				lr = LinearRegression().fit(X, y)
				list_d = [lr.coef_[n] * df[independent_value[n]] for n in range(len(independent_value))]
				df[dv] = df[dv] - sum(list_d)
		df = self.pseudo_data_backward(df)
		str_save_path = self.str_working_directory + "03Revised\\"
		if not os.path.exists(str_save_path):
			os.mkdir(str_save_path)
		df.to_excel(str_save_path + f"Preprocessed{self.num_date}.xlsx", index=False)


if __name__ == "__main__":
	today = datetime.now()
	formatted_date = today.strftime("%Y%m%d")
	# formatted_date = 20250716
	N_SPLITS = 4
	iter_times = 500
	log_transform_method = {
		"10": lambda x: np.log10(x),
		"2": lambda x: np.log2(x),
		"e": lambda x: np.log(x),
	}
	warnings.filterwarnings("ignore", category=FutureWarning)
	working_directory = f"E:\\学习\\研二上\\重金属\\"	# Change with your working directory

	pp = Preprocessing(str_dir=working_directory,
	                   str_df_path="C:\\Users\\admin\\Desktop\\繁殖1.xlsx",	# Change with the raw data path
	                   str_sheetname="Sheet1",
	                   fill_method="PhyloKNN",
	                   delete_abnormal_data=True,
	                   num_date=formatted_date,
	                   num_iter=iter_times)

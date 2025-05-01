import csv
import matplotlib.pyplot as plt
import statistics as stat
import numpy as np
import scipy as scipy
import math
import os
from matplotlib.ticker import (AutoMinorLocator, FormatStrFormatter)
from matplotlib.patches import FancyArrowPatch

# Define Dixon's Q Tables
# dixon_q_table = {
#   "q90" : {"3" : 0.941, "4" : 0.765, "5" : 0.642, "6" : 0.560, "7" : 0.507, "8" : 0.468, "9" : 0.437, "10" : 0.412},
#   "q95" : {"3" : 0.970, "4" : 0.829, "5" : 0.710, "6" : 0.625, "7" : 0.568, "8" : 0.526, "9" : 0.493, "10" : 0.466},
#   "q99" : {"3" : 0.994, "4" : 0.926, "5" : 0.821, "6" : 0.740, "7" : 0.680, "8" : 0.634, "9" : 0.598, "10" : 0.568}
# }

# # Function which uses Dixon's Q Test to test for and remove outliers from a list
# def remove_outliers(data_set: list, confidence_level: str) -> list:
#   data_set.sort()
#   output = data_set
#   while True:
#     found_outlier = False
#     # Test if smallest value is outlier
#     if(len(data_set) <= 2):
#       break
#     data_range = data_set[-1] - data_set[0]
#     if(data_range > 0):
#       q_min = (data_set[1] - data_set[0]) / data_range
#       if q_min > dixon_q_table[confidence_level][str(len(data_set))]:
#         # Smallest value is outlier, remove from list
#         data_set.pop(0)
#         found_outlier = True
    
#     # Test if largest value is outlier
#     if(len(data_set) <= 2):
#       break
#     data_range = data_set[-1] - data_set[0]
#     if(data_range > 0):
#       q_max = (data_set[-1] - data_set[-2]) / data_range
#       if q_max > dixon_q_table[confidence_level][str(len(data_set))]:
#         # Largest value is outlier, remove from list
#         data_set.pop(-1)
#         found_outlier = True
#     # If no outlier found, end loop and return data set
#     if found_outlier == False:
#       break
#   return data_set


# Define a class (data structure) for each age-region gene expression CT values
# Telling my coputer that I will be creating a data set in this structure. It's like a bluprint.
class age_region_data:
  def __init__(self, age, region):
    self.age = age
    self.region = region
    self.ct = []                                                  #list of CTs?
    self.average_ct = 0.0
    self.outlier = False 

# Define a key-value structure. Each gene is a key, with its associated value being many (a list of) `age_region_data`s
extractedData_N1 = {
  "-RT RPL13" : [],
  "EWS" : [],
  "FUS" : [],
  "RPL13" : [],
  "TAF15" : [],
  "TBP" : [],
  "Ube2d2a" : []
}
extractedData_N2 = {
  "-RT RPL13" : [],
  "EWS" : [],
  "FUS" : [],
  "RPL13" : [],
  "TAF15" : [],
  "TBP" : [],
  "Ube2d2a" : []
}
extractedData_N3 = {
  "-RT RPL13" : [],
  "EWS" : [],
  "FUS" : [],
  "RPL13" : [],
  "TAF15" : [],
  "TBP" : [],
  "Ube2d2a" : []
}
extractedData_N4 = {
  "-RT RPL13" : [],
  "EWS" : [],
  "FUS" : [],
  "RPL13" : [],
  "TAF15" : [],
  "TBP" : [],
  "Ube2d2a" : []
}

# Declare common variables that may be needed
genes = ["-RT RPL13", "EWS", "FUS", "RPL13", "TAF15", "TBP", "Ube2d2a"] # Must match names at start of file names
# `ages` and `regions` names MUST match (inc. capitalisation) the name from CSV file 
ages = ["3m", "12m", "24m"]
regions = ["Brain Ctrl",
           "Spinal Cord Ctrl",
           "Cortex",
           "Hippocampus",
           "Striatum",
           "Cerebellum",
           "Spinal Cord",
           "ROB"]

Nx_File_Dir = "~/Desktop/python_qPCR_analysis/"
Nx_FileNames = ["processed_N1_data.csv",
                "processed_N2_data.csv",
                "processed_N3_data.csv",
                "processed_N4_data.csv"]

#--------------------------- Read In All Data ---------------------------------------

# Define quick function to add record to correct dictionary
def addToDic(N, record, gene):
  if N == 1:
    extractedData_N1[gene].append(record)
  elif N == 2:
    extractedData_N2[gene].append(record)
  elif N == 3:
    extractedData_N3[gene].append(record)
  else: # N == 4
    extractedData_N4[gene].append(record)
  # EoF

N = 0 # variable to track which N set we are in
for file in Nx_FileNames:
  N = N + 1
  # Open each file in Nx_FileNames list
  with open(Nx_File_Dir + file, mode='r') as inputFile:
    csvReader = csv.reader(inputFile)
    # Go through all rows and save data into correct Dictionary
    for row in csvReader:
      # Format of a row : [Gene, Age, Region, ct_0, ct_1, ct_2, ct_avg, Outlier]
      # Discard Header row
      if row[0] == "Gene":
        continue # Continue = done with this iteration of for loop, start next one
      ar_record = age_region_data(row[1], row[2]) # (age, region)
      ar_record.ct = [float(row[3]), float(row[4]), float(row[5])]
      # CHANGE MADE - re-calculate average ct so that all ct values are always used. Outliers are still marked as such (i.e. sets where all three values are ±0.5 from each other)
      # ar_record.average_ct = float(row[6])
      ar_record.average_ct = sum(ar_record.ct) / len(ar_record.ct)
      if row[7] == "True":
        ar_record.outlier = True
      else :
        ar_record.outlier = False
      addToDic(N, ar_record, row[0])

#------------------------------------------------------------------------------------


#----------------------------- Analysis 1 - FET Regional Expression ---------------------



 # N_number_data {
 #      3_month_data {
 #          RPL13 {
 #              Brain_Total {
 #                  Average_CT
 #                  Delta_CT (THIS-GENE.Brain_Total.AverageCT - AverageCT)
 #                  Delta_Delta_CT (DeltaCT - HOUSEKEEPING-GENE.THIS-REGION.DeltaCT)
 #                  2_^_Delta_Delta_CT (2 ^ Delta_Delta_CT)                 <------- PLOT ALL THESE ONES
 #              }
 #              Cortex
 #              Hippocampus
 #              Striatum
 #              Cerebellum
 #              ROB
 #          }
 #          FUS
 #          EWS
 #          TAF15
 #      }
 #      12_month_data {
 # 
 #      }
 #      24_month_data {
 # 
 #      }
 # }
class region_values_1:
    def __init__(self, region):
      self.region_name = region
      self.average_ct = 0.0
      self.d_ct = 0.0
      self.dd_ct = {
        "RPL13" : 0.0,
        "TBP" : 0.0,
        "Ube2d2a" : 0.0
      }
      self.two_pow_ddct = {
        "RPL13" : 0.0,
        "TBP" : 0.0,
        "Ube2d2a" : 0.0
      }
      self.outlier = False

class gene_regions_1:
    def __init__(self, gene):
       self.gene_name = gene
       self.brain_total = region_values_1("Brain Total")
       self.cortex = region_values_1("Cortex")
       self.hippocampus = region_values_1("Hippocampus")
       self.striatum = region_values_1("Striatum")
       self.cerebellum = region_values_1("Cerebellum")
       self.rob = region_values_1("ROB")

class N_Num_Data_1:
    def __init__(self):
        self.data = {
          "3m" : {
            "RPL13" : gene_regions_1("RPL13"),
            "FUS" : gene_regions_1("FUS"),
            "EWS" : gene_regions_1("EWS"),
            "TAF15" : gene_regions_1("TAF15"),
            "TBP" : gene_regions_1("TBP"),
            "Ube2d2a" : gene_regions_1("Ube2d2a")
          },
          "12m" : {
            "RPL13" : gene_regions_1("RPL13"),
            "FUS" : gene_regions_1("FUS"),
            "EWS" : gene_regions_1("EWS"),
            "TAF15" : gene_regions_1("TAF15"),
            "TBP" : gene_regions_1("TBP"),
            "Ube2d2a" : gene_regions_1("Ube2d2a")
          },
          "24m" : {
            "RPL13" : gene_regions_1("RPL13"),
            "FUS" : gene_regions_1("FUS"),
            "EWS" : gene_regions_1("EWS"),
            "TAF15" : gene_regions_1("TAF15"),
            "TBP" : gene_regions_1("TBP"),
            "Ube2d2a" : gene_regions_1("Ube2d2a")
          }
        }


Ns_1 = [N_Num_Data_1(), N_Num_Data_1(), N_Num_Data_1(), N_Num_Data_1()]
#    [ N1 , N2 , N3 , N4 ]

outputDir_1 = "ANALYSIS-1_Regional_FET_expression/"
# make sure outputDir_1 exists
if not os.path.exists(outputDir_1):
  os.makedirs(outputDir_1)

target_genes = ["FUS", "EWS", "TAF15", "RPL13", "TBP", "Ube2d2a"]
target_regions = ["Brain Total", "Cortex", "Hippocampus", "Striatum", "Cerebellum"]#, "ROB"]
housekeeping_genes = ["RPL13", "TBP", "Ube2d2a"]
FET_genes = ["FUS", "EWS", "TAF15"]

def init_avg_ct_1(Nx, gene, age, region, record):
  # Must -1 from Nx when we use it so that it goes to the correct N
  # i.e. Ns_1[0] corresponds to N1 --> therefore a Nx=1 corresponding to N1, we must -1
  if region == "Brain Ctrl":
    Ns_1[Nx - 1].data[age][gene].brain_total.average_ct = record.average_ct
    Ns_1[Nx - 1].data[age][gene].brain_total.outlier = record.outlier
  if region == "Cortex":
    Ns_1[Nx - 1].data[age][gene].cortex.average_ct = record.average_ct
    Ns_1[Nx - 1].data[age][gene].cortex.outlier = record.outlier
  if region == "Hippocampus":
    Ns_1[Nx - 1].data[age][gene].hippocampus.average_ct = record.average_ct
    Ns_1[Nx - 1].data[age][gene].hippocampus.outlier = record.outlier
  if region == "Striatum":
    Ns_1[Nx - 1].data[age][gene].striatum.average_ct = record.average_ct
    Ns_1[Nx - 1].data[age][gene].striatum.outlier = record.outlier
  if region == "Cerebellum":
    Ns_1[Nx - 1].data[age][gene].cerebellum.average_ct = record.average_ct
    Ns_1[Nx - 1].data[age][gene].cerebellum.outlier = record.outlier
  if region == "ROB":
    Ns_1[Nx - 1].data[age][gene].rob.average_ct = record.average_ct
    Ns_1[Nx - 1].data[age][gene].rob.outlier = record.outlier


# Initialise all average_ct values
for g in target_genes:
   for r1 in extractedData_N1[g]:
      init_avg_ct_1(1, g, r1.age, r1.region, r1)
   for r2 in extractedData_N2[g]:
      init_avg_ct_1(2, g, r2.age, r2.region, r2)
   for r3 in extractedData_N3[g]:
      init_avg_ct_1(3, g, r3.age, r3.region, r3)
   for r4 in extractedData_N4[g]:
      init_avg_ct_1(4, g, r4.age, r4.region, r4) 



def calc_delta_ct_1(Nx, gene):
  # Must -1 from Nx when we use it so that it goes to the correct N
  # i.e. Ns_1[0] corresponds to N1 --> therefore a Nx=1 corresponding to N1, we must -1
  for a in ages:
    brain_total_avg_ct = Ns_1[Nx - 1].data[a][gene].brain_total.average_ct # Get brain_total average_ct used to normalise delta_CT
    Ns_1[Nx - 1].data[a][gene].brain_total.d_ct = brain_total_avg_ct - Ns_1[Nx - 1].data[a][gene].brain_total.average_ct
    Ns_1[Nx - 1].data[a][gene].cortex.d_ct = brain_total_avg_ct - Ns_1[Nx - 1].data[a][gene].cortex.average_ct
    Ns_1[Nx - 1].data[a][gene].hippocampus.d_ct = brain_total_avg_ct - Ns_1[Nx - 1].data[a][gene].hippocampus.average_ct
    Ns_1[Nx - 1].data[a][gene].striatum.d_ct = brain_total_avg_ct - Ns_1[Nx - 1].data[a][gene].striatum.average_ct
    Ns_1[Nx - 1].data[a][gene].cerebellum.d_ct = brain_total_avg_ct - Ns_1[Nx - 1].data[a][gene].cerebellum.average_ct
    Ns_1[Nx - 1].data[a][gene].rob.d_ct = brain_total_avg_ct - Ns_1[Nx - 1].data[a][gene].rob.average_ct

# Calculate all Delta_CT values --> Delta_CT = (AverageCT - THIS-GENE.Brain_Total.AverageCT)
for g in target_genes:
  for n in [1,2,3,4]:
    calc_delta_ct_1(n, g)



def calc_dd_ct_for_region_1(Nx, region, gene, age, housekeepGene):
  if region == "Brain Total":
    Ns_1[Nx - 1].data[age][gene].brain_total.dd_ct[housekeepGene] = Ns_1[Nx - 1].data[age][gene].brain_total.d_ct - Ns_1[Nx - 1].data[age][housekeepGene].brain_total.d_ct
    Ns_1[Nx - 1].data[age][gene].brain_total.two_pow_ddct[housekeepGene] = pow(2.0, Ns_1[Nx - 1].data[age][gene].brain_total.dd_ct[housekeepGene])
  if region == "Cortex":
    Ns_1[Nx - 1].data[age][gene].cortex.dd_ct[housekeepGene] = Ns_1[Nx - 1].data[age][gene].cortex.d_ct - Ns_1[Nx - 1].data[age][housekeepGene].cortex.d_ct
    Ns_1[Nx - 1].data[age][gene].cortex.two_pow_ddct[housekeepGene] = pow(2.0, Ns_1[Nx - 1].data[age][gene].cortex.dd_ct[housekeepGene])
  if region == "Hippocampus":
    Ns_1[Nx - 1].data[age][gene].hippocampus.dd_ct[housekeepGene] = Ns_1[Nx - 1].data[age][gene].hippocampus.d_ct - Ns_1[Nx - 1].data[age][housekeepGene].hippocampus.d_ct
    Ns_1[Nx - 1].data[age][gene].hippocampus.two_pow_ddct[housekeepGene] = pow(2.0, Ns_1[Nx - 1].data[age][gene].hippocampus.dd_ct[housekeepGene])
  if region == "Striatum":
    Ns_1[Nx - 1].data[age][gene].striatum.dd_ct[housekeepGene] = Ns_1[Nx - 1].data[age][gene].striatum.d_ct - Ns_1[Nx - 1].data[age][housekeepGene].striatum.d_ct
    Ns_1[Nx - 1].data[age][gene].striatum.two_pow_ddct[housekeepGene] = pow(2.0, Ns_1[Nx - 1].data[age][gene].striatum.dd_ct[housekeepGene])
  if region == "Cerebellum":
    Ns_1[Nx - 1].data[age][gene].cerebellum.dd_ct[housekeepGene] = Ns_1[Nx - 1].data[age][gene].cerebellum.d_ct - Ns_1[Nx - 1].data[age][housekeepGene].cerebellum.d_ct
    Ns_1[Nx - 1].data[age][gene].cerebellum.two_pow_ddct[housekeepGene] = pow(2.0, Ns_1[Nx - 1].data[age][gene].cerebellum.dd_ct[housekeepGene])
  if region == "ROB":
    Ns_1[Nx - 1].data[age][gene].rob.dd_ct[housekeepGene] = Ns_1[Nx - 1].data[age][gene].rob.d_ct - Ns_1[Nx - 1].data[age][housekeepGene].rob.d_ct
    Ns_1[Nx - 1].data[age][gene].rob.two_pow_ddct[housekeepGene] = pow(2.0, Ns_1[Nx - 1].data[age][gene].rob.dd_ct[housekeepGene])


# Calculate all Delta_Delta_CT values & 2^Delta_Delta_CT vlaues
# --> Delta_Delta_CT = (DeltaCT - RPL13.THIS-REGION.DeltaCT)
for hkg in housekeeping_genes:
  for g in target_genes:
    for n in [1,2,3,4]:
      for r in target_regions:
        for a in ages:
          calc_dd_ct_for_region_1(n, r, g, a, hkg)


# Export data to single CSV
def make_row_1(Nx, gene, age, region):
  row = []
  row.append(Nx)
  row.append(gene)
  row.append(age)
  row.append(region)
  if region == "Brain Total":
    row.append(Ns_1[Nx - 1].data[age][gene].brain_total.average_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].brain_total.d_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].brain_total.dd_ct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].brain_total.two_pow_ddct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].brain_total.dd_ct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].brain_total.two_pow_ddct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].brain_total.dd_ct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].brain_total.two_pow_ddct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].brain_total.outlier)
  if region == "Cortex":
    row.append(Ns_1[Nx - 1].data[age][gene].cortex.average_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].cortex.d_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].cortex.dd_ct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].cortex.two_pow_ddct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].cortex.dd_ct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].cortex.two_pow_ddct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].cortex.dd_ct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].cortex.two_pow_ddct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].cortex.outlier)
  if region == "Hippocampus":
    row.append(Ns_1[Nx - 1].data[age][gene].hippocampus.average_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].hippocampus.d_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].hippocampus.dd_ct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].hippocampus.two_pow_ddct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].hippocampus.dd_ct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].hippocampus.two_pow_ddct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].hippocampus.dd_ct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].hippocampus.two_pow_ddct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].hippocampus.outlier)
  if region == "Striatum":
    row.append(Ns_1[Nx - 1].data[age][gene].striatum.average_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].striatum.d_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].striatum.dd_ct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].striatum.two_pow_ddct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].striatum.dd_ct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].striatum.two_pow_ddct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].striatum.dd_ct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].striatum.two_pow_ddct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].striatum.outlier)
  if region == "Cerebellum":
    row.append(Ns_1[Nx - 1].data[age][gene].cerebellum.average_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].cerebellum.d_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].cerebellum.dd_ct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].cerebellum.two_pow_ddct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].cerebellum.dd_ct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].cerebellum.two_pow_ddct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].cerebellum.dd_ct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].cerebellum.two_pow_ddct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].cerebellum.outlier)
  if region == "ROB":
    row.append(Ns_1[Nx - 1].data[age][gene].rob.average_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].rob.d_ct)
    row.append(Ns_1[Nx - 1].data[age][gene].rob.dd_ct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].rob.two_pow_ddct["RPL13"])
    row.append(Ns_1[Nx - 1].data[age][gene].rob.dd_ct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].rob.two_pow_ddct["TBP"])
    row.append(Ns_1[Nx - 1].data[age][gene].rob.dd_ct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].rob.two_pow_ddct["Ube2d2a"])
    row.append(Ns_1[Nx - 1].data[age][gene].rob.outlier)
  return row

headers = ["N Replicate", "Gene", "Age", "Region", "Average CT", "Delta CT", "Delta Delta CT (Norm. against RPL13)", "2^Delta_Delta_CT (Norm. against RPL13)", "Delta Delta CT (Norm. against TBP)", "2^Delta_Delta_CT (Norm. against TBP)", "Delta Delta CT (Norm. against Ube2d2a)", "2^Delta_Delta_CT (Norm. against Ube2d2a)", "Outlier"]
with open(outputDir_1 + "ANALYSIS-1_Regional_FET_expression_in_brain_raw_data.csv", mode='w') as outfile:
  writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)  # .writer function is specific for csv file writing- allows to write rows to csv files
  writer.writerow(headers)                         # this places the defined headers lists into the a new rown in the csv file - each element in headers list corrisponds to a single cell of a single new row
  for n in [1, 2, 3, 4]:
    for g in target_genes:
      for r in target_regions:
        for a in ages:
          writer.writerow(make_row_1(n, g, a, r))


# Create plots - FUS=Blue, EWS=Greem, TAF15=Orange
#              - One plot for each gene-age-houseKeeping combination
#              - Plot scatter points for 2^ddCT of each N
#              - Plot mean of 2^ddCT across Ns_1
#              - Plot standard deviation using N 2^ddCT values and calculated means
#                         - sd = sqrt( ( sum( abs(x - mu)^2 ) ) / N )     where x is 2^ddCT from each N and mu is the mean
#                         - Plotted as error bars as ±Standard-Deviation
#              - Calculate p value for each non-Brain-Total region and plot:
#                         - If p < 0.05 then put * above region's bar
#                         - If p < 0.01 then put ** above region's bar
#                         - If p < 0.001 then put *** above region's bar
#                         - Also print actual p-value above region's bar to 4dp
#                         - p-value = two-tailed one-sample T-test

for hkg in housekeeping_genes:
  for g in FET_genes:
    for a in ages:
      # Get all data points we will need for plot
      data_points_ddCT = {
        "Brain Total" : [],
        "Cortex" : [],
        "Hippocampus" : [],
        "Striatum" : [],
        "Cerebellum" : [],
        "ROB" : []
      }
      data_points_2pow_ddCT = {
        "Brain Total" : [],
        "Cortex" : [],
        "Hippocampus" : [],
        "Striatum" : [],
        "Cerebellum" : [],
        "ROB" : []
      }
      twoPow_ddct_mean = {
        "Brain Total" : 0.0,
        "Cortex" : 0.0,
        "Hippocampus" : 0.0,
        "Striatum" : 0.0,
        "Cerebellum" : 0.0,
        "ROB" : 0.0
      }
      twoPow_ddct_sd = {
        "Brain Total" : 0.0,
        "Cortex" : 0.0,
        "Hippocampus" : 0.0,
        "Striatum" : 0.0,
        "Cerebellum" : 0.0,
        "ROB" : 0.0
      }
      significance_data = {
        "Cortex" : {
          "p_value" : 0.0,
          "p_bonferroni" : 0.0,
          "sig" : ""
        },
        "Hippocampus" : {
          "p_value" : 0.0,
          "p_bonferroni" : 0.0,
          "sig" : ""
        },
        "Striatum" : {
          "p_value" : 0.0,
          "p_bonferroni" : 0.0,
          "sig" : ""
        },
        "Cerebellum" : {
          "p_value" : 0.0,
          "p_bonferroni" : 0.0,
          "sig" : ""
        },
        "ROB" : {
          "p_value" : 0.0,
          "p_bonferroni" : 0.0,
          "sig" : ""
        }
      }
      for n in [1,2,3,4]:
        # Brain total
        data_points_ddCT["Brain Total"].append(Ns_1[n-1].data[a][g].brain_total.dd_ct[hkg])
        data_points_2pow_ddCT["Brain Total"].append(Ns_1[n-1].data[a][g].brain_total.two_pow_ddct[hkg])
        # Cortex
        data_points_ddCT["Cortex"].append(Ns_1[n-1].data[a][g].cortex.dd_ct[hkg])
        data_points_2pow_ddCT["Cortex"].append(Ns_1[n-1].data[a][g].cortex.two_pow_ddct[hkg])
        # Hippocampus
        data_points_ddCT["Hippocampus"].append(Ns_1[n-1].data[a][g].hippocampus.dd_ct[hkg])
        data_points_2pow_ddCT["Hippocampus"].append(Ns_1[n-1].data[a][g].hippocampus.two_pow_ddct[hkg])
        # Striatum
        data_points_ddCT["Striatum"].append(Ns_1[n-1].data[a][g].striatum.dd_ct[hkg])
        data_points_2pow_ddCT["Striatum"].append(Ns_1[n-1].data[a][g].striatum.two_pow_ddct[hkg])
        # Cerebellum
        data_points_ddCT["Cerebellum"].append(Ns_1[n-1].data[a][g].cerebellum.dd_ct[hkg])
        data_points_2pow_ddCT["Cerebellum"].append(Ns_1[n-1].data[a][g].cerebellum.two_pow_ddct[hkg])
        # ROB
        data_points_ddCT["ROB"].append(Ns_1[n-1].data[a][g].rob.dd_ct[hkg])
        data_points_2pow_ddCT["ROB"].append(Ns_1[n-1].data[a][g].rob.two_pow_ddct[hkg])
      # Calculate means + standard deviations
      for r in target_regions:
        # Remove outliers : possible confidence levels = {"q90", "q95", "q99"}
        # data_points_2pow_ddCT[r] = remove_outliers(data_points_2pow_ddCT[r], "q95")
        # data_points_ddCT[r] = remove_outliers(data_points_ddCT[r], "q95")
        twoPow_ddct_mean[r] = sum(data_points_2pow_ddCT[r]) / len(data_points_2pow_ddCT[r])
        twoPow_ddct_sd[r] = stat.stdev(data_points_2pow_ddCT[r])
        # Calculate significance (p-value) via One Sample T Test (two tailed)      
        if r == "Brain Total": # Ignore Brain Total
          continue
        SD_ddct = stat.stdev(data_points_ddCT[r])
        num_replicates = len(data_points_ddCT[r])
        std_err_mean = SD_ddct / math.sqrt(num_replicates)
        x_bar = sum(data_points_ddCT["Brain Total"]) / len(data_points_ddCT["Brain Total"])
        alt_hypothesis = sum(data_points_ddCT[r]) / len(data_points_ddCT[r]) # M
        degs_of_freedom = len(target_regions) - 1
        num_tails = 1
        t_val = (x_bar - alt_hypothesis) / std_err_mean
        p_val = (1 - scipy.stats.t.cdf(x=abs(t_val), df=degs_of_freedom)) * num_tails
        p_bonf = p_val * len(target_regions)
        significance_data[r]["p_value"] = p_val
        significance_data[r]["p_bonferroni"] = p_bonf
        sig_str = ""
        if p_val < 0.001:
          sig_str = "***"
        elif p_val < 0.01:
          sig_str = "**"
        elif p_val < 0.05:
          sig_str = "*"
        significance_data[r]["sig"] = sig_str



      w = 0.6 # bar width
      x = list(range(1,len(target_regions)+1)) # create list so each region has a x coord

      colour = "#99b4cb" # Colour (blue) for FUS
      ecolour = "#2b4357"
      if g == "TAF15" :
        colour = "#c7d98b" # orange
        ecolour = "#414d1e"
      if g == "EWS":
        colour = "#8ed6b7" # green
        ecolour = "#385749"

      fig, ax = plt.subplots()
      allMeans = [twoPow_ddct_mean["Brain Total"], twoPow_ddct_mean["Cortex"], twoPow_ddct_mean["Hippocampus"], twoPow_ddct_mean["Striatum"], twoPow_ddct_mean["Cerebellum"]] # twoPow_ddct_mean["ROB"]]
      allSD = [twoPow_ddct_sd["Brain Total"], twoPow_ddct_sd["Cortex"], twoPow_ddct_sd["Hippocampus"], twoPow_ddct_sd["Striatum"], twoPow_ddct_sd["Cerebellum"]] # twoPow_ddct_sd["ROB"]]
      ax.bar(x, 
             height=allMeans,
             yerr = allSD,
             capsize = 6, # error bar cap width in points #NEW SIZE
             error_kw = {"elinewidth":0.8, "capthick":0.8}, ###NEW
             width = w,
             tick_label = target_regions,
             color = colour,
             edgecolor = ecolour
             )
      
      for i in range(len(x)):
        # distribute scatter randomly across whole width of bar
        bar_x_coords = np.linspace(x[i] - (w/2.5), x[i] + w/2.5, len(data_points_2pow_ddCT[target_regions[i]]))
        ax.scatter(bar_x_coords, data_points_2pow_ddCT[target_regions[i]], color="black", s=8)
      
      # Do some editing of aesthetics
      ax.spines[['right', 'top']].set_visible(False)
      ax.yaxis.set_minor_locator(AutoMinorLocator())
      ax.tick_params(which='major', length=7)
      ax.tick_params(which='minor', length=4, color='black')
      ax.tick_params(axis='both',labelsize=8)##
      plt.xticks(rotation=45,ha='right')##

      ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

      # Add comparison caps
      for i in range(len(x)):
        if i == 0:
          continue
        region = target_regions[i]
        # if any stars, print error bar between brain total and current region
        if significance_data[region]["sig"] != "":
          xStart = x[0]
          xEnd = x[i]
          y = ax.get_ylim()[-1]
          xs = [xStart, xStart, xEnd, xEnd]
          ys = [y, y*1.1, y*1.1, y]
          plt.plot(xs, ys, lw=1, color='black')
          plt.text(xStart + ((xEnd - xStart)/2), y*1.1, significance_data[region]["sig"], ha='center', va='bottom', fontsize=12) #NEW

      # Add p-value text
      extraArt = [] ######NEW
      for i in range(len(x)):
        if i == 0:
          continue
        height = ax.get_ylim()
        region = target_regions[i]
        #               plt.text(x[i], height[-1]*1.05, "p = " + str(round(significance_data[region]["p_value"], 4)), ha='center', va='bottom', fontsize=12)
        extraArt.append(plt.text(x[i], height[-1]*1.05, "p = " + str(round(significance_data[region]["p_value"], 3)), ha='center', va='bottom', fontsize=6)) ####NEW 

      #               plt.ylabel(g + " mRNA Levels Relative to\nBrain Total Normalised to " + hkg, fontsize=14)
      extraArt.append(plt.ylabel(g + " mRNA Levels Relative to\nBrain Total Normalised to " + hkg, fontsize=10)) ## NEW
      xlabelText = "3 Months"
      if a == "12m":
        xlabelText = "12 Months"
      elif a == "24m":
        xlabelText = "24 Months"
        #             plt.xlabel(xlabelText, fontsize=14
      extraArt.append(plt.xlabel(xlabelText, fontsize=8))   #NEW       # x label added for powerpoint figures remove or change to 3 months, 12 months or 24 months
      plt.axhline(y=0, linestyle='-', color='black', lw=0.5)
      plt.tight_layout()
      fig.set_size_inches(2.8,4)
  #   plt.savefig(outputDir_1 + hkg + "_" + g + "_" + a + "_graph.png"
      plt.savefig(outputDir_1 + hkg + "_" + g + "_" + a + "_graph.png", dpi=500, bbox_inches='tight', bbox_extra_artists=tuple(extraArt)) #NEW
      plt.close()



#------------------------------------------------------------------------------------





#----------------------------- Analysis 2 - FET Expression Across Ages  ---------------------



 # N_number_data {
 #      3_month_data {
 #          RPL13 {
 #              Brain_Total {
 #                  Average_CT
 #                  Delta_CT (AverageCT - HOUSEKEEPING-GENE.THIS_REGION.AverageCT)
 #                  Delta_Delta_CT (DeltaCT - 3-MONTH.HOUSEKEEPING-GENE.THIS_REGION.DeltaCT)
 #                  2_^_Delta_Delta_CT (2 ^ Delta_Delta_CT)                 <------- PLOT ALL THESE ONES
 #              }
 #              Cortex
 #              Hippocampus
 #              Striatum
 #              Cerebellum
 #              ROB
 #          }
 #          FUS
 #          EWS
 #          TAF15
 #      }
 #      12_month_data {
 # 
 #      }
 #      24_month_data {
 # 
 #      }
 # }
class region_values_2:
    def __init__(self, region):
      self.region_name = region
      self.average_ct = 0.0
      self.d_ct = {
        "RPL13" : 0.0,
        "TBP" : 0.0,
        "Ube2d2a" : 0.0
      }
      self.dd_ct = {
        "RPL13" : 0.0,
        "TBP" : 0.0,
        "Ube2d2a" : 0.0
      }
      self.two_pow_ddct = {
        "RPL13" : 0.0,
        "TBP" : 0.0,
        "Ube2d2a" : 0.0
      }
      self.outlier = False

class gene_regions_2:
    def __init__(self, gene):
       self.gene_name = gene
       self.brain_total = region_values_2("Brain Total")
       self.cortex = region_values_2("Cortex")
       self.hippocampus = region_values_2("Hippocampus")
       self.striatum = region_values_2("Striatum")
       self.cerebellum = region_values_2("Cerebellum")
       self.rob = region_values_2("ROB")

class N_Num_Data_2:
    def __init__(self):
        self.data = {
          "3m" : {
            "RPL13" : gene_regions_2("RPL13"),
            "FUS" : gene_regions_2("FUS"),
            "EWS" : gene_regions_2("EWS"),
            "TAF15" : gene_regions_2("TAF15"),
            "TBP" : gene_regions_2("TBP"),
            "Ube2d2a" : gene_regions_2("Ube2d2a")
          },
          "12m" : {
            "RPL13" : gene_regions_2("RPL13"),
            "FUS" : gene_regions_2("FUS"),
            "EWS" : gene_regions_2("EWS"),
            "TAF15" : gene_regions_2("TAF15"),
            "TBP" : gene_regions_2("TBP"),
            "Ube2d2a" : gene_regions_2("Ube2d2a")
          },
          "24m" : {
            "RPL13" : gene_regions_2("RPL13"),
            "FUS" : gene_regions_2("FUS"),
            "EWS" : gene_regions_2("EWS"),
            "TAF15" : gene_regions_2("TAF15"),
            "TBP" : gene_regions_2("TBP"),
            "Ube2d2a" : gene_regions_2("Ube2d2a")
          }
        }


Ns_2 = [N_Num_Data_2(), N_Num_Data_2(), N_Num_Data_2(), N_Num_Data_2()]
#    [ N1 , N2 , N3 , N4 ]

outputDir_2 = "ANALYSIS-2_FET_expression_Across_Ages/"
# make sure outputDir_2 exists
if not os.path.exists(outputDir_2):
  os.makedirs(outputDir_2)

target_genes = ["FUS", "EWS", "TAF15", "RPL13", "TBP", "Ube2d2a"]
target_regions = ["Brain Total", "Cortex", "Hippocampus", "Striatum", "Cerebellum"]#, "ROB"]
housekeeping_genes = ["RPL13", "TBP", "Ube2d2a"]
FET_genes = ["FUS", "EWS", "TAF15"]

def init_avg_ct_2(Nx, gene, age, region, record):
  # Must -1 from Nx when we use it so that it goes to the correct N
  # i.e. Ns_2[0] corresponds to N1 --> therefore a Nx=1 corresponding to N1, we must -1
  if region == "Brain Ctrl":
    Ns_2[Nx - 1].data[age][gene].brain_total.average_ct = record.average_ct
    Ns_2[Nx - 1].data[age][gene].brain_total.outlier = record.outlier
  if region == "Cortex":
    Ns_2[Nx - 1].data[age][gene].cortex.average_ct = record.average_ct
    Ns_2[Nx - 1].data[age][gene].cortex.outlier = record.outlier
  if region == "Hippocampus":
    Ns_2[Nx - 1].data[age][gene].hippocampus.average_ct = record.average_ct
    Ns_2[Nx - 1].data[age][gene].hippocampus.outlier = record.outlier
  if region == "Striatum":
    Ns_2[Nx - 1].data[age][gene].striatum.average_ct = record.average_ct
    Ns_2[Nx - 1].data[age][gene].striatum.outlier = record.outlier
  if region == "Cerebellum":
    Ns_2[Nx - 1].data[age][gene].cerebellum.average_ct = record.average_ct
    Ns_2[Nx - 1].data[age][gene].cerebellum.outlier = record.outlier
  if region == "ROB":
    Ns_2[Nx - 1].data[age][gene].rob.average_ct = record.average_ct
    Ns_2[Nx - 1].data[age][gene].rob.outlier = record.outlier


# Initialise all average_ct values
for g in target_genes:
   for r1 in extractedData_N1[g]:
      init_avg_ct_2(1, g, r1.age, r1.region, r1)
   for r2 in extractedData_N2[g]:
      init_avg_ct_2(2, g, r2.age, r2.region, r2)
   for r3 in extractedData_N3[g]:
      init_avg_ct_2(3, g, r3.age, r3.region, r3)
   for r4 in extractedData_N4[g]:
      init_avg_ct_2(4, g, r4.age, r4.region, r4) 



def calc_delta_ct_2(Nx, gene, housekeepGene):
  # Must -1 from Nx when we use it so that it goes to the correct N
  # i.e. Ns_2[0] corresponds to N1 --> therefore a Nx=1 corresponding to N1, we must -1
  for a in ages:
    Ns_2[Nx - 1].data[a][gene].brain_total.d_ct[housekeepGene] = Ns_2[Nx - 1].data[a][gene].brain_total.average_ct - Ns_2[Nx - 1].data[a][housekeepGene].brain_total.average_ct
    Ns_2[Nx - 1].data[a][gene].cortex.d_ct[housekeepGene] = Ns_2[Nx - 1].data[a][gene].cortex.average_ct - Ns_2[Nx - 1].data[a][housekeepGene].cortex.average_ct
    Ns_2[Nx - 1].data[a][gene].hippocampus.d_ct[housekeepGene] = Ns_2[Nx - 1].data[a][gene].hippocampus.average_ct - Ns_2[Nx - 1].data[a][housekeepGene].hippocampus.average_ct
    Ns_2[Nx - 1].data[a][gene].striatum.d_ct[housekeepGene] = Ns_2[Nx - 1].data[a][gene].striatum.average_ct - Ns_2[Nx - 1].data[a][housekeepGene].striatum.average_ct
    Ns_2[Nx - 1].data[a][gene].cerebellum.d_ct[housekeepGene] = Ns_2[Nx - 1].data[a][gene].cerebellum.average_ct - Ns_2[Nx - 1].data[a][housekeepGene].cerebellum.average_ct
    Ns_2[Nx - 1].data[a][gene].rob.d_ct[housekeepGene] = Ns_2[Nx - 1].data[a][gene].rob.average_ct - Ns_2[Nx - 1].data[a][housekeepGene].rob.average_ct

# Calculate all Delta_CT values --> Delta_CT = (AverageCT - HOUSEKEEPING-GENE.THIS_REGION.AverageCT)
for g in target_genes:
  for hkg in housekeeping_genes:
    for n in [1,2,3,4]:
      calc_delta_ct_2(n, g, hkg)



def calc_dd_ct_for_region_2(Nx, region, gene, age, housekeepGene):
  if region == "Brain Total":
    Ns_2[Nx - 1].data[age][gene].brain_total.dd_ct[housekeepGene] = Ns_2[Nx - 1].data[age][gene].brain_total.d_ct[housekeepGene] - Ns_2[Nx - 1].data["3m"][gene].brain_total.d_ct[housekeepGene]
    Ns_2[Nx - 1].data[age][gene].brain_total.two_pow_ddct[housekeepGene] = pow(2.0, Ns_2[Nx - 1].data[age][gene].brain_total.dd_ct[housekeepGene])
  if region == "Cortex":
    Ns_2[Nx - 1].data[age][gene].cortex.dd_ct[housekeepGene] = Ns_2[Nx - 1].data[age][gene].cortex.d_ct[housekeepGene] - Ns_2[Nx - 1].data["3m"][gene].cortex.d_ct[housekeepGene]
    Ns_2[Nx - 1].data[age][gene].cortex.two_pow_ddct[housekeepGene] = pow(2.0, Ns_2[Nx - 1].data[age][gene].cortex.dd_ct[housekeepGene])
  if region == "Hippocampus":
    Ns_2[Nx - 1].data[age][gene].hippocampus.dd_ct[housekeepGene] = Ns_2[Nx - 1].data[age][gene].hippocampus.d_ct[housekeepGene] - Ns_2[Nx - 1].data["3m"][gene].hippocampus.d_ct[housekeepGene]
    Ns_2[Nx - 1].data[age][gene].hippocampus.two_pow_ddct[housekeepGene] = pow(2.0, Ns_2[Nx - 1].data[age][gene].hippocampus.dd_ct[housekeepGene])
  if region == "Striatum":
    Ns_2[Nx - 1].data[age][gene].striatum.dd_ct[housekeepGene] = Ns_2[Nx - 1].data[age][gene].striatum.d_ct[housekeepGene] - Ns_2[Nx - 1].data["3m"][gene].striatum.d_ct[housekeepGene]
    Ns_2[Nx - 1].data[age][gene].striatum.two_pow_ddct[housekeepGene] = pow(2.0, Ns_2[Nx - 1].data[age][gene].striatum.dd_ct[housekeepGene])
  if region == "Cerebellum":
    Ns_2[Nx - 1].data[age][gene].cerebellum.dd_ct[housekeepGene] = Ns_2[Nx - 1].data[age][gene].cerebellum.d_ct[housekeepGene] - Ns_2[Nx - 1].data["3m"][gene].cerebellum.d_ct[housekeepGene]
    Ns_2[Nx - 1].data[age][gene].cerebellum.two_pow_ddct[housekeepGene] = pow(2.0, Ns_2[Nx - 1].data[age][gene].cerebellum.dd_ct[housekeepGene])
  if region == "ROB":
    Ns_2[Nx - 1].data[age][gene].rob.dd_ct[housekeepGene] = Ns_2[Nx - 1].data[age][gene].rob.d_ct[housekeepGene] - Ns_2[Nx - 1].data["3m"][gene].rob.d_ct[housekeepGene]
    Ns_2[Nx - 1].data[age][gene].rob.two_pow_ddct[housekeepGene] = pow(2.0, Ns_2[Nx - 1].data[age][gene].rob.dd_ct[housekeepGene])


# Calculate all Delta_Delta_CT values & 2^Delta_Delta_CT vlaues
# --> Delta_Delta_CT = (DeltaCT - 3-MONTH.HOUSEKEEPING-GENE.THIS_REGION.DeltaCT)
for hkg in housekeeping_genes:
  for g in target_genes:
    for n in [1,2,3,4]:
      for r in target_regions:
        for a in ages:
          calc_dd_ct_for_region_2(n, r, g, a, hkg)


# Export data to single CSV
def make_row_2(Nx, gene, age, region):
  row = []
  row.append(Nx)
  row.append(gene)
  row.append(age)
  row.append(region)
  if region == "Brain Total":
    row.append(Ns_2[Nx - 1].data[age][gene].brain_total.average_ct)
    row.append(Ns_2[Nx - 1].data[age][gene].brain_total.d_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].brain_total.dd_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].brain_total.two_pow_ddct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].brain_total.d_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].brain_total.dd_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].brain_total.two_pow_ddct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].brain_total.d_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].brain_total.dd_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].brain_total.two_pow_ddct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].brain_total.outlier)
  if region == "Cortex":
    row.append(Ns_2[Nx - 1].data[age][gene].cortex.average_ct)
    row.append(Ns_2[Nx - 1].data[age][gene].cortex.d_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].cortex.dd_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].cortex.two_pow_ddct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].cortex.d_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].cortex.dd_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].cortex.two_pow_ddct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].cortex.d_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].cortex.dd_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].cortex.two_pow_ddct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].cortex.outlier)
  if region == "Hippocampus":
    row.append(Ns_2[Nx - 1].data[age][gene].hippocampus.average_ct)
    row.append(Ns_2[Nx - 1].data[age][gene].hippocampus.d_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].hippocampus.dd_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].hippocampus.two_pow_ddct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].hippocampus.d_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].hippocampus.dd_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].hippocampus.two_pow_ddct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].hippocampus.d_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].hippocampus.dd_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].hippocampus.two_pow_ddct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].hippocampus.outlier)
  if region == "Striatum":
    row.append(Ns_2[Nx - 1].data[age][gene].striatum.average_ct)
    row.append(Ns_2[Nx - 1].data[age][gene].striatum.d_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].striatum.dd_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].striatum.two_pow_ddct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].striatum.d_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].striatum.dd_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].striatum.two_pow_ddct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].striatum.d_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].striatum.dd_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].striatum.two_pow_ddct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].striatum.outlier)
  if region == "Cerebellum":
    row.append(Ns_2[Nx - 1].data[age][gene].cerebellum.average_ct)
    row.append(Ns_2[Nx - 1].data[age][gene].cerebellum.d_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].cerebellum.dd_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].cerebellum.two_pow_ddct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].cerebellum.d_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].cerebellum.dd_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].cerebellum.two_pow_ddct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].cerebellum.d_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].cerebellum.dd_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].cerebellum.two_pow_ddct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].cerebellum.outlier)
  if region == "ROB":
    row.append(Ns_2[Nx - 1].data[age][gene].rob.average_ct)
    row.append(Ns_2[Nx - 1].data[age][gene].rob.d_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].rob.dd_ct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].rob.two_pow_ddct["RPL13"])
    row.append(Ns_2[Nx - 1].data[age][gene].rob.d_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].rob.dd_ct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].rob.two_pow_ddct["TBP"])
    row.append(Ns_2[Nx - 1].data[age][gene].rob.d_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].rob.dd_ct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].rob.two_pow_ddct["Ube2d2a"])
    row.append(Ns_2[Nx - 1].data[age][gene].rob.outlier)
  return row

headers = ["N Replicate", "Gene", "Age", "Region", "Average CT", "Delta CT (Norm. against RPL13)", "Delta Delta CT (Norm. against RPL13)", "2^Delta_Delta_CT (Norm. against RPL13)", "Delta CT (Norm. against TBP)", "Delta Delta CT (Norm. against TBP)", "2^Delta_Delta_CT (Norm. against TBP)", "Delta CT (Norm. against Ube2d2a)", "Delta Delta CT (Norm. against Ube2d2a)", "2^Delta_Delta_CT (Norm. against Ube2d2a)", "Outlier"]
with open(outputDir_2 + "ANALYSIS-2_FET_expression-across_ages_in_brain_raw_data.csv", mode='w') as outfile:
  writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)  # .writer function is specific for csv file writing- allows to write rows to csv files
  writer.writerow(headers)                         # this places the defined headers lists into the a new rown in the csv file - each element in headers list corrisponds to a single cell of a single new row
  for n in [1, 2, 3, 4]:
    for g in target_genes:
      for r in target_regions:
        for a in ages:
          writer.writerow(make_row_2(n, g, a, r))


# Create plots - FUS=Blue, EWS=Greem, TAF15=Orange
#              - One plot for each gene-region-houseKeeping combination
#              - Plot scatter points for 2^ddCT of each N
#              - Plot mean of 2^ddCT across Ns_2
#              - Plot standard deviation using N 2^ddCT values and calculated means
#                         - sd = sqrt( ( sum( abs(x - mu)^2 ) ) / N )     where x is 2^ddCT from each N and mu is the mean
#                         - Plotted as error bars as Standard-Deviation
#              - Calculate p value for each non-Brain-Total region and plot:
#                         - If p < 0.05 then put * above region's bar
#                         - If p < 0.01 then put ** above region's bar
#                         - If p < 0.001 then put *** above region's bar
#                         - Also print actual p-value above region's bar to 4dp
#                         - p-value = two-tailed one-sample T-test


for hkg in housekeeping_genes:
  for g in FET_genes:
    for r in target_regions:
      data_points_ddCT = {
        "3m" : [],
        "12m" : [],
        "24m" : []
      }
      data_points_2pow_ddCT = {
        "3m" : [],
        "12m" : [],
        "24m" : []
      }
      twoPow_ddct_mean = {
        "3m" : 0.0,
        "12m" : 0.0,
        "24m" : 0.0
      }
      twoPow_ddct_sd = {
        "3m" : 0.0,
        "12m" : 0.0,
        "24m" : 0.0
      }
      significance_data = {
        "3m" : {
          "p_value" : 0.0,
          "p_bonferroni" : 0.0,
          "sig" : ""
        },
        "12m" : {
          "p_value" : 0.0,
          "p_bonferroni" : 0.0,
          "sig" : ""
        },
        "24m" : {
          "p_value" : 0.0,
          "p_bonferroni" : 0.0,
          "sig" : ""
        }
      }
      # !!!!!!! OUTLIERS IGNORED !!!!!!!!
      for n in [1,2,3,4]:
        for a in ages:
          if r == "Brain Total":
            data_points_ddCT[a].append(Ns_2[n-1].data[a][g].brain_total.dd_ct[hkg])
            data_points_2pow_ddCT[a].append(Ns_2[n-1].data[a][g].brain_total.two_pow_ddct[hkg])
          if r == "Cortex":
            data_points_ddCT[a].append(Ns_2[n-1].data[a][g].cortex.dd_ct[hkg])
            data_points_2pow_ddCT[a].append(Ns_2[n-1].data[a][g].cortex.two_pow_ddct[hkg])
          if r == "Hippocampus":
            data_points_ddCT[a].append(Ns_2[n-1].data[a][g].hippocampus.dd_ct[hkg])
            data_points_2pow_ddCT[a].append(Ns_2[n-1].data[a][g].hippocampus.two_pow_ddct[hkg])
          if r == "Striatum":
            data_points_ddCT[a].append(Ns_2[n-1].data[a][g].striatum.dd_ct[hkg])
            data_points_2pow_ddCT[a].append(Ns_2[n-1].data[a][g].striatum.two_pow_ddct[hkg])
          if r == "Cerebellum":
            data_points_ddCT[a].append(Ns_2[n-1].data[a][g].cerebellum.dd_ct[hkg])
            data_points_2pow_ddCT[a].append(Ns_2[n-1].data[a][g].cerebellum.two_pow_ddct[hkg])
          if r == "ROB":
            data_points_ddCT[a].append(Ns_2[n-1].data[a][g].rob.dd_ct[hkg])
            data_points_2pow_ddCT[a].append(Ns_2[n-1].data[a][g].rob.two_pow_ddct[hkg])
      # Calculate means + standard deviations
      for a in ages:
        # Remove outliers : possible confidence levels = {"q90", "q95", "q99"}
        # data_points_2pow_ddCT[a] = remove_outliers(data_points_2pow_ddCT[a], "q95")
        # data_points_ddCT[a] = remove_outliers(data_points_ddCT[a], "q95")
        twoPow_ddct_mean[a] = sum(data_points_2pow_ddCT[a]) / len(data_points_2pow_ddCT[a])
        twoPow_ddct_sd[a] = stat.stdev(data_points_2pow_ddCT[a])
        # Calculate significance (p-value) via One Sample T Test  
        if a == "3m": # Ignore 3m
          continue
        SD_ddct = stat.stdev(data_points_ddCT[a])
        num_replicates = len(data_points_ddCT[a])
        std_err_mean = SD_ddct / math.sqrt(num_replicates)
        x_bar = sum(data_points_ddCT["3m"]) / len(data_points_ddCT["3m"])
        alt_hypothesis = sum(data_points_ddCT[a]) / len(data_points_ddCT[a]) # M
        degs_of_freedom = len(ages) - 1
        num_tails = 1
        t_val = (x_bar - alt_hypothesis) / std_err_mean
        p_val = (1 - scipy.stats.t.cdf(x=abs(t_val), df=degs_of_freedom)) * num_tails
        p_bonf = p_val * len(ages) # IGNORE
        significance_data[a]["p_value"] = p_val
        significance_data[a]["p_bonferroni"] = p_bonf
        sig_str = ""
        if p_val < 0.001:
          sig_str = "***"
        elif p_val < 0.01:
          sig_str = "**"
        elif p_val < 0.05:
          sig_str = "*"
        significance_data[a]["sig"] = sig_str



      w = 0.5 # bar width #### NEW BAR WIDTH
      x = list(range(1,len(ages)+1)) # create list so each region has a x coord

      colour = ["#2b4357","#99b4cb","#D3D3D3"] # Colour (blue) for FUS ##NEW TWO COLOURS
      if g == "TAF15" :
        colour = ["#414d1e","#c7d98b","#D3D3D3"] # orange    ##NEW TWO COLOURS
      if g == "EWS":
        colour = ["#385749","#8ed6b7","#D3D3D3"] # green     ##NEW TWO COLOURS 


      fig, ax = plt.subplots()
      tickLabels_ages = ["3 Months", "12 Months", "24 Months"]
      ax.bar(x, 
             height=[twoPow_ddct_mean["3m"], twoPow_ddct_mean["12m"], twoPow_ddct_mean["24m"]],
             yerr = [twoPow_ddct_sd["3m"], twoPow_ddct_sd["12m"], twoPow_ddct_sd["24m"]],
             capsize = 12, # error bar cap width in points
             width = w,
             tick_label = tickLabels_ages,
             color = colour
             )
      
      for i in range(len(x)):
        # distribute scatter randomly across whole width of bar
        bar_x_coords = np.linspace(x[i] - (w/2.5), x[i] + w/2.5, len(data_points_2pow_ddCT[ages[i]]))
        ax.scatter(bar_x_coords, data_points_2pow_ddCT[ages[i]], color="black", s=15)

      # Do some editing of aesthetics
      ax.spines[['right', 'top']].set_visible(False)
      ax.yaxis.set_minor_locator(AutoMinorLocator())
      ax.tick_params(which='major', length=7)
      ax.tick_params(which='minor', length=4, color='black')
      ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

      # Add comparison caps
      for i in range(len(x)):
        if i == 0:
          continue
        region = ages[i]
        # if any stars, print error bar between brain total and current region
        if significance_data[region]["sig"] != "":
          xStart = x[0]
          xEnd = x[i]
          y = ax.get_ylim()[-1]
          xs = [xStart, xStart, xEnd, xEnd]
          ys = [y, y*1.1, y*1.1, y]
          plt.plot(xs, ys, lw=1, color='black')
          plt.text(xStart + ((xEnd - xStart)/2), y*1.1, significance_data[region]["sig"], ha='center', va='bottom', fontsize=12) ######## NEW FONTSIZE

      # Add significance + p_value under bars (exclude 3m (i.e. age 1))
      extraArt = []   ##### NEW
      for i in range(len(x)):
        if i == 0:
          continue
        height = ax.get_ylim()
        region = ages[i]
        #               plt.text(x[i], height[-1]*1.05, "p = " + str(round(significance_data[region]["p_value"], 4)), ha='center', va='bottom', fontsize=12)
        extraArt.append(plt.text(x[i], height[-1]*1.05, "p = " + str(round(significance_data[region]["p_value"], 4)), ha='center', va='bottom', fontsize=8))     ##### NEW


      #               plt.ylabel(g + " mRNA Levels Relative to\n3 months Normalised to " + hkg, fontsize=14)
      extraArt.append(plt.ylabel(g + " mRNA Levels Relative to\n3 months Normalised to " + hkg, fontsize=12))   ##### NEW
      #               plt.xlabel(r, fontsize=14)
      extraArt.append(plt.xlabel(r, fontsize=14))        ##### NEW
      plt.axhline(y=0, linestyle='-', color='black', lw=0.5)
      plt.tight_layout()
      fig.set_size_inches(3,4)                    ##### NEW
    # plt.savefig(outputDir_2 + hkg + "_" + g + "_" + r + "_graph.png")
      plt.savefig(outputDir_2 + hkg + "_" + g + "_" + r + "_graph.png", dpi=500, bbox_inches='tight', bbox_extra_artists=tuple(extraArt))     ##### NEW
      plt.close()

#------------------------------------------------------------------------------------

#-------------- Analysis 3 - Spinal Cord FET Expression Across Ages  ----------------



 # N_number_data {
 #      3_month_data {
 #          RPL13 {
 #              Spinal_cord {
 #                  Average_CT
 #                  Delta_CT (AverageCT - HOUSEKEEPING-GENE.THIS_REGION.AverageCT)
 #                  Delta_Delta_CT (DeltaCT - 3-MONTH.HOUSEKEEPING-GENE.THIS_REGION.DeltaCT)
 #                  2_^_Delta_Delta_CT (2 ^ Delta_Delta_CT)                 <------- PLOT ALL THESE ONES
 #              }
 #              
 #          }
 #          FUS
 #          EWS
 #          TAF15
 #      }
 #      12_month_data {
 # 
 #      }
 #      24_month_data {
 # 
 #      }
 # }
class region_values_3:
    def __init__(self, region):
      self.region_name = region
      self.average_ct = 0.0
      self.d_ct = {
        "RPL13" : 0.0,
        "TBP" : 0.0,
        "Ube2d2a" : 0.0
      }
      self.dd_ct = {
        "RPL13" : 0.0,
        "TBP" : 0.0,
        "Ube2d2a" : 0.0
      }
      self.two_pow_ddct = {
        "RPL13" : 0.0,
        "TBP" : 0.0,
        "Ube2d2a" : 0.0
      }
      self.outlier = False

class gene_regions_3:
    def __init__(self, gene):
       self.gene_name = gene
       self.spinal_cord = region_values_3("Spinal Cord")

class N_Num_Data_3:
    def __init__(self):
        self.data = {
          "3m" : {
            "RPL13" : gene_regions_3("RPL13"),
            "FUS" : gene_regions_3("FUS"),
            "EWS" : gene_regions_3("EWS"),
            "TAF15" : gene_regions_3("TAF15"),
            "TBP" : gene_regions_3("TBP"),
            "Ube2d2a" : gene_regions_3("Ube2d2a")
          },
          "12m" : {
            "RPL13" : gene_regions_3("RPL13"),
            "FUS" : gene_regions_3("FUS"),
            "EWS" : gene_regions_3("EWS"),
            "TAF15" : gene_regions_3("TAF15"),
            "TBP" : gene_regions_3("TBP"),
            "Ube2d2a" : gene_regions_3("Ube2d2a")
          },
          "24m" : {
            "RPL13" : gene_regions_3("RPL13"),
            "FUS" : gene_regions_3("FUS"),
            "EWS" : gene_regions_3("EWS"),
            "TAF15" : gene_regions_3("TAF15"),
            "TBP" : gene_regions_3("TBP"),
            "Ube2d2a" : gene_regions_3("Ube2d2a")
          }
        }


Ns_3 = [N_Num_Data_3(), N_Num_Data_3(), N_Num_Data_3(), N_Num_Data_3(), N_Num_Data_3(), N_Num_Data_3(), N_Num_Data_3(), N_Num_Data_3()]
#    [ N1 , N2 , N3 , N4 , N5, N6, N7, N8]

outputDir_3 = "ANALYSIS-3_Spinal_Cord_FET_expression_Across_Ages/"
# make sure outputDir_3 exists
if not os.path.exists(outputDir_3):
  os.makedirs(outputDir_3)

target_genes = ["FUS", "EWS", "TAF15", "RPL13", "TBP", "Ube2d2a"]
target_regions = ["Spinal Cord"]
housekeeping_genes = ["RPL13", "TBP", "Ube2d2a"]
FET_genes = ["FUS", "EWS", "TAF15"]

def init_avg_ct_3(Nx, gene, age, region, record):
  # Must -1 from Nx when we use it so that it goes to the correct N
  # i.e. Ns_3[0] corresponds to N1 --> therefore a Nx=1 corresponding to N1, we must -1
  # Spinal Cord Ctrl = N1-N4
  # Spinal Cord = N5-N8
  if region == "Spinal Cord":
    Nx = Nx + 4
  if region == "Spinal Cord" or region == "Spinal Cord Ctrl":
    Ns_3[Nx - 1].data[age][gene].spinal_cord.average_ct = record.average_ct
    Ns_3[Nx - 1].data[age][gene].spinal_cord.outlier = record.outlier


# Initialise all average_ct values
for g in target_genes:
   for r1 in extractedData_N1[g]:
      init_avg_ct_3(1, g, r1.age, r1.region, r1)
   for r2 in extractedData_N2[g]:
      init_avg_ct_3(2, g, r2.age, r2.region, r2)
   for r3 in extractedData_N3[g]:
      init_avg_ct_3(3, g, r3.age, r3.region, r3)
   for r4 in extractedData_N4[g]:
      init_avg_ct_3(4, g, r4.age, r4.region, r4) 


# Calculate all Delta_CT values --> Delta_CT = (AverageCT - HOUSEKEEPING-GENE.THIS_REGION.AverageCT)
for g in target_genes:
  for hkg in housekeeping_genes:
    for n in [1,2,3,4,5,6,7,8]:
      # Must -1 from n when we use it so that it goes to the correct N
      # i.e. Ns_3[0] corresponds to N1 --> therefore a n=1 corresponding to N1, we must -1
      for a in ages:
        Ns_3[n - 1].data[a][g].spinal_cord.d_ct[hkg] = Ns_3[n - 1].data[a][g].spinal_cord.average_ct - Ns_3[n - 1].data[a][hkg].spinal_cord.average_ct


# Calculate all Delta_Delta_CT values & 2^Delta_Delta_CT vlaues
# --> Delta_Delta_CT = (DeltaCT - 3-MONTH.HOUSEKEEPING-GENE.THIS_REGION.DeltaCT)
for hkg in housekeeping_genes:
  for g in target_genes:
    for n in [1,2,3,4,5,6,7,8]:
        for a in ages:
          Ns_3[n - 1].data[a][g].spinal_cord.dd_ct[hkg] = Ns_3[n - 1].data[a][g].spinal_cord.d_ct[hkg] - Ns_3[n - 1].data["3m"][g].spinal_cord.d_ct[hkg]
          Ns_3[n - 1].data[a][g].spinal_cord.two_pow_ddct[hkg] = pow(2.0, Ns_3[n - 1].data[a][g].spinal_cord.dd_ct[hkg])



headers = ["N Replicate", "Gene", "Age", "Region", "Average CT", "Delta CT (Norm. against RPL13)", "Delta Delta CT (Norm. against RPL13)", "2^Delta_Delta_CT (Norm. against RPL13)", "Delta CT (Norm. against TBP)", "Delta Delta CT (Norm. against TBP)", "2^Delta_Delta_CT (Norm. against TBP)", "Delta CT (Norm. against Ube2d2a)", "Delta Delta CT (Norm. against Ube2d2a)", "2^Delta_Delta_CT (Norm. against Ube2d2a)", "Outlier"]
with open(outputDir_3 + "ANALYSIS-2_FET_expression-across_ages_in_brain_raw_data.csv", mode='w') as outfile:
  writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)  # .writer function is specific for csv file writing- allows to write rows to csv files
  writer.writerow(headers)                         # this places the defined headers lists into the a new rown in the csv file - each element in headers list corrisponds to a single cell of a single new row
  for n in [1,2,3,4,5,6,7,8]:
    for g in target_genes:
        for a in ages:
            # Export data to single CSV
            row = []
            row.append(n)
            row.append(g)
            row.append(a)
            row.append("Spinal Cord")
            row.append(Ns_3[n - 1].data[a][g].spinal_cord.average_ct)
            row.append(Ns_3[n - 1].data[a][g].spinal_cord.d_ct["RPL13"])
            row.append(Ns_3[n - 1].data[a][g].spinal_cord.dd_ct["RPL13"])
            row.append(Ns_3[n - 1].data[a][g].spinal_cord.two_pow_ddct["RPL13"])
            row.append(Ns_3[n - 1].data[a][g].spinal_cord.d_ct["TBP"])
            row.append(Ns_3[n - 1].data[a][g].spinal_cord.dd_ct["TBP"])
            row.append(Ns_3[n - 1].data[a][g].spinal_cord.two_pow_ddct["TBP"])
            row.append(Ns_3[n - 1].data[a][g].spinal_cord.d_ct["Ube2d2a"])
            row.append(Ns_3[n - 1].data[a][g].spinal_cord.dd_ct["Ube2d2a"])
            row.append(Ns_3[n - 1].data[a][g].spinal_cord.two_pow_ddct["Ube2d2a"])
            row.append(Ns_3[n - 1].data[a][g].spinal_cord.outlier)
            writer.writerow(row)



# Create plots - FUS=Blue, EWS=Greem, TAF15=Orange
#              - One plot for each gene-region-houseKeeping combination
#              - Plot scatter points for 2^ddCT of each N
#              - Plot mean of 2^ddCT across Ns_3
#              - Plot standard deviation using N 2^ddCT values and calculated means
#                         - sd = sqrt( ( sum( abs(x - mu)^2 ) ) / N )     where x is 2^ddCT from each N and mu is the mean
#                         - Plotted as error bars as Standard-Deviation
#              - Calculate p value for each non-Brain-Total region and plot:
#                         - If p < 0.05 then put * above region's bar
#                         - If p < 0.01 then put ** above region's bar
#                         - If p < 0.001 then put *** above region's bar
#                         - Also print actual p-value above region's bar to 4dp
#                         - p-value = two-tailed one-sample T-test


for hkg in housekeeping_genes:
  for g in FET_genes:
    for r in target_regions:
      data_points_ddCT = {
        "3m" : [],
        "12m" : [],
        "24m" : []
      }
      data_points_2pow_ddCT = {
        "3m" : [],
        "12m" : [],
        "24m" : []
      }
      twoPow_ddct_mean = {
        "3m" : 0.0,
        "12m" : 0.0,
        "24m" : 0.0
      }
      twoPow_ddct_sd = {
        "3m" : 0.0,
        "12m" : 0.0,
        "24m" : 0.0
      }
      significance_data = {
        "3m" : {
          "p_value" : 0.0,
          "p_bonferroni" : 0.0,
          "sig" : ""
        },
        "12m" : {
          "p_value" : 0.0,
          "p_bonferroni" : 0.0,
          "sig" : ""
        },
        "24m" : {
          "p_value" : 0.0,
          "p_bonferroni" : 0.0,
          "sig" : ""
        }
      }
      # !!!!!!! OUTLIERS IGNORED !!!!!!!!
      for n in [1,2,3,4,5,6,7,8]:
        for a in ages:
          if r == "Spinal Cord":
            data_points_ddCT[a].append(Ns_3[n-1].data[a][g].spinal_cord.dd_ct[hkg])
            data_points_2pow_ddCT[a].append(Ns_3[n-1].data[a][g].spinal_cord.two_pow_ddct[hkg])
      # Calculate means + standard deviations
      for a in ages:
        # Remove outliers : possible confidence levels = {"q90", "q95", "q99"}
        # data_points_2pow_ddCT[a] = remove_outliers(data_points_2pow_ddCT[a], "q95")
        # data_points_ddCT[a] = remove_outliers(data_points_ddCT[a], "q95")
        twoPow_ddct_mean[a] = sum(data_points_2pow_ddCT[a]) / len(data_points_2pow_ddCT[a])
        twoPow_ddct_sd[a] = stat.stdev(data_points_2pow_ddCT[a])
        # Calculate significance (p-value) via One Sample T Test      
        if a == "3m": # Ignore 3m
          continue
        SD_ddct = stat.stdev(data_points_ddCT[a])
        num_replicates = len(data_points_ddCT[a])
        std_err_mean = SD_ddct / math.sqrt(num_replicates)
        x_bar = sum(data_points_ddCT["3m"]) / len(data_points_ddCT["3m"])
        alt_hypothesis = sum(data_points_ddCT[a]) / len(data_points_ddCT[a]) # M
        degs_of_freedom = len(ages) - 1 
        num_tails = 1
        t_val = (x_bar - alt_hypothesis) / std_err_mean
        p_val = (1 - scipy.stats.t.cdf(x=abs(t_val), df=degs_of_freedom)) * num_tails
        p_bonf = p_val * len(ages) # IGNORE
        significance_data[a]["p_value"] = p_val
        significance_data[a]["p_bonferroni"] = p_bonf
        sig_str = ""
        if p_val < 0.001:
          sig_str = "***"
        elif p_val < 0.01:
          sig_str = "**"
        elif p_val < 0.05:
          sig_str = "*"
        significance_data[a]["sig"] = sig_str



      w = 0.5 # bar width
      x = list(range(1,len(ages)+1)) # create list so each region has a x coord

      colour = ["#2b4357","#99b4cb","#D3D3D3"] # Colour (blue) for FUS ##NEW TWO COLOURS
      if g == "TAF15" :
        colour = ["#414d1e","#c7d98b","#D3D3D3"] # orange     ##NEW TWO COLOURS
      if g == "EWS":
        colour = ["#385749","#8ed6b7","#D3D3D3"] # green      ##NEW TWO COLOURS

      fig, ax = plt.subplots()
      tickLabels_ages = ["3 Months", "12 Months", "24 Months"]
      ax.bar(x, 
             height=[twoPow_ddct_mean["3m"], twoPow_ddct_mean["12m"], twoPow_ddct_mean["24m"]],
             yerr = [twoPow_ddct_sd["3m"], twoPow_ddct_sd["12m"], twoPow_ddct_sd["24m"]],
             capsize = 12, # error bar cap width in points
             width = w,
             tick_label = tickLabels_ages,
             color = colour
             )
      
      for i in range(len(x)):
        # distribute scatter randomly across whole width of bar
        bar_x_coords = np.linspace(x[i] - (w/2.5), x[i] + w/2.5, len(data_points_2pow_ddCT[ages[i]]))
        ax.scatter(bar_x_coords, data_points_2pow_ddCT[ages[i]], color="black", s=15)

      # Do some editing of aesthetics
      ax.spines[['right', 'top']].set_visible(False)
      ax.yaxis.set_minor_locator(AutoMinorLocator())
      ax.tick_params(which='major', length=7)
      ax.tick_params(which='minor', length=4, color='black')
      ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

      # Add comparison caps
      for i in range(len(x)):
        if i == 0:
          continue
        region = ages[i]
        # if any stars, print error bar between brain total and current region
        if significance_data[region]["sig"] != "":
          xStart = x[0]
          xEnd = x[i]
          y = ax.get_ylim()[-1]
          xs = [xStart, xStart, xEnd, xEnd]
          ys = [y, y*1.1, y*1.1, y]
          plt.plot(xs, ys, lw=1, color='black')
          plt.text(xStart + ((xEnd - xStart)/2), y*1.1, significance_data[region]["sig"], ha='center', va='bottom', fontsize=12)

      # Add significance + p_value under bars (exclude 3m (i.e. age 1))
      extraArt = []   ##### NEW
      for i in range(len(x)):
        if i == 0:
          continue
        height = ax.get_ylim()
        region = ages[i]
        # NEW           plt.text(x[i], height[-1]*1.05, "p = " + str(round(significance_data[region]["p_value"], 4)), ha='center', va='bottom', fontsize=12))
        extraArt.append(plt.text(x[i], height[-1]*1.05, "p = " + str(round(significance_data[region]["p_value"], 4)), ha='center', va='bottom', fontsize=8)) ##### NEW

      # NEW            plt.ylabel(g + " mRNA Levels Relative to\n3 months Normalised to " + hkg, fontsize=14)
      extraArt.append (plt.ylabel(g + " mRNA Levels Relative to\n3 months Normalised to " + hkg, fontsize=12))
      # NEW            plt.xlabel(r, fontsize=14)
      extraArt.append (plt.xlabel(r, fontsize=14))
      plt.axhline(y=0, linestyle='-', color='black', lw=0.5)
      plt.tight_layout()
      fig.set_size_inches(3,4) ##### NEW 
     #plt.savefig(outputDir_3 + hkg + "_" + g + "_" + r + "_graph.png") 
      plt.savefig(outputDir_3 + hkg + "_" + g + "_" + r + "_graph.png", dpi=500, bbox_inches='tight', bbox_extra_artists=tuple(extraArt)) ###### NEW
      plt.close()

#------------------------------------------------------------------------------------
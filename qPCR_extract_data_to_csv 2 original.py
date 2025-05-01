import csv

# ---------------------------------------------- Change Below for Each Analysis ----------------------------------------------------------- #

# ADD BELOW THE FULL PATH TO THE CSV FILE DIRECTORY
#   - Makes it easier to access the csv data files

csvDir = "~/Desktop/PhD/Experiments/RT - qPCR/Taqman assays C57BL6/FET in aging_C57BL6 Mice/N1/"
# csvDir = "/Users/fw17231/Desktop/Laura_qPCR_analysis/"

# ADD BELOW ALL ROTOR-GENE qPCR GENERATED CSV FILES
#   - Files containing the data which we want to use 
csvFiles = ["-RT RPL13 in 3m_12m_24mTakyone Probe TaqMan assay 2023-06-13 (1).csv", 
            "EWS in 3m_12m_24mTakyone Probe TaqMan assay 2023-06-13 (1).csv", 
            "FUS in 3m_12m_24mTakyone Probe TaqMan assay 2023-06-13 (1).csv", 
            "RPL13 in 3m_12m_24mTakyone Probe TaqMan assay 2023-06-13 (1).csv", 
            "TAF15 in 3m_12m_24mTakyone Probe TaqMan assay 2023-06-13 (1).csv", 
            "TBP in 3m_12m_24mTakyone Probe TaqMan assay 2023-06-13 (1).csv", 
            "Ube2d2a in 3m_12m_24mTakyone Probe TaqMan assay 2023-06-13 (1).csv"]


# ADD OUTPUT DIRECTORY PATH BELOW
#   - Filename of data that will be save | # naming file into working directory
output_Filename = "processed_N1_data.csv"      


# ----------------------------------------------------------------------------------------------------------------------------------------- #

# Define a class (data structure) for each age-region gene expression CT values
# Telling my coputer that I will be creating a data set in this structure. It's like a bluprint.
class age_region_data:
  def __init__(self, age, region):
    self.age = age
    self.region = region
    self.ct = []                                                  #list of CTs?
    self.average_ct = 0.0
    self.outlier = False # all data points are more than 0.5_delta from each other will be labelled as false


# Declare common variables that will be needed
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
    
# Define a key-value structure. Each gene is a key, with its associated value being many (a list of) `age_region_data`s
extractedData = {
  "-RT RPL13" : [],
  "EWS" : [],
  "FUS" : [],
  "RPL13" : [],
  "TAF15" : [],
  "TBP" : [],
  "Ube2d2a" : []
}

# Function to easily get which region the csv row relates to
def getRegion(row):
  name = row[2]
  for i in regions:
    if i in name:
      return i
  return "NULL"
#-----------------------------------------------------------

# Function to easily get which age the csv row relates to
def getAge(row):
  name = row[2]
  for i in ages:
    if i in name:
      return i
  return "NULL"
#-----------------------------------------------------------

# Function to easily get which gene we are reading data for
def getGene(filename):
  for g in genes:
    if g in filename:
      return g
  return "NULL"
#-----------------------------------------------------------

# Function to easily add CT value to correct data field in `extractedData`
def addCT(gene, region, age, ct):
  # Convert CT to float
  ct_float = 0.0
  if ct != "":
    ct_float = float(ct)
  
  # See if entry for age-region already exists for target gene
  exists = False
  index = -1
  for i in range(0, len(extractedData[gene])):
    entries = extractedData[gene]
    if entries[i].age == age and entries[i].region == region:
      exists = True
      index = i
      break
  
  # If entry exists, append CT value
  if exists == True:
    extractedData[gene][index].ct.append(ct_float)
  # Else, create new entry
  else:
    newEntry = age_region_data(age, region)
    newEntry.ct.append(ct_float)
    extractedData[gene].append(newEntry)

#-----------------------------------------------------------


# Read in data from all files
for f in csvFiles:
  gene = getGene(f)
  # Check that gene was properly extracted
  if gene == "NULL":
    print("Gene not recognised in filename {n}".format(n=f))
  # Open CSV file in read mode as variable `file`
  with open(csvDir + f, mode='r') as file:
    csvreader = csv.reader(file)
    for row in csvreader:
      # Quick check to make sure row has minimum number of columns
      if len(row) < 5:
        continue
      # Get the region and age we are looking at
      reg = getRegion(row)
      age = getAge(row)
      # Check region and age could be extracted
      if(reg == "NULL" or age == "NULL"):
        continue
      # CT filed is in the 4th column (starting from 0th)
      addCT(gene, reg, age, row[4])


# Calculate average for all CT values 
# - if value is 0, don't include in average
# - if value is more than 0.5 away from other values don't include in average
for g in genes:
  for i in extractedData[g]:
    ct_data = i.ct
    ct_data.sort()
    total = ct_data[1]
    num_elems = 1
    if abs(ct_data[0]-ct_data[1]) <= 0.5:
      num_elems = num_elems + 1
      total = total + ct_data[0]
    if abs(ct_data[2]-ct_data[1]) <= 0.5:
      num_elems = num_elems + 1
      total = total + ct_data[2]
  
    if num_elems == 1:
      # All data outlier, mark entry 
      i.outlier = True
    i.average_ct = total / num_elems



for g in genes:
  print(g)
  for i in extractedData[g]:
    print("\t{a} {r} : ct = {ct}, avg_ct = {act}, outlier = {b}".format(a=i.age, r=i.region, ct=i.ct, act=i.average_ct, b=i.outlier))


# Save data to CSV file
headers = ["Gene", "Age", "Region", "ct_0", "ct_1", "ct_2", "Average CT", "Outlier"]
with open(output_Filename, mode='w') as outfile:
  writer = csv.writer(outfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)  # .writer function is specific for csv file writing- allows to write rows to csv files
  writer.writerow(headers)                         # this places the defined headers lists into the a new rown in the csv file - each element in headers list corrisponds to a single cell of a single new row
  for g in genes:
    for i in extractedData[g]:
      row = []
      row.append(g)
      row.append(i.age)
      row.append(i.region)
      row.append(i.ct[0])
      row.append(i.ct[1])
      row.append(i.ct[2])
      row.append(i.average_ct)
      row.append(i.outlier)
      writer.writerow(row)
      
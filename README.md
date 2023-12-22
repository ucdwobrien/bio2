# Bio Principles & Cellular Org - Assignment 2
### Explanation
The explanation of the code will be done on a step to step level using the instructions of the assignment as the reference for each step.


### Step 2
- The code extracts files from tar zip and sets working directory.

### Steps 3, 4 & 5
- Each dataset that is to be used fo the investigation is read in.

### Step 6
- The RNA sequencing dataset has duplicates filtered out.
- The rows are labelled by HUGO values.
- Explanatory information is removed from datasets.
- CNA & RNA datasets are filtered to have the same patients.

### Step 7
- Metadata object created and values added to indicate HER2+ patients.
- Histogram plot to show distribution.

### Step 8
- Data is normalized using DESeq.

### Step 9
- Significantly expressed genes filtered and ordered by log2fold change.
- Upregulated & downregulated gene lists created.

### Step 10
- Upregulated significantly expressed genes filtered.
- Entrez IDs of these are then retrieved from RNA dataset.

### Steps 11 & 12
- Variance stabilised transformed values calculated.
- Performed PCA, printed out variance of PCs.
- Plotting of PCA

### Step 13
- Kmeans clustering performed to cluster the data.
- PCA performed on each cluster.
- First 2 PCs of each cluster plotted.

### Steps 14 & 15
- Transposed VST values for gene dataset to be used in cox regression.
- Replaced patient IDs in patient dataset to match naming in other data.
- Filter out patients that do not match.
- Change survival status to numeric.
- Create survival object with patient survival status and survival months.
- Cox regression calculated for each gene and placed in list.
- Coeffecients of each gene extracted from list.
- Coeffecients sorted and top 5 lowest values extracted.
- Plot of survival rate between HER2+ and HER2-

---
title: "Starting with data in R"
teaching: 30
exercises: 20
questions:
- How can I use data stored as a table of values (a data frame) for analysis?
objectives:
- Load external data from a .csv file into a data frame.
- Install and load packages.
- Describe what a data frame (tibble) is.
- Summarize the contents of a data frame.
- Describe the purpose of the **`dplyr`** package.
- Select certain columns in a data frame with the **`dplyr`** function `select`.
- Extract certain rows in a data frame according to logical (boolean) conditions with the **`dplyr`** function `filter` .
- Link the output of one **`dplyr`** function to the input of another function with the 'pipe' operator `%>%`.
- Add new columns to a data frame that are functions of existing columns with `mutate`.

keypoints:
- The **`dplyr`** package provides functions for manipulating data frames, such as `select`, `filter`, `mutate`, `summarize`, `group_by`, and `count`.
- The **`tidyr`** package provides functions for reshaping data frames, such as `gather` and `spread` (see https://lessons.datacarpentry.org/R-ecology-lesson/03-dplyr.html).
---

### Loading the Hoptspot data

We are investigating single AA residue recurrently mutated sites in cancer genes, also known as mutation hotspots. These include 1,165 hotspots reported in [Chang et al. 2017] and [Chang et al. 2016], as well as 164 hotspots reported in [Bandlamudi et al., 2026]. The dataset is stored as a comma separated value (CSV) file.
Each row holds information for a AA change relative to a reference AA, and the columns represent specific attributes of the mutation, such as the gene, the position of the change, the reference amino acid, the mutated amino acid, and the frequency of the mutation in the dataset (see https://www.cancerhotspots.org/#/home):

```bash
$ mkdir workshop/data
$ cd workshop/data/
$ wget https://raw.githubusercontent.com/JenniferShelton/JenniferShelton.github.io/refs/heads/gh-pages/data/SNV-hotspots.csv
$ wget https://raw.githubusercontent.com/JenniferShelton/JenniferShelton.github.io/refs/heads/gh-pages/data/inflammation-01.csv
$ ls -thl

```

In an R session you can load the CSV file into a data frame using the `read_csv` function:

```r
library("tidyverse")

hotspots <- read_csv("~/workshop/data/SNV-hotspots.csv")
inflammation <- read_csv("~/workshop/data/inflammation-01.csv")
```

When you execute `read_csv` on a data file, it looks through the first 1000 rows
of each column and guesses its data type. For example, in this dataset,
`read_csv()` reads `Hugo_Symbol` as `chr` (a string data type), and `Mutability`
as `dbl` (or double-precision floating-point number). The column `inOncokb` is 
read as `lgl` or logical. You have the option to specify the data type for a column
manually by using the `col_types` argument in `read_csv`.

Next, we're going to learn some of the most common **`dplyr`** functions:

- `select()`: subset columns
- `filter()`: subset rows on conditions
- `mutate()`: create new columns by using information from other columns
- `group_by()` and `summarize()`: create summary statistics on grouped data
- `arrange()`: sort results
- `count()`: count discrete values



```r
dim(hotspots) # returns a vector with the number of rows in the first element, and the number of columns as the second element (the dimensions of the object)
nrow(hotspots) # returns the number of rows
ncol(hotspots) # returns the number of columns

str(hotspots) # provides the structure of the data frame, including the type of each column
summary(hotspots) # provides a summary of each column in the data frame
```

```r
head(hotspots) # shows the first 6 rows
```

Shows all columns (useful for wide tibbles)
```r
options(tibble.width = Inf)
head(hotspots, n=2)
```

Show the unique values in a specific column of the data frame, for example, the `Hugo_Symbol` or `inOncokb` column:

```r
distinct(hotspots, inOncokb)
hotspots %>% count(inOncokb)
```

Here "mutability" refers to the underlying rate/context where mutations tend to happen (like in specific DNA sequences or repair-prone regions), and the "log10 p-value" represents the statistical significance from the original recurrence calculation.


Select specific columns from the data frame:

```r
cols(hotspots)
cols <- c("Hugo_Symbol", "Reference_Amino_Acid", "Variant_Amino_Acid")
hotspots[cols]
```


Find specific rows in the data frame based on the values in the `Hugo_Symbol` column.

```r
hotspots[cols] %>%
  filter(Hugo_Symbol == "HIST1H3B")

genes <- c("HIST1H3B", "TP53")
hotspots[cols] %>%
  filter(Hugo_Symbol %in% genes)
```

Modify a column in the data frame using the `mutate()` function (note that you haven't saved your changes to the dataframe 
so they print to the screen):


```r
x <- "H:386"
sub(":.*", "", x)
hotspots %>%
  mutate(Mutability_Var_AA = sub(":.*", "", Variant_Amino_Acid))
```

Save a new dataframe with your customized (mutated) column.

```r
hotspots <- hotspots %>%
  mutate(Variant_Amino_Acid = sub(":.*", "", Variant_Amino_Acid)) %>%
  mutate(Reference_Amino_Acid = sub(":.*", "", Reference_Amino_Acid))
```

Additional R lessons are available at: https://lessons.datacarpentry.org/R-ecology-lesson/01-intro-to-r.html








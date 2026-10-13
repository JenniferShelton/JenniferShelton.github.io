---
title: "Manipulating, analyzing and exporting data with tidyverse"
teaching: 30
exercises: 20
questions:
- How can I use data stored as a table of values (a data frame) for analysis?
objectives:
- Select certain columns in a data frame with the **`dplyr`** functions.
- Produce scatter plots.
- Set universal plot settings.
- Modify the aesthetics of an existing ggplot plot (including axis size and color).
- Build complex and customized plots from data in a data frame.

keypoints:
- The **`dplyr`** package provides functions for manipulating data frames, such as `select`, `filter`, `mutate`, `summarize`, `group_by`, and `count`.
- The **`tidyr`** package provides functions for reshaping data frames, such as `gather` and `spread` (see https://lessons.datacarpentry.org/R-ecology-lesson/03-dplyr.html).
---

### Data manipulation using dplyr and tidyr

R bracket subsetting for row and column is handy, but it can be cumbersome and difficult to read,
especially for complicated operations. Enter **`dplyr`**. **`dplyr`** is a package for
helping with tabular data manipulation. It pairs nicely with **`tidyr`** which enables you to swiftly
convert between different data formats for plotting and analysis.

```r
library("tidyverse")
library("plotly")
```

The **`tidyverse`** package is an
"umbrella-package" that installs **`tidyr`**, **`dplyr`**, and several other useful packages for 
data analysis, such as  **`ggplot2`**, **`tibble`**, etc.

The **`tidyverse`** package tries to address 3 common issues that arise when
doing data analysis in R:

1. The results from a base R function sometimes depend on the type of data.
2. R expressions are used in a non standard way, which can be confusing for new
  learners.
3. The existence of hidden arguments having default operations that new learners are not aware
  of.

If you are working in an environment where the **`tidyverse`** package is not already installed, you can type `install.packages("tidyverse")` straight into the console. Then, type `library(tidyverse)` to load the package.

## What are **`dplyr`** and **`tidyr`**?

The package **`dplyr`** provides helper tools for the most common data manipulation
tasks. It is built to work directly with data frames, with many common tasks
optimized by being written in a compiled language (C++). An additional feature is the
ability to work directly with data stored in an external database. The benefits of
doing this are that the data can be managed natively in a relational database,
queries can be conducted on that database, and only the results of the query are
returned.

This addresses a common problem with R in that all operations are conducted
in-memory and thus the amount of data you can work with is limited by available
memory. The database connections essentially remove that limitation in that you
can connect to a database of many hundreds of GB, conduct queries on it directly, and pull
back into R only what you need for analysis.

## Plotting with **`ggplot2`**

**`ggplot2`** is a plotting package that provides helpful commands to create complex plots
from data in a data frame. It provides a more programmatic interface for
specifying what variables to plot, how they are displayed, and general visual
properties. Therefore, we only need minimal changes if the underlying data
change or if we decide to change from a bar plot to a scatterplot. This helps in
creating publication quality plots with minimal amounts of adjustments and
tweaking.

**`ggplot2`** refers to the name of the package itself. When using the package we use the
function **`ggplot()`** to generate the plots, and so references to using the function will
be referred to as **`ggplot()`** and the package as a whole as **`ggplot2`**

**`ggplot2`** plots work best with data in the 'long' format, i.e., a column for every variable,
and a row for every observation. Well-structured data will save you lots of time
when making figures with **`ggplot2`**

ggplot graphics are built layer by layer by adding new elements. Adding layers in
this fashion allows for extensive flexibility and customization of plots.

```r
ggplot(data = hotspots)
```

### Here we called ggplot and told it what data we want to show on our figure. This is not enough information for ggplot to actually draw anything. 
It only creates a blank slate for other elements to be added to.

```r
ggplot(data = hotspots,  mapping = aes(x = Mutability, y = log10_pvalue))
```

### Here we told ggplot we want to plot the “Mutability” column of the hotspots data frame on the x-axis, 
and the “log10_pvalue” column on the y-axis. Notice that we didn’t need to explicitly pass `aes` 
these columns (e.g. x = hotspots[, "Mutability"]), this is because ggplot is smart enough to 
know to look in the data for that column!

### The final part of making our plot is to tell ggplot how we want to visually represent the data. 
We do this by adding a new layer to the plot using one of the geom functions.

```r
ggplot(data = hotspots,  mapping = aes(x = Mutability, y = log10_pvalue)) +
    geom_point()
```

- add 'geoms' – graphical representations of the data in the plot (points,
  lines, bars). **`ggplot2`** offers many different geoms; we will use some
  common ones today, including:
  
  - `geom_point()` for scatter plots, dot plots, etc.
  - `geom_boxplot()` for, well, boxplots!
  - `geom_line()` for trend lines, time series, etc.

To add a geom to the plot use `+` operator. Because we have two continuous
variables, let's use `geom_point()` to create a scatter plot first:

```r
ggplot(data = hotspots,  mapping = aes(x = Mutability, y = log10_pvalue, size = Mutation_Count, color=inOncokb)) +
    geom_point()
```
Here we added marker size and color to `aes`.

Now let's add a theme:

```r
ggplot(data = hotspots,  mapping = aes(x = Mutability, y = log10_pvalue, size = Mutation_Count, color=inOncokb)) +
    geom_point() +
    theme_bw()
```

**Notes**

- Anything you put in the `ggplot()` function can be seen by any geom layers
  that you add (i.e., these are universal plot settings). This includes the x-
  and y-axis you set up in `aes()`.
- You can also specify aesthetics for a given geom independently of the
  aesthetics defined globally in the `ggplot()` function.
- The `+` sign used to add layers must be placed at the end of each line
  containing a layer. If, instead, the `+` sign is added in the line before the
  other layer, **`ggplot2`** will not add the new layer and will return an error
  message.
- You may notice that we sometimes reference 'ggplot2' and sometimes 'ggplot'.
  To clarify, 'ggplot2' is the name of the most recent version of the package.
  However, any time we call the function itself, it's just called 'ggplot'.
- The previous version of the **`ggplot2`** package, called **`ggplot`**,
  which also contained the `ggplot()` function is now unsupported and has
  been removed from CRAN in order to reduce accidental installations
  and further confusion.


You can save a figure created with ggplot and pass it to ggplotly to add/customize hover text (to see exact values for a point)
```r
p <- ggplot(data = hotspots,  mapping = aes(x = Mutability, y = log10_pvalue, size = Mutation_Count, color=inOncokb)) +
  geom_point(aes(text = paste0("Hugo_Symbol: ", Hugo_Symbol, "<br>",
                               "Mutability: ", Mutability, "<br>",
                               "log10 p-value: ", log10_pvalue, "<br>",
                               "Mutation Count: ", Mutation_Count, "<br>"
  )))  +
  theme_bw()

ggplotly(p, tooltip = "text")
```

## Exporting data

Now that you have learned how to use **`dplyr`** to extract information from
or summarize your raw data, you may want to export these new data sets to share
them with your collaborators or for archival.

Similar to the `read_csv()` function used for reading CSV files into R, there is
a `write_csv()` function that generates CSV files from data frames.

Before using `write_csv()`, we are going to create a new folder, `output`,
in our working directory (if it does not already exist) that will store this generated dataset. We don't want
to write generated datasets in the same directory as our raw data. It's good
practice to keep them separate. The `data` folder should only contain the raw,
unaltered data, and should be left alone to make sure we don't delete or modify
it. In contrast, our script will generate the contents of the `output`
directory, so even if the files it contains are deleted, we can always
re-generate them.

In preparation for our next lesson on plotting, we are going to prepare a
cleaned up version of the data set that doesn't include any missing data.

Let's start by removing observations of animals for which `Hugo_Symbol` is missing:

```{r, purl=FALSE}

hotspots <- hotspots %>%
  mutate(Variant_Amino_Acid = sub(":.*", "", Variant_Amino_Acid)) %>%
  mutate(Reference_Amino_Acid = sub(":.*", "", Reference_Amino_Acid)) %>%
  filter(!is.na(Hugo_Symbol))

```

Because we are interested in plotting recurrent mutations we may want to set a minimum  
`Mutation_Count`:

```{r, purl=FALSE}
## Only keep the most common mutations
hotspots <- hotspots %>%
    filter(Mutation_Count >= 50)
```


Now that our data set is ready, we can save it as a CSV file in our `data`
folder.

```{r, purl=FALSE, eval=FALSE}
write_csv(hotspots, file = "~/workshop/output/hotspot_results.csv")
```

```{r, purl=FALSE, eval=TRUE, echo=FALSE}
if (!dir.exists("~/workshop/output/")) dir.create("~/workshop/output/")
write_csv(hotspots, file = "~workshop/output/hotspot_results.csv")
```




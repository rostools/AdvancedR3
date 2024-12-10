# AdvancedR3: Learning R in a course

This project contains code and files from an R course.

# Brief description of folder and file contents

The following folders contain:

-   `data/`: The data ready for an analysis
-   `doc/`: Quarto files
-   `R/`: R scripts.

# Installing project R package dependencies

If dependencies have been managed by using
`usethis::use_package("packagename")` through the `DESCRIPTION` file,
installing dependencies is as easy as opening the `AdvancedR3.Rproj`
file and running this command in the console:

```         
# install.packages("remotes")
remotes::install_deps()
```

You'll need to have remotes installed for this to work.

# Resource

For more information on this folder and file workflow and setup, check
out the [prodigenr](https://rostools.github.io/prodigenr) online
documentation.


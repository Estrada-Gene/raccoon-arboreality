# Create the GitHub repository from RStudio

The preparation utility creates a local project; it does not create or upload a
GitHub repository. The following commands perform those later steps from your
own account.

1. Open `raccoon-arboreality.Rproj` and run `source("scripts/00_setup.R")`.
2. Install Git if needed (https://git-scm.com/downloads), then restart RStudio.
   Check Tools > Global Options > Git/SVN if the Git executable is not detected.
3. Configure your Git identity and credentials:

```r
install.packages(c("usethis", "gitcreds"))
usethis::use_git_config(user.name = "Gene Estrada", user.email = "YOUR_GITHUB_EMAIL")
usethis::git_default_branch_configure("main")
usethis::create_github_token(scopes = "repo", description = "Raccoon repository")
gitcreds::gitcreds_set()
```

Generate the token in the browser and enter it at the `gitcreds_set()` prompt;
it belongs in the credential store, not in an R script. If working GitHub
credentials are already configured, reuse them.

4. Initialize the local repository:

```r
usethis::use_git(message = "Add raccoon analysis scripts and input data")
```

Accept the initial commit when prompted. RStudio may reopen the project. If you
make any later edits, commit them through the Git pane before the next command.

5. Create the GitHub repository and push the committed files:

```r
usethis::use_github(private = TRUE, protocol = "https")
```

The repository name is taken from the local project folder. Use `private = FALSE`
instead if the files are ready to be public. An existing remote named `origin`
or an already-created repository with the same name requires a different setup;
this walkthrough assumes a new local project and a new remote repository.

6. On GitHub, confirm that the scripts, input data, README, and metadata are
present, while generated model files and the large download are absent. Add a
code license, replace the repository URL placeholder, and make the repository
public when appropriate for the manuscript.

For later changes: save, review the diff in RStudio's Git pane, stage the intended
files, Commit with a short description, then Push. Generated results are ignored.
To demonstrate reproducibility, clone the repository into a separate folder and
run scripts 00-05 there; an RStudio restart alone does not check for undeclared
dependencies on files in the original working directory.

Documentation:

- https://docs.posit.co/ide/user/ide/guide/tools/version-control.html
- https://usethis.r-lib.org/articles/articles/git-credentials.html
- https://usethis.r-lib.org/reference/use_git.html
- https://usethis.r-lib.org/reference/use_github.html
- https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github


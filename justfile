#!/usr/bin/env just
# shellcheck shell=bash

package_name := 'BLASTr'

# github_org := 'heronoh'

@default:
  just --choose

@test:
  #!/usr/bin/env -vS bash -i
  \builtin set -euxo pipefail;
  R -q -e 'devtools::load_all();styler::style_pkg();';
  R -q -e 'devtools::load_all();usethis::use_tidy_description();';
  R -q -e 'devtools::load_all();devtools::document();';
  R -q -e 'devtools::load_all();devtools::run_examples();';
  R -q -e 'devtools::load_all();devtools::test();';
  R -q -e 'devtools::load_all();if(file.exists("README.Rmd"))rmarkdown::render("README.Rmd", encoding = "UTF-8")';
  just check;

@test-all-examples:
  #!/usr/bin/env -vS bash -i
  \builtin set -euxo pipefail;
  R -q -e 'devtools::load_all();devtools::run_examples(run_dontrun = TRUE, run_donttest = TRUE);';

@check:
  #!/usr/bin/env -vS bash -i
  \builtin set -euxo pipefail;
  R -q -e 'rcmdcheck::rcmdcheck(args = c("--as-cran"), repos = c(CRAN = "https://cloud.r-project.org"));';

# Use R package version on the Description file to tag latest commit of the git repo
@git-tag:
  #!/usr/bin/env -vS bash -i
  \builtin set -euxo pipefail;
  __r_pkg_version="$(R -q --no-echo --silent -e 'suppressMessages({pkgload::load_all()});cat(as.character(utils::packageVersion("{{ package_name }}")));')";
  \builtin echo -ne "Tagging version: ${__r_pkg_version}\n";
  git tag -a "v${__r_pkg_version}" HEAD -m "Version ${__r_pkg_version} released";
  git push --tags;

# Things to run before releasing a new version
@pre-release:
  #!/usr/bin/env -vS bash -i
  \builtin set -euxo pipefail;
  R -q -e 'urlchecker::url_check()';
  R -q -e 'devtools::build_readme()';
  R -q -e 'withr::with_options(list(repos = c(CRAN = "https://cloud.r-project.org")), {devtools::check(remote = TRUE, manual = TRUE)})';
  R -q -e 'devtools::check_win_devel()';
  # Update CRAN comments
  # usethis::use_version('patch')
  # devtools::submit_cran()
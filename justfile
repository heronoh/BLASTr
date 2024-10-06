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
  R -q -e 'devtools::load_all();devtools::document();';
  R -q -e 'devtools::load_all();devtools::test();';

@check:
  #!/usr/bin/env -vS bash -i
  \builtin set -euxo pipefail;
  R -q -e 'rcmdcheck::rcmdcheck(args="--as-cran");';

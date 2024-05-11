1. test if unit and system tests pass: `make test; make integration_test; make pipeline_test; make system_test`
2. check if docker image works: `make build; make all`
3. update version number: `poetry version <major/minor/patch>`
4. add information to changelog on master `git log --pretty="%s"`
5. commit and push to master on github
6. create a release on github, use the tag naming convention `rel-<major>-<minor>-<patch>`, e.g. `rel-0-01-1`
   this automatically triggers a GitHub actions workflow that builds and publishes the packages to TestPyPI and PyPI
7. create a branch on github, use the naming convention `branch-<major>-<minor>-<patch>`, e.g. `branch-0-01-1`
8. merge `develop` branch into `main`

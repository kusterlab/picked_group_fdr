1. test if unit and system tests pass: `make test; make integration_test; make pipeline_test; make system_test`
2. check if docker image works: `make build; make all`
3. check if the GUI works: `python gui.py`
4. update version number: `poetry version <major/minor/patch>`
5. add information to changelog on master `git log --pretty="%s"`
6. commit and push to master on github
   this automatically triggers a GitHub actions workflow that builds and publishes the packages to TestPyPI and PyPI
7. merge `develop` branch into `main`

1. test if unit and system tests pass: `make test_for_release`
2. test if the docker image and container work: `make build; make all`
3. check if the GUI works: `python gui.py`
4. update version number: `poetry version <major/minor/patch>`
5. add information to changelog on master `git log --pretty="%s"`
6. commit and push to master on github using the commit message `Ready for v<x.x.x>`
   this automatically triggers a GitHub actions workflow that builds and publishes the packages to TestPyPI and PyPI and merges the `develop` branch into `main`

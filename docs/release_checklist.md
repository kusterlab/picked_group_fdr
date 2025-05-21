1. test if unit and system tests pass: `make test_for_release`
2. check if the GUI works: `python gui.py`
3. update version number: `poetry version <major/minor/patch>`
4. add information to changelog on master `git log --pretty="%s"`
5. commit and push to master on github using the commit message `Ready for v<x.x.x>`
   this automatically triggers a GitHub actions workflow that builds and publishes the packages to TestPyPI and PyPI and merges the `develop` branch into `main`

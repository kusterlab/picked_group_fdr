from pathlib import Path
import shutil

dist_dir = Path.cwd() / 'dist' / 'PickedGroupFDR'

lib_dir = dist_dir / 'lib'
lib_dir.mkdir(parents=True, exist_ok=True)

lib_bin_dir = dist_dir / 'Library' / 'bin'
lib_bin_dir.mkdir(parents=True, exist_ok=True)

print("Moving libraries to lib directory")
for f in dist_dir.glob('*'):
    if f.name.endswith(".egg-info'") or f.name in ["pytz", "matplotlib", "sqlalchemy", "tcl8", "PIL", "greenlet", "certifi"]:
        shutil.rmtree(f)
    #elif f.is_file() and f.name not in ['base_library.zip', 'python36.dll', 'python37.dll', 'python38.dll', 'python39.dll', 'PickedGroupFDR.exe']:
    #    f.rename(lib_dir / f.name)

print("Creating zip archive")
shutil.make_archive(Path.cwd() / 'dist' / 'PickedGroupFDR_GUI_windows', 'zip', dist_dir)

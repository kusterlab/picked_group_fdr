from pathlib import Path
import shutil

simsi_dist_dir = Path.cwd() / 'dist' / 'PickedGroupFDR'

lib_dir = simsi_dist_dir / 'lib'
lib_dir.mkdir(parents=True, exist_ok=True)

print("Moving libraries to lib directory")
for f in simsi_dist_dir.glob('*'):
    if f.name.endswith(".egg-info'") or f.name in ["pytz", "matplotlib", "sqlalchemy", "tcl8", "PIL", "greenlet", "certifi"]:
        shutil.rmtree(f)
    elif f.is_file() and f.name not in ['base_library.zip', 'python38.dll', 'python39.dll', 'SIMSI-Transfer.exe']:
        f.rename(lib_dir / f.name)

print("Creating zip archive")
shutil.make_archive(Path.cwd() / 'dist' / 'PickedGroupFDR_GUI_windows', 'zip', simsi_dist_dir)

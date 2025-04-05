from pathlib import Path
import shutil

dist_dir = Path.cwd() / 'dist' / 'PickedGroupFDR'

print("Creating zip archive")
shutil.make_archive(Path.cwd() / 'dist' / 'PickedGroupFDR_GUI_windows', 'zip', dist_dir)

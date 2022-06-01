:::: pyinstaller gui.py --noconfirm --onedir --name="PickedGroupFDR"
pyinstaller --upx-dir ../../../Downloads/upx-3.96-win64/upx-3.96-win64 --noconfirm ../PickedGroupFDR.spec && python create_gui_zip.py

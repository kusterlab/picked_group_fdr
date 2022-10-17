:::: pyinstaller gui.py --noconfirm --onedir --name="PickedGroupFDR"
pyinstaller --upx-dir ..\UPX\upx-3.96\ --noconfirm ../PickedGroupFDR.spec && python create_gui_zip.py

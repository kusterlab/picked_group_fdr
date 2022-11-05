# -*- mode: python ; coding: utf-8 -*-

import llvmlite

llvmlite_dll = os.path.join(os.path.dirname(llvmlite.__file__), 'binding', 'llvmlite.dll')

block_cipher = None

a = Analysis(['gui.py'],
             pathex=[],
             binaries=[(llvmlite_dll, '.')],
             datas=[('matplotlib_mock', 'matplotlib_mock'), ('picked_group_fdr/methods', 'picked_group-fdr/methods'), ('./pyproject.toml', '.')],
             hiddenimports=['sklearn.utils._cython_blas',
                            'sklearn.utils._typedefs',
                            'sklearn.utils._heap',
                            'sklearn.utils._sorting',
                            'sklearn.utils._vector_sentinel',
                            'sklearn.utils._weight_vector',
                            'lxml._elementpath'],
             hookspath=[],
             runtime_hooks=['add_lib.py'],
             excludes=['matplotlib'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)

a.binaries = [x for x in a.binaries if not x[0].startswith("opengl")]
a.binaries = [x for x in a.binaries if not x[0].startswith("openssl")]
a.binaries = [x for x in a.binaries if not x[0].startswith("sqlite")]
a.binaries = [x for x in a.binaries if not x[0].startswith("Qt5DBus")]
a.binaries = [x for x in a.binaries if not x[0].startswith("Qt5Network")]
a.binaries = [x for x in a.binaries if not x[0].startswith("Qt5Qml")]
a.binaries = [x for x in a.binaries if not x[0].startswith("Qt5Quick")]
a.binaries = [x for x in a.binaries if not x[0].startswith("Qt5Svg")]
a.binaries = [x for x in a.binaries if not x[0].startswith("Qt5VirtualKeyboard")]
a.binaries = [x for x in a.binaries if not x[0].startswith("Qt5WebSockets")]
a.binaries = [x for x in a.binaries if not x[0].endswith("qwebgl.dll")]
a.binaries = [x for x in a.binaries if not x[0].endswith("qofscreen.dll")]
a.binaries = [x for x in a.binaries if not x[0].endswith("qoffscreen.dll")]
a.binaries = [x for x in a.binaries if not x[0].endswith("qminimal.dll")]
a.binaries = [x for x in a.binaries if not x[0].endswith("qdirect2d.dll")]
a.binaries = [x for x in a.binaries if not "imageformats" in x[0]]
#a.binaries = [x for x in a.binaries if not "gfortran-win" in x[0]]

a.datas = [x for x in a.datas if not x[0].startswith("libcrypto")]
a.datas = [x for x in a.datas if not x[0].startswith("libssl")]
a.datas = [x for x in a.datas if not x[0].startswith("openssl")]
a.datas = [x for x in a.datas if not x[0].startswith("sqlite")]
a.datas = [x for x in a.datas if not "translations" in x[0]]

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='PickedGroupFDR',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='PickedGroupFDR')

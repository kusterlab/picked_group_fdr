import PyInstaller.utils.hooks

print("Trying to collect submodules for lxml")
PyInstaller.utils.hooks.collect_submodules('lxml')
print("Trying to collect submodules for lxml.isoschematron")
PyInstaller.utils.hooks._collect_submodules('lxml.isoschematron', 'warn once')
print("Success??")
# -*- mode: python -*-
import PyInstaller.config
PyInstaller.config.CONF['distpath'] = "..\dist"

block_cipher = None


a = Analysis(['GSadjust.py'],
             pathex=['E:\\Shared\\current\\python\\GSadjust\\main_code'],
             binaries=None,
             datas=None,
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='GSadjust',
          debug=False,
          strip=False,
          upx=True,
          console=True )

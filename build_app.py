import PyInstaller.__main__
import os

# Get the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

PyInstaller.__main__.run([
    'app.py',
    '--onefile',
    '--windowed',
    '--name=Primerly',
    '--add-data=templates;templates',
    '--add-data=static;static',
    '--icon=icon.ico',  # You'll need to add an icon file
    '--distpath=dist',
    '--workpath=build',
    '--specpath=build'
]) 
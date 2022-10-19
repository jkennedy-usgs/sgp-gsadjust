These files aren't bundled with GSadjust directly, they're compiled into the resources.py file.

Something like this on the command line:
Pyrcc5 resources.qrc -o resources.py
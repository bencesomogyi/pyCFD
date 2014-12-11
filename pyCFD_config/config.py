"""
module for setting global variables

* | use sparse the scipy.sparse module instead of numpy arrays
  | __sparse__        = False
  
* | directory for saving field files
  | __FIELDDIR__      = ""

* | directory of output files
  | __OUTDIR__        = ""

* | directory to save vtk files of iterations
  | __OUTITERDIR__    = ""

* | relative directory of vtk files within the output directory
  | __OUTITERDIRREL__ = ""

* | use IMEX schame
  | __IMEX__          = False
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

# configuration settings:
__sparse__        = False
"""use sparse the scipy.sparse module instead of numpy arrays"""
__FIELDDIR__      = ""
"""directory for saving field files"""
__OUTDIR__        = ""
"""directory of output files"""
__OUTITERDIR__    = ""
"""directory to save vtk files of iterations"""
__OUTITERDIRREL__ = ""
"""relative directory of vtk files within the output directory"""
__IMEX__          = False
"""use the IMEX schame"""
__runSettingsFile__ = "squareCylinderNonDim.conf"
"""name of run settings file"""
__runSettingsModTime__ = 0
"""modification time of run settings file"""
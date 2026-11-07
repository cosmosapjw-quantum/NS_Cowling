from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"python_frontend"))
sys.path.insert(1,str(ROOT/"pipeline"))

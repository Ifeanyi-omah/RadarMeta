import sys
from pathlib import Path

# Make the pipeline's bin/ scripts importable as modules in tests.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))

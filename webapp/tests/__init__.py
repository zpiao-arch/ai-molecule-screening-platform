"""Test package path setup for root-level test execution."""

from pathlib import Path
import sys


WEBAPP_ROOT = Path(__file__).resolve().parents[1]
if str(WEBAPP_ROOT) not in sys.path:
    sys.path.insert(0, str(WEBAPP_ROOT))

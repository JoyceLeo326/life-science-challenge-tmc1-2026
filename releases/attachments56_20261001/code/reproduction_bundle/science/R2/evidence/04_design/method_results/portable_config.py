"""Explicit replay locations; missing configuration is an error."""
import json
import os
from pathlib import Path

def required_path(name):
    config = Path(__file__).with_name("portable_paths.json")
    values = json.loads(config.read_text(encoding="utf-8")) if config.is_file() else {}
    value = os.environ.get(name) or values.get(name)
    if not value or any(token in str(value) for token in ("<PACKAGE_ROOT>", "<", ">")):
        raise RuntimeError(f"Set {name} explicitly or use reproduction/prepare_portable.py")
    result = Path(value).expanduser()
    if not result.is_absolute():
        raise ValueError(f"{name} must be an absolute location")
    return result.resolve()

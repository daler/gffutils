import re
from importlib.metadata import PackageNotFoundError, version as distribution_version
from pathlib import Path


def _detect_version():
    """
    v0.14 migrated to pyproject.toml format, and the version is now only stored
    there. If this package is installed, resolve the installed version.
    Otherwise, inspect pyproject.toml.
    """
    try:
        return distribution_version("gffutils")
    except PackageNotFoundError:
        pyproject = Path(__file__).resolve().parent.parent / "pyproject.toml"
        try:
            contents = pyproject.read_text(encoding="utf-8")
        except OSError:
            return "0+unknown"

        # tomllib is in py3.11+ and we're supporting earlier versions, so rely
        # on regex here. Add "+unknown" to indicate possible divergence from
        # the cloned checkout.
        match = re.search(r'^version = "([^"]+)"$', contents, re.MULTILINE)
        if match:
            return match.group(1) + "+unknown"
        return "0+unknown"


version = _detect_version()

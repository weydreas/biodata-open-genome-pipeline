#!/usr/bin/env python3
"""Dante knowledge enrichment entrypoint — v9.0.0 hotfix8 / rules-v8."""
from __future__ import annotations

import sys
import knowledge_enrichment_base_v7 as _base
from knowledge_enrichment_hotfix8 import install as _install

_install(_base)

if __name__ == "__main__":
    raise SystemExit(_base.main())

# Make normal imports resolve to the patched base module object itself. This
# preserves monkey-patching/introspection behavior used by the regression suite.
sys.modules[__name__] = _base

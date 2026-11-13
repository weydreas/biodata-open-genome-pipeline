#!/usr/bin/env python3
"""Layout hardening for disease-chaptered Stage-24 PDF reports.

Hotfix12 changes only ReportLab pagination behavior. The scientific Stage-24 data,
clinical chapter taxonomy, common-major filtering, legal split and output versioning
remain unchanged.

The disease-chapter summary table can contain variants linked to very large condition
sets. ReportLab Table rows are not splittable by default, so one unusually tall cell
can exceed the page frame and raise LayoutError. This layer enables in-row splitting
for all report tables during PDF rendering, preserving the complete condition list
rather than truncating it.
"""
from __future__ import annotations

from contextlib import contextmanager

import genome_report_hotfix11 as h11

PUBLICATION_LAYOUT_HOTFIX = "dante-report-layout-v12"


@contextmanager
def splittable_report_tables():
    """Temporarily make ReportLab tables splittable inside oversized rows.

    Unit tests may substitute a fake PDF renderer and intentionally run without
    ReportLab installed. In that case this context is a no-op; production Stage 25/26
    already preflights and installs pinned ReportLab 4.2.5 before the generator runs.
    """
    try:
        from reportlab import platypus
    except ModuleNotFoundError:
        yield
        return

    original_table = platypus.Table

    def table_with_row_splitting(*args, **kwargs):
        kwargs.setdefault("splitInRow", 1)
        return original_table(*args, **kwargs)

    platypus.Table = table_with_row_splitting
    try:
        yield
    finally:
        platypus.Table = original_table


def main() -> int:
    # Preserve the existing renderer seam used by regression tests. Production
    # normally points render_report_pdf to h11.render_report_pdf; callers may
    # monkeypatch this module-level symbol before invoking main().
    current_renderer = h11.render_report_pdf
    h11.render_report_pdf = globals().get("render_report_pdf", current_renderer)
    try:
        with splittable_report_tables():
            return h11.main()
    finally:
        h11.render_report_pdf = current_renderer


for _name in dir(h11):
    if _name.startswith("__") or _name == "main":
        continue
    globals().setdefault(_name, getattr(h11, _name))


if __name__ == "__main__":
    raise SystemExit(main())

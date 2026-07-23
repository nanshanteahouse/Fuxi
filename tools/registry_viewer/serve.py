#!/usr/bin/env python3
"""Fuxi Registry Viewer launcher.

Serves viewer.html from this directory while routing the three YAML data files
(papers.yaml / datasets.yaml / links.yaml) to a user-specified data directory.

Usage:
    python tools/registry_viewer/serve.py
    python tools/registry_viewer/serve.py --data /abs/path/to/yaml
    python tools/registry_viewer/serve.py --port 8000 --open
"""

from __future__ import annotations

import argparse
import functools
import http.server
import os
import sys
import webbrowser

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
DEFAULT_DATA = os.path.join(REPO, "projects", "registry")
DATA_FILES = ("papers.yaml", "datasets.yaml", "links.yaml")


class _Handler(http.server.SimpleHTTPRequestHandler):
    """Serve viewer.html from HERE, but route the three YAML data files
    to a separate data directory."""

    def __init__(self, *args, data_dir: str, **kwargs):
        self._data_dir = data_dir
        super().__init__(*args, directory=HERE, **kwargs)

    def translate_path(self, path: str) -> str:
        name = os.path.basename(path)
        if name in DATA_FILES:
            return os.path.join(self._data_dir, name)
        return super().translate_path(path)

    def log_message(self, fmt, *args):
        pass


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--data", default=DEFAULT_DATA, help=f"YAML data directory (default: {DEFAULT_DATA})"
    )
    ap.add_argument("--port", type=int, default=8000, help="listen port (default: 8000)")
    ap.add_argument("--host", default="127.0.0.1", help="listen host (default: 127.0.0.1)")
    ap.add_argument("--open", action="store_true", help="open the viewer in the default browser")
    args = ap.parse_args()

    data_dir = os.path.abspath(args.data)
    if not os.path.isdir(data_dir):
        print(f"ERROR: data dir not found: {data_dir}", file=sys.stderr)
        sys.exit(1)

    missing = [f for f in DATA_FILES if not os.path.isfile(os.path.join(data_dir, f))]
    if missing:
        print(f"WARNING: missing in {data_dir}: {', '.join(missing)}", file=sys.stderr)

    handler_factory = functools.partial(_Handler, data_dir=data_dir)

    url = f"http://{args.host}:{args.port}/viewer.html"
    print("Fuxi Registry Viewer")
    print(f"  data : {data_dir}")
    print(f"  url  : {url}")
    print("  (Ctrl-C to stop)")

    if args.open:
        webbrowser.open(url)

    try:
        http.server.HTTPServer((args.host, args.port), handler_factory).serve_forever()
    except KeyboardInterrupt:
        print("\nstopped", file=sys.stderr)


if __name__ == "__main__":
    main()

# Fuxi Registry Viewer

Browser-based viewer for the Fuxi paper / dataset / link registry. Renders three
YAML files (`papers.yaml`, `datasets.yaml`, `links.yaml`) as tabbed, searchable
tables.

The viewer ships with the framework; your YAML data stays in your private
`projects/` tree and is read-only at runtime.

## Quick start

```bash
python tools/registry_viewer/serve.py --open
```

Starts an HTTP server on http://127.0.0.1:8000/viewer.html, pulling data from
`projects/registry/` (the default).

## Options

| Flag | Default | Purpose |
|------|---------|---------|
| `--data <dir>` | `projects/registry` | Directory containing the three YAML files |
| `--port <n>` | `8000` | Listen port |
| `--host <h>` | `127.0.0.1` | Listen host |
| `--open` | off | Open the viewer in the default browser on start |

Point at any data directory:

```bash
python tools/registry_viewer/serve.py \
  --data /mnt/e/my-fuxi-data/registry \
  --port 8765 --open
```

## Expected data layout

```
<data-dir>/
├── papers.yaml
├── datasets.yaml
└── links.yaml
```

Missing files are reported at startup but do not block it (the viewer will show
a load error for the missing tables).

## How it works

`serve.py` subclasses `http.server.SimpleHTTPRequestHandler` and overrides
`translate_path` to route the three YAML filenames to your `--data` directory,
while serving `viewer.html` and its assets from this directory. No files are
copied or symlinked; the data directory is read-only.

## Features

- Three tabs: Papers / Datasets / Links
- Global search across all tables (count summary follows the active tab)
- Per-column sort (click headers)
- Modality filter pills (Datasets tab only)
- Paper detail expansion (click the triangle)
- Light / dark mode (follows OS preference)

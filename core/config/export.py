"""Export resolved configuration — merge global + project → single YAML file.

Call ``export_resolved_config(cfg, merged_dict, source_map, global_path, project_path)``
after ``resolve_config`` to produce the self-contained ``config_resolved.yaml``
that external researchers can reproduce results from directly.
"""

import datetime
import os

import yaml


def export_resolved_config(
    cfg,
    merged_dict: dict,
    source_map: dict,
    global_path: str | None,
    project_path: str,
) -> str | None:
    """Export the resolved config to ``{results_dir}/config_resolved.yaml``.

    Returns the output path, or None if writing fails.
    """
    results_dir = cfg.results_dir
    os.makedirs(results_dir, exist_ok=True)

    # Build field-level source map
    field_sources: dict[str, str] = {}
    for path, origin in source_map.items():
        if origin == "override":
            field_sources[path] = os.path.basename(project_path)
        elif origin == "base" and global_path:
            field_sources[path] = os.path.basename(global_path)
        else:
            field_sources[path] = "schema_defaults"

    # Build meta section
    meta = {
        "_config_meta": {
            "generated_at": datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
            "tool_version": "fuxi/1.0",
            "merged_priority": "project > global > schema_default",
            "sources": [],
            "field_sources": dict(sorted(field_sources.items())),
        }
    }

    # Add source entries
    if global_path:
        meta["_config_meta"]["sources"].append(
            {
                "type": "global",
                "path": os.path.relpath(global_path),
            }
        )
    meta["_config_meta"]["sources"].append(
        {
            "type": "project",
            "path": os.path.relpath(project_path),
        }
    )

    # Build output: serialize cfg as dict, add meta
    output = cfg.model_dump(mode="json", by_alias=False)
    output = {**meta, **output}  # meta first, then all fields

    output_path = os.path.join(results_dir, "config_resolved.yaml")
    try:
        with open(output_path, "w") as f:
            yaml.safe_dump(
                output, f, sort_keys=False, allow_unicode=True, default_flow_style=False
            )
        return output_path
    except Exception:
        return None


def resolve_cli():
    """CLI entry point: ``python -m core.config resolve --project <path>``."""
    import argparse

    parser = argparse.ArgumentParser(description="Resolve Fuxi configuration")
    parser.add_argument("--project", required=True, help="Path to project config YAML")
    parser.add_argument("--output", "-o", default=None, help="Output file (default: stdout)")
    args = parser.parse_args()

    from core.utils._config import resolve_config

    cfg = resolve_config(args.project)

    output = cfg.model_dump(mode="json", by_alias=False)

    yaml_str = yaml.safe_dump(
        output, sort_keys=False, allow_unicode=True, default_flow_style=False
    )

    if args.output:
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        with open(args.output, "w") as f:
            f.write(yaml_str)
    else:
        print(yaml_str)


if __name__ == "__main__":
    resolve_cli()

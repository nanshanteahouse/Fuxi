"""CLI entry point for ``python -m core.config <command> [args]``.

Supported commands:
  resolve   — export resolved config to stdout/file
  scaffold  — render starter config YAML from specs (--list/--format/--out)
"""

import sys


def main():
    if len(sys.argv) < 2:
        print(
            "Usage: python -m core.config <command> [args]\n",
            "Commands:\n",
            "  resolve   --project <path>  Resolve and export config to stdout/file\n",
            "  scaffold [--list|--format KEY [--out PATH]]  Render starter configs from specs",
        )
        sys.exit(1)

    command = sys.argv[1]

    if command == "resolve":
        # Strip "resolve" so resolve_cli()'s argparse sees --project directly
        sys.argv.pop(1)
        from core.config.export import resolve_cli

        resolve_cli()
    elif command == "scaffold":
        sys.argv.pop(1)
        from core.config.scaffold import main as scaffold_main

        sys.exit(scaffold_main())
    else:
        print(f"Unknown command: {command}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

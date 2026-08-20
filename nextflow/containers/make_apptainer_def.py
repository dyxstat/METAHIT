#!/usr/bin/env python3
import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", required=True)
    parser.add_argument("--lock-dir", required=True)
    parser.add_argument("--pip-requirements", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    template = Path(args.template)
    lock_dir = Path(args.lock_dir)
    pip_requirements = Path(args.pip_requirements)
    out = Path(args.out)

    for name in [
        "metahict_env.explicit.txt",
        "checkm2.explicit.txt",
        "gtdbtk-2.4.0.explicit.txt",
        "genomad.explicit.txt",
        "ccfind_env.explicit.txt",
    ]:
        if not (lock_dir / name).is_file():
            raise FileNotFoundError(f"Missing Conda lock: {lock_dir / name}")
    if not pip_requirements.is_file():
        raise FileNotFoundError(f"Missing Pip requirements file: {pip_requirements}")

    out.write_text(template.read_text())


if __name__ == "__main__":
    main()

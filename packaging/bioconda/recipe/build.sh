#!/usr/bin/env bash
set -euo pipefail

install -d "${PREFIX}/share/metahict" "${PREFIX}/bin"
cp -R . "${PREFIX}/share/metahict/"
install -m 0755 metahict "${PREFIX}/bin/metahict"

cat > "${PREFIX}/bin/metahict-nextflow" <<'EOF'
#!/usr/bin/env python3
import os
from pathlib import Path
import sys

prefix_dir = Path(__file__).resolve().parents[1]
workflow = prefix_dir / "share" / "metahict" / "nextflow" / "main_dsl2.nf"
os.execvp("nextflow", ["nextflow", "run", str(workflow), *sys.argv[1:]])
EOF

chmod +x "${PREFIX}/bin/metahict-nextflow"

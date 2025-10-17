#!/usr/bin/env bash
# Fix hydrogen naming for 5MC residues in a PDB:
#   H41  -> HN41   (delete leading space if present)
#   H42  -> HN42
#   H61  -> H6     (write back as " H6 " in cols 13–16)

set -euo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: $0 <input.pdb> [output.pdb]" >&2
  exit 1
fi

in="$1"
out="${2:-${in%.pdb}_hnfix.pdb}"

awk '
function trim(s){ gsub(/^ +| +$/, "", s); return s }

{
  # Only touch 5MC residues; atom name field is columns 13–16
  res  = substr($0, 18, 3)
  atom = substr($0, 13, 4)
  at   = trim(atom)

  if (res == "5MC") {
    if (at == "H41") {
      # Write exactly "HN41" in cols 13–16 (4 chars, no spaces)
      $0 = substr($0,1,12) "HN41" substr($0,17)
    } else if (at == "H42") {
      # Write exactly "HN42"
      $0 = substr($0,1,12) "HN42" substr($0,17)
    } else if (at == "H61") {
      # Write right-justified hydrogen name " H6 " (space,H,6,space)
      $0 = substr($0,1,12) " H6 " substr($0,17)
    }
  }
  print $0
}
' "$in" > "$out"

echo "Wrote fixed file: $out"

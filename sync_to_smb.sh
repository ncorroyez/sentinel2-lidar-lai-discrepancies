#!/usr/bin/env bash
# sync_to_smb.sh — synchronise uniquement les données nécessaires à la révision
#                  vers le partage INRAE (pnas3 / mo-mtd-pulse).
#
# Lancer depuis la racine du projet (NC_Full/) :
#   bash sync_to_smb.sh            # tout
#   bash sync_to_smb.sh --outputs  # sorties de révision uniquement
#   bash sync_to_smb.sh --inputs   # données d'entrée uniquement
#
# Volume : ~3.8 Go (03_RESULTS/Metrics+PROSAIL+LiDAR, 01_DATA, PROSAIL-Opt, revision/output)

set -euo pipefail
cd "$(dirname "$0")"

# ── Détection du montage SMB ───────────────────────────────────────────────────

SMB_SERVER="pnas3.stockage.inrae.fr"
SMB_SHARE="mo-mtd-pulse"
SMB_SUBPATH="root/_PROJETS/2023_2026_These_Nathan_Corroyez/NC_Full"

UID_VAL=$(id -u)
LINUX_MOUNT="/run/user/${UID_VAL}/gvfs/smb-share:server=${SMB_SERVER},share=${SMB_SHARE}/${SMB_SUBPATH}"
MACOS_MOUNT="/Volumes/${SMB_SHARE}/${SMB_SUBPATH}"

if   [[ -d "$LINUX_MOUNT" ]]; then DEST="$LINUX_MOUNT"
elif [[ -d "$MACOS_MOUNT" ]]; then DEST="$MACOS_MOUNT"
else
  echo "ERREUR : NC_Full/ non détecté sur le partage SMB."
  echo "  Attendu : $LINUX_MOUNT"
  echo "Monter le partage d'abord (Nautilus → Autres emplacements)."
  exit 1
fi

echo "Destination : $DEST"
echo ""

# ── Mode ───────────────────────────────────────────────────────────────────────

MODE="${1:-}"
do_inputs=true
do_outputs=true
[[ "$MODE" == "--outputs" ]] && do_inputs=false
[[ "$MODE" == "--inputs"  ]] && do_outputs=false

# ── Fonction rsync ─────────────────────────────────────────────────────────────

sync_dir() {
  local src="$1" dest_rel="$2"
  if [[ ! -d "$src" ]]; then
    echo "  ⚠ absent en local, ignoré : $src"
    return
  fi
  local dest_abs="$DEST/$dest_rel"
  mkdir -p "$dest_abs"
  echo "▶ $(date '+%H:%M:%S')  $src"
  rsync --recursive --links --times --update --inplace \
    --no-perms --no-owner --no-group \
    --info=progress2 --human-readable \
    "$src/" "$dest_abs/"
  echo ""
}

SITES=(Aigoual Blois Mormal)

# ── Données d'entrée ───────────────────────────────────────────────────────────
# Source locale : revision/data/{results,raw,prosail_lidar}/
# Destination   : NC_Full/{03_RESULTS,01_DATA,PROSAIL-Optimization/01_DATA}/

if $do_inputs; then
  echo "══ Données d'entrée ════════════════════════════════════════════════════"

  for site in "${SITES[@]}"; do
    echo "── $site ──"

    # Métriques ALS (LAI, fCover, max height, LAD stack…)
    sync_dir "revision/data/results/$site/Metrics/Deciduous_Only" \
             "03_RESULTS/$site/Metrics/Deciduous_Only"

    # Optimisation PROSAIL (réflectances S2, géométries d'acquisition…)
    sync_dir "revision/data/results/$site/PROSAIL_Optimization" \
             "03_RESULTS/$site/PROSAIL_Optimization"

    # CHM 1 m (normalisation hauteurs)
    sync_dir "revision/data/results/$site/LiDAR/chm/res_1_m" \
             "03_RESULTS/$site/LiDAR/chm/res_1_m"

    # DSM 1 m (hétérogénéité)
    sync_dir "revision/data/results/$site/LiDAR/dsm" \
             "03_RESULTS/$site/LiDAR/dsm"

    # Fichiers géographiques (emprise, vecteurs)
    sync_dir "revision/data/raw/$site/Geo_Files" \
             "01_DATA/$site/Geo_Files"

    # Rasters LiDAR bruts utilisés (lskew, mean)
    sync_dir "revision/data/raw/$site/LiDAR" \
             "01_DATA/$site/LiDAR"

    # Nuages de points LiDAR pour PROSAIL (LAD profiles)
    sync_dir "revision/data/prosail_lidar/$site/LiDAR" \
             "PROSAIL-Optimization/01_DATA/$site/LiDAR"
  done

  # Fonction utilitaire angles S2 (seul fichier code externe nécessaire)
  sync_dir "PROSAIL-Optimization/02_CODES/libraries" \
           "PROSAIL-Optimization/02_CODES/libraries"
fi

# ── Sorties de révision ────────────────────────────────────────────────────────

if $do_outputs; then
  echo "══ Sorties de révision ══════════════════════════════════════════════════"
  sync_dir "revision/output" "revision/output"
fi

# ── Résumé ────────────────────────────────────────────────────────────────────

echo "Synchronisation terminée — $(date '+%H:%M:%S')"
echo "Destination : $DEST"

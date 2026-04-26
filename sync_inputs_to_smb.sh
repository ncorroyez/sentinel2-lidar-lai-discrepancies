#!/usr/bin/env bash
# sync_inputs_to_smb.sh — copie les données d'entrée vers le partage INRAE.
#
# Prérequis : partage monté via Nautilus (GNOME) ou Finder (macOS).
# Lancer depuis la racine du projet (NC_Full/) :
#   bash sync_inputs_to_smb.sh
#
# Utilise rsync --update : seuls les fichiers absents ou plus récents en local
# sont transférés. Le script est donc relançable sans tout recopier.

set -euo pipefail
cd "$(dirname "$0")"

# ── Détection du point de montage SMB ─────────────────────────────────────────

SMB_SERVER="pnas3.stockage.inrae.fr"
SMB_SHARE="mo-mtd-pulse"
SMB_SUBPATH="root/_PROJETS/2023_2026_These_Nathan_Corroyez/NC_Full"

# Linux (GNOME/GVfs)
UID_VAL=$(id -u)
LINUX_MOUNT="/run/user/${UID_VAL}/gvfs/smb-share:server=${SMB_SERVER},share=${SMB_SHARE}/${SMB_SUBPATH}"

# macOS
MACOS_MOUNT="/Volumes/${SMB_SHARE}/${SMB_SUBPATH}"

if   [[ -d "$LINUX_MOUNT" ]]; then
  DEST="$LINUX_MOUNT"
elif [[ -d "$MACOS_MOUNT" ]]; then
  DEST="$MACOS_MOUNT"
else
  echo "ERREUR : partage SMB non détecté."
  echo "  Attendu (Linux) : $LINUX_MOUNT"
  echo "  Attendu (macOS) : $MACOS_MOUNT"
  echo "Monter le partage d'abord, puis relancer."
  exit 1
fi

echo "Destination : $DEST"
echo ""

# ── Répertoires à synchroniser ─────────────────────────────────────────────────
# Lecture seule en local — ne jamais modifier ces dossiers.

SOURCES=(
  "01_DATA"
  "03_RESULTS"
  "PROSAIL-Optimization/01_DATA"
  "PROSAIL-Optimization/02_CODES"
)

# ── Transfert ─────────────────────────────────────────────────────────────────

t_total=$SECONDS

for SRC in "${SOURCES[@]}"; do
  echo "▶ $(date '+%H:%M:%S') — $SRC"

  # Créer le répertoire parent sur la destination si besoin
  DEST_PARENT="$DEST/$(dirname "$SRC")"
  mkdir -p "$DEST_PARENT"

  rsync \
    --archive \
    --update \
    --info=progress2 \
    --human-readable \
    "$SRC/" \
    "$DEST/$SRC/"

  echo "  ✓ $SRC synchronisé"
  echo ""
done

echo "Synchronisation terminée — $(date '+%H:%M:%S') — total $(( SECONDS - t_total ))s"
echo "Destination : $DEST"

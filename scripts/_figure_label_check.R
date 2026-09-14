# ==============================================================================
# _figure_label_check.R — empêcher qu'une étiquette de figure contredise le texte
# du chapitre.
#
# Motivation. La valeur brute du R² LiDAR en canopée dense vaut exactement 0.595.
# `sprintf("%.2f", 0.595)` rend "0.59", parce que 0.595 vaut 0.59499999... en
# double ; `round(0.595, 2)` rend 0.60. Le script de figure utilisait le premier,
# le texte du chapitre le second, et la figure a affiché pendant des mois un
# chiffre que le manuscrit contredisait de quatre endroits différents.
#
# Ce fichier fournit :
#   fmt2(x)                   formateur canonique 2 décimales, arrondi au demi
#                             supérieur — à utiliser dans TOUTES les figures
#   check_labels(x, md, what) échoue si le chapitre énonce l'autre arrondi
#
# Usage :
#   source("scripts/_figure_label_check.R")
#   Jm[, lab := fmt2(R2)]
#   check_labels(Jm$R2, "manuscripts/ch3/Chapter3_article_standalone_EN.md",
#                what = "Fig J — R² stratifiés")
# ==============================================================================

#' Arrondi au demi supérieur, deux décimales
#'
#' @param x numérique
#' @return chaîne de caractères, toujours deux décimales
fmt2 <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.2f", floor(x * 100 + 0.5) / 100))
}

#' Vérifier que les étiquettes d'une figure ne contredisent pas le chapitre
#'
#' Pour chaque valeur brute, on calcule l'étiquette retenue (fmt2, demi
#' supérieur) et l'étiquette concurrente (sprintf, arrondi IEEE). Quand les deux
#' diffèrent — c'est-à-dire quand la valeur tombe exactement sur un demi — et que
#' le texte du chapitre contient la concurrente sans contenir la retenue, la
#' figure est en train d'imprimer un chiffre que le manuscrit dément : on
#' s'arrête.
#'
#' L'absence pure et simple d'une étiquette dans le texte n'est PAS une erreur :
#' une figure montre légitimement des valeurs que la prose ne cite pas. Elle est
#' seulement rapportée, pour information.
#'
#' @param values numérique, les valeurs brutes portées par les étiquettes
#' @param md chemin du markdown du chapitre
#' @param what libellé court identifiant la figure, pour les messages
#' @param quiet logique, TRUE pour taire le récapitulatif de couverture
#' @return invisible(TRUE) si aucun conflit
check_labels <- function(values, md, what = "figure", quiet = FALSE) {
  if (!file.exists(md)) {
    warning(sprintf("[%s] chapitre introuvable (%s) — contrôle ignoré", what, md))
    return(invisible(FALSE))
  }
  txt <- paste(readLines(md, warn = FALSE, encoding = "UTF-8"), collapse = " ")
  v <- unique(values[is.finite(values)])

  kept <- fmt2(v)
  rival <- sprintf("%.2f", v)
  amb <- which(kept != rival)          # valeurs pile sur un demi

  conflicts <- character(0)
  for (i in amb) {
    has_kept  <- grepl(kept[i],  txt, fixed = TRUE)
    has_rival <- grepl(rival[i], txt, fixed = TRUE)
    if (has_rival && !has_kept)
      conflicts <- c(conflicts, sprintf(
        "  valeur %.4f : la figure afficherait %s, le chapitre écrit %s",
        v[i], kept[i], rival[i]))
  }
  if (length(conflicts))
    stop(sprintf("[%s] étiquette(s) en contradiction avec %s :\n%s\n",
                 what, basename(md), paste(conflicts, collapse = "\n")),
         "  -> harmoniser l'arrondi, ou corriger le chapitre.", call. = FALSE)

  if (!quiet) {
    found <- vapply(kept, function(s) grepl(s, txt, fixed = TRUE), logical(1))
    message(sprintf("[%s] %d/%d étiquettes retrouvées dans %s%s",
                    what, sum(found), length(found), basename(md),
                    if (length(amb)) sprintf(" ; %d valeur(s) sur un demi vérifiée(s)",
                                             length(amb)) else ""))
  }
  invisible(TRUE)
}

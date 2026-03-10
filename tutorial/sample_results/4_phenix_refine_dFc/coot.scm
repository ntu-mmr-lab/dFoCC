;;;; Coot script file in scheme
(set-recentre-on-read-pdb 0)
(set-nomenclature-errors-on-read "ignore")

;; load the dark and the calculated light structure
(read-pdb "./input/8Z1J-without-H.pdb")
(read-pdb "4_phenix_refine_dFc/DEMO_cycle000_initial_struct.pdb")

;; read the dFo region-extracted sigma-cutoff map
;; and show almost all signal on the map
;; (handle-read-ccp4-map "4_phenix_refine_dFc/dFo.DEMO_cycle000_initial_struct.masked_sigma3.0.ccp4" 1)
;; (set-last-map-contour-level-by-sigma 0.01)

;; read the entire dFo struct factor
;; and apply suitable sigma cutoff on it
(make-and-draw-map "./input/8Z3X-dFo.mtz" "FoFo" "PHFc" "" 0 1)
(set-last-map-contour-level-by-sigma 3.00)
(set-last-map-colour 0.25 1.0 0.25)
; (export-map-fragment 2 (rotation-centre-position 0) (rotation-centre-position 1) (rotation-centre-position 2) 20.0 "4_phenix_refine_dFc/dFo_for_pymol.map")

;; read the dFc region-extracted sigma-cutoff map
;; show almost all signal on the map, and change its colour into cyan-yellow pair
;; (handle-read-ccp4-map "4_phenix_refine_dFc/dFc.DEMO_cycle000_initial_struct.masked_sigma3.0.ccp4" 1)
;; (set-last-map-contour-level-by-sigma 0.01)

;; read the entire dFo struct factor
;; apply suitable sigma cutoff on it
;; and change its colour into cyan-yellow pair
(make-and-draw-map "4_phenix_refine_dFc/dFc.DEMO_cycle000_initial_struct.mtz" "dFc" "dPHIc" "" 0 1)
(set-last-map-contour-level-by-sigma 3.00)
(set-last-map-colour 0.0 1.0 1.0)
; (export-map-fragment 3 (rotation-centre-position 0) (rotation-centre-position 1) (rotation-centre-position 2) 20.0 "4_phenix_refine_dFc/dFc_for_pymol.map")

;; read the extrapolated struct factor
;; (auto-read-make-and-draw-maps "./input/nestimation_11.0.mtz")
(make-and-draw-map "4_phenix_refine/refined__001.mtz" "2FOFCWT_no_fill" "PH2FOFCWT_no_fill" "" 0 0)

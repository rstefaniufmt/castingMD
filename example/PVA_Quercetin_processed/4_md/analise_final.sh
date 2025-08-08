#!/bin/bash

TPR="step0.tpr"
NDX="index.ndx"
OUTDIR="analises_filme"
mkdir -p "$OUTDIR"

# Loop de 0 a 19 com passo 5
for STEP in 0 4 9 1 19 24; do
    echo "=== Análises para step${STEP}.xtc ==="
    TRAJ="step${STEP}.xtc"
    TPR="step${STEP}.tpr"
    PREFIX="step${STEP}" 
    STEP_OUT="${OUTDIR}"
    mkdir -p "$STEP_OUT"

    if [ "$STEP" -eq 0 ]; then
        TRAJ="npt.xtc"
        TPR="npt.tpr"
    fi    

   echo "--- Perfil de densidade separado (espessura do filme) ---"
 
  gmx density -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/density_PVA_${STEP}.xvg -d Z -center -sl 200 -symm -dens number << EOF
PVA
PVA
EOF

    echo "--- Perfil de densidade (espessura do filme) ---"
    gmx density -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/densidade_z_${STEP}.xvg -d Z -sl 200 <<EOF
PVA
EOF

    echo "--- RDF PVA–PEG ---"
    gmx rdf -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/rdf_PVA${STEP}.xvg -ref PVA -sel PVA

    echo "--- Densidade do filme (massa) ---"
  gmx density -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/densidade_filme_${STEP}.xvg -dens mass <<EOF
PVA
EOF

    echo "--- RMSD dos polímeros ---"
    gmx rms -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/rmsd_polimeros_${STEP}.xvg <<EOF
PVA
PVA
EOF

    echo "--- RMSF dos polímeros ---"
    gmx rmsf -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/rmsf_polimeros_${STEP}.xvg <<EOF
PVA
EOF

    echo "--- RG dos polímeros ---"

    echo "--- RG do PVA ---"
    gmx gyrate -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/rg_PVA_${STEP}.xvg <<EOF
PVA
EOF

    echo "--- Fim das análises do step${STEP} ---"
    echo
done

echo "==== Analises para Filme apos Resfriamento ====="

   mkdir -p "${OUTDIR}/final"
   
   TRAJ="filme_final.xtc"
   TPR="filme_final.tpr"
   STEP_OUT="${OUTDIR}/final"


  gmx density -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/density_PVA_final.xvg -d Z -center -sl 200 -symm -dens number << EOF
PVA
PVA
EOF


    echo "--- RDF PVA ---"
    gmx rdf -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/rdf_PVA_FINAL.xvg -ref PVA -sel PVA

    echo "--- Densidade do filme (massa) ---"
  gmx density -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/densidade_filme_FINAL.xvg -dens mass <<EOF
PVA
EOF

    echo "--- RMSD dos polímeros ---"
    gmx rms -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/rmsd_polimeros_FINAL.xvg <<EOF
PVA
PVA
EOF

    echo "--- RMSF dos polímeros ---"
    gmx rmsf -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/rmsf_polimeros_FINAL.xvg <<EOF
PVA
EOF

    echo "--- RG do PVA ---"
    gmx gyrate -f $TRAJ -s $TPR -n $NDX -o ${STEP_OUT}/rg_PVA_FINAL.xvg <<EOF
PVA
EOF

echo "=== Todas as análises finalizadas ==="



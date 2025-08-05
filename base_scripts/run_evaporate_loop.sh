#!/usr/bin/env bash

# ------------ CONFIGURATION ------------
source evaporate.conf

#if not exists configuration, set default values

SEG_NS="${SEG_NS:-1}"                      # Duration of each cycle (ns)
CYCLES="${CYCLES:-20}"                     # Total number of cycles
MDP="${MDP:-md.mdp}"                    # Production .mdp file
TOP="${TOP:-film.top}"                  # Initial topology (adjusted after solvatation)
GRO="${GRO:-npt.gro}"                  # Final structure from NPT and solvatation
TPR="${TPR:-md.tpr}" 
CPT="${CPT:-npt.cpt}"                   # Temporary CPT file
PYTHON=python3                # Python with MDAnalysis installed
LOG_FILE="${LOG_FILE:-evaporate.log}"

plane=1.00 #ONLY TO START THE VARIABLE

# ------------ MAIN LOOP ------------
for ((i=0; i<CYCLES; i++)); do
    echo "==== Starting cycle $i ====" | tee -a "$LOG_FILE"
    
    tag=$(printf "%d" ${i})
    
    if [[ -f "step${tag}_last.gro" ]]; then
        echo "File step${tag}_last.gro already exists. Skipping cycle $i." | tee -a "$LOG_FILE"
        continue
    fi
    
    # If not the first cycle, update input files
    if [[ $i -ne 0 ]]; then
        prev_tag=$(printf "%d" $((i-1)))
        CPT="step${prev_tag}_npt.cpt"
        GRO="step${prev_tag}_npt.gro"
        TOP="step${prev_tag}_post.top"
        
        # --- new NVT and MDP
        if [[ ! -f $GRO ]]; then
           GRO=step${prev_tag}_post.gro
           TOP=step${prev_tag}_post.top
           
           gmx grompp -f nvt_eq.mdp -c ${GRO} -r ${GRO} -p ${TOP} -o step${prev_tag}_nvt.tpr -maxwarn 10
           gmx mdrun -deffnm step${prev_tag}_nvt -v || { echo "Error in NVT during cycle $i"; exit 1; }

           gmx grompp -f npt_eq.mdp -c step${prev_tag}_nvt.gro -t step${prev_tag}_nvt.cpt -r step${prev_tag}_nvt.gro -p ${TOP} -o step${prev_tag}_npt.tpr -maxwarn 10
           gmx mdrun -deffnm step${prev_tag}_npt -v || { echo "Error in NPT during cycle $i"; exit 1; }
           
           CPT="step${prev_tag}_npt.cpt"
           GRO="step${prev_tag}_npt.gro"
           TOP="step${prev_tag}_post.top"
        fi   
                 
        # 3. Determine Z-height of the system
        echo ">>> Calculating new box height..." | tee -a "$LOG_FILE"
        read boxX boxY boxZ <<< $(awk 'END {print $1, $2, $3}' "$GRO")
        zmax=$(tail -n 1 "$GRO" | awk '{print $3}')

        # 4. Safety check: ensure Z >= 2.89 nm (for default rlist of 1.2 nm)
        rlist_real=$(grep "rlist from" step${prev_tag}_npt.log | tail -1 | awk '{print $11}')
        buffer=0.3
        
        if [ -z "$rlist_real" ]; then
          echo "[WARNING] Could not find rlist adjustment in log. Using TPR default value." | tee -a "$LOG_FILE"
          rlist_real=$(gmx dump -s step${prev_tag}_npt.tpr | grep "rlist" | head -1 | awk '{print $3}')
        fi
        
        minZ=$(echo "$rlist_real * 2 + $buffer" | bc -l)
        
        zsafe=$(echo "$zmax < $minZ" | bc -l)
        if [ "$zsafe" -eq 1 ]; then
            echo ">>> Calculated Z-height ($zmax nm) is too small. Correcting to ${minZ} nm" | tee -a "$LOG_FILE"
            zmax=$minZ
            cp ${GRO} old_${GRO}
            
            # 5. Redefine box with new height
            echo ">>> Applying new Z box height = ${zmax} nm" | tee -a "$LOG_FILE"
            gmx editconf -f step${prev_tag}_npt.gro -o newbox.gro -c -box ${boxX} ${boxY} ${minZ}
            # 6. Prepare for next cycle
            cp newbox.gro box_fixed_${prev_tag}.gro
            cp newbox.gro ${GRO}
            rm -rf newbox.gro
            
            echo "[INFO] Checking consistency after cycle $i" | tee -a "$LOG_FILE"
            n_atoms_gro=$(tail -n +3 ${GRO} | wc -l)
            n_sol_top=$(grep "^SOL" ${TOP} | awk '{print $2}')
            echo "  Atoms in .gro: $n_atoms_gro" | tee -a "$LOG_FILE"
            echo "  SOL molecules in .top: $n_sol_top" | tee -a "$LOG_FILE"
        fi 
    fi

    # ------ Gradually adjust plane_frac ------
    if   [ $i -le 5 ]; then plane=0.95
    elif [ $i -le 10 ]; then plane=0.85
    elif [ $i -le 12 ]; then plane=0.80
    elif [ $i -le 14 ]; then plane=0.75
    elif [ $i -le 16 ]; then plane=0.70
    elif [ $i -le 18 ]; then plane=0.60
    else                     plane=0.50
    fi

    # 1. Prepare TPR
    if [[ $i -eq 0 ]]; then
        gmx grompp -f ${MDP} -c ${GRO} -t ${CPT} -p ${TOP} -o step${tag}.tpr -maxwarn 10
    else
        gmx grompp -f ${MDP} -c ${GRO} -t ${CPT} -p ${TOP} -o step${tag}.tpr -maxwarn 10
    fi
   
    # 2. ------ Run MD ------
    gmx mdrun -deffnm step${tag} -v || { echo "Error in mdrun during cycle $i" | tee -a "$LOG_FILE"; exit 1; }

    # 3. Extract the last frame
    echo 0 | gmx trjconv -f step${tag}.trr -s step${tag}.tpr -dump $((SEG_NS*1000)) -o step${tag}_last.gro -pbc mol -ur compact
    
    rm -rf step${tag}.trr 
    # ------ Evaporation ------
    echo "Removing solvent (plane_frac=$plane)..." | tee -a "$LOG_FILE"
    python3 evaporate.py step${tag}_last.gro ${TOP} step${tag}_post.gro step${tag}_post.top ${plane} || { echo "Error during evaporation in cycle $i" | tee -a "$LOG_FILE"; exit 1; }

    # ------ Update for next cycle ------
    GRO=step${tag}_post.gro
    TOP=step${tag}_post.top
    
       
    # 3. Determine Z-height of the system
    echo ">>> Calculating new box height..." | tee -a "$LOG_FILE"
    
    read boxX boxY boxZ <<< $(awk 'END {print $1, $2, $3}' "$GRO")
    zmax=$(tail -n 1 "$GRO" | awk '{print $3}')

    # 4. Safety check: ensure Z >= 2.75 nm (for default rlist of 1.2 nm)
    
    rlist_real=$(grep "rlist from" step${prev_tag}_npt.log | tail -1 | awk '{print $11}')
    buffer=0.3
        
   if [ -z "$rlist_real" ]; then
       echo "[WARNING] Could not find rlist adjustment in log. Using TPR default value." | tee -a "$LOG_FILE"
       rlist_real=$(gmx dump -s step${prev_tag}_npt.tpr | grep "rlist" | head -1 | awk '{print $3}')
   fi
        
    minZ=$(echo "$rlist_real * 2 + $buffer" | bc -l)
    
    
    zsafe=$(echo "$zmax < $minZ" | bc -l)
    if [ "$zsafe" -eq 1 ]; then
        echo ">>> Calculated Z-height ($zmax nm) is too small. Correcting to ${minZ} nm" | tee -a "$LOG_FILE"
        zmax=$minZ
        cp ${GRO} old_${GRO}
        
        # 5. Redefine box with new height
        echo ">>> Applying new Z box height = ${zmax} nm" | tee -a "$LOG_FILE"
        gmx editconf -f step${tag}_post.gro -o newbox.gro -c -box ${boxX} ${boxY} ${minZ}
        
        # 6. Prepare for next cycle
        cp newbox.gro box_fixed_${tag}.gro
        cp newbox.gro ${GRO}
        rm -rf newbox.gro
        
        echo "[INFO] Checking consistency after cycle $i" | tee -a "$LOG_FILE"
        n_atoms_gro=$(tail -n +3 ${GRO} | wc -l)
        n_sol_top=$(grep "^SOL" ${TOP} | awk '{print $2}')
        echo "  Atoms in .gro: $n_atoms_gro" | tee -a "$LOG_FILE"
        echo "  SOL molecules in .top: $n_sol_top" | tee -a "$LOG_FILE"
    fi 
    
    # --- Minimization
    gmx grompp -f em.mdp -c ${GRO} -p ${TOP} -o step${tag}_min.tpr -maxwarn 10
    gmx mdrun -deffnm step${tag}_min -v || { echo "Error in EM during cycle $i" | tee -a "$LOG_FILE"; exit 1; }
    
    # --- new NVT and MDP
    gmx grompp -f nvt_eq.mdp -c step${tag}_min.gro -r step${tag}_min.gro -p ${TOP} -o step${tag}_nvt.tpr -maxwarn 10
    gmx mdrun -deffnm step${tag}_nvt -v || { echo "Error in NVT during cycle $i" | tee -a "$LOG_FILE"; exit 1; }

    gmx grompp -f npt_eq.mdp -c step${tag}_nvt.gro -t step${tag}_nvt.cpt -r step${tag}_nvt.gro -p ${TOP} -o step${tag}_npt.tpr -maxwarn 10
    gmx mdrun -deffnm step${tag}_npt -v || { echo "Error in NPT during cycle $i" | tee -a "$LOG_FILE"; exit 1; }
    
    CPT=step${tag}_npt.cpt
    GRO=step${tag}_npt.gro
   
done

rm -rf *#*

# Final annealing
echo "==== Simulation completed up to cycle $CYCLES, cooling system... ====" | tee -a "$LOG_FILE"

CPT=step19_npt.cpt
GRO=step19_npt.gro
TOP=step19_post.top

# 3. Determine Z-height of the system
echo ">>> Calculating new box height..."
zmax=$(tail -n 1 "$GRO" | awk '{print $3}')

# 4. Safety check: ensure Z >= 2.75 nm (for default rlist of 1.2 nm)
minZ=2.75
zsafe=$(echo "$zmax < $minZ" | bc -l)
if [ "$zsafe" -eq 1 ]; then
    echo ">>> Calculated Z-height ($zmax nm) is too small. Correcting to ${minZ} nm"
    read boxX boxY boxZ <<< $(awk 'END {print $1, $2, $3}' "$GRO")
    zmax=$minZ
    cp ${GRO} old_${GRO}
    
    # 5. Redefine box with new height
    echo ">>> Applying new Z box height = ${zmax} nm"
    gmx editconf -f step19_post.gro -o newbox.gro -c -box ${boxX} ${boxY} ${minZ}
    
    # 6. Prepare for next cycle
    cp newbox.gro box_fixed_${tag}.gro
    cp newbox.gro ${GRO}
    rm -rf newbox.gro
    
    echo "[INFO] Checking consistency after cycle $i"
    n_atoms_gro=$(tail -n +3 ${GRO} | wc -l)
    n_sol_top=$(grep "^SOL" ${TOP} | awk '{print $2}')
    echo "  Atoms in .gro: $n_atoms_gro"
    echo "  SOL molecules in .top: $n_sol_top"
fi 

gmx grompp -f md_anneling.mdp -c step19_npt.gro -t step19_npt.cpt -p step19_post.top -o filme_final.tpr -maxwarn 10 || { echo "Error during annealing"; exit 1; }
gmx mdrun -deffnm filme_final -v || { echo "Error during annealing" | tee -a "$LOG_FILE"; exit 1; }

echo "==== Simulation finished successfully!! ====" | tee -a "$LOG_FILE"


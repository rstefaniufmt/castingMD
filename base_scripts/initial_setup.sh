#!/bin/bash

# ==============================================================================
# WARNING: This script generates the initial system for evaporation (npt.gro)
# as mentioned extensively in the documentation.
#
# NOTE ON SOLVENTS:
# This script is designed for SINGLE SOLVENT systems (defaulting to SPC water).
# For binary or ternary solvent systems, you must manually edit this script
# and adapt Step 3. Using mixed solvents requires a pre-equilibrated solvent
# box, which demands intermediate to advanced knowledge of GROMACS.
# ==============================================================================

# ==========================================
# Reading setup.conf
# ==========================================
CONF_FILE="setup.conf"

if [ ! -f "$CONF_FILE" ]; then
    echo "Error: File $CONF_FILE not found in the current directory!"
    exit 1
fi

keys=()
declare -A molecules
box_dims=""

echo "=> Reading $CONF_FILE..."
while IFS='=' read -r key value; do
    # Remove leading and trailing whitespaces
    key=$(echo "$key" | xargs)
    value=$(echo "$value" | xargs)
    
    # Ignore empty lines or comments (starting with # or ;)
    if [ -z "$key" ] || [[ "$key" == \#* ]] || [[ "$key" == \;* ]]; then
        continue
    fi
    
    if [ "$key" == "BOX" ]; then
        box_dims="$value"
    else
        keys+=("$key")
        molecules["$key"]="$value"
    fi
done < "$CONF_FILE"

if [ -z "$box_dims" ]; then
    echo "Error: BOX directive not defined in $CONF_FILE!"
    exit 1
fi

# ==========================================
# Step 1 and 2 - Create box, insert molecules and generate topology
# ==========================================
first=true
output_gro="filme_inicial.gro"

echo "=> Generating film.top file..."
cat <<EOF > film.top
; Include forcefield parameters
#include "gromos53a6.ff/forcefield.itp"

; Include chain topologies
EOF

for mol in "${keys[@]}"; do
    count="${molecules[$mol]}"
    
    # Search for the file with the exact name or all lowercase
    if [ -f "${mol}.gro" ]; then
        mol_gro="${mol}.gro"
        mol_itp="${mol}.itp"
    elif [ -f "${mol,,}.gro" ]; then
        mol_gro="${mol,,}.gro"
        mol_itp="${mol,,}.itp"
    else
        echo "Error: Corresponding .gro file for molecule $mol not found!"
        exit 1
    fi
    
    # Add the molecule's .itp to film.top
    if [ -f "$mol_itp" ]; then
        echo "#include \"$mol_itp\"" >> film.top
    else
        echo "; Warning: File $mol_itp not found. Please check if the name matches." >> film.top
    fi
    
    # Insert molecules into the box
    if [ "$first" = true ]; then
        echo "=> Creating box with $box_dims nm and inserting the first molecule of $mol..."
        gmx editconf -f "$mol_gro" -o box.gro -c -box $box_dims
        
        insert_count=$((count - 1))
        if [ "$insert_count" -gt 0 ]; then
            echo "=> Inserting $insert_count additional molecules of $mol (Total = $count)..."
            gmx insert-molecules -f box.gro -o "$output_gro" -ci "$mol_gro" -nmol "$insert_count"
        else
            cp box.gro "$output_gro"
        fi
        first=false
    else
        echo "=> Inserting $count molecules of $mol..."
        gmx insert-molecules -f "$output_gro" -o temp_filme.gro -ci "$mol_gro" -nmol "$count"
        mv temp_filme.gro "$output_gro"
    fi
done

# Add solvent and ions to the topology
cat <<EOF >> film.top

; Include water topology
#include "gromos53a6.ff/spc.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include topology for ions
#include "gromos53a6.ff/ions.itp"

[ system ]
; Name
Polymeric Film System

[ molecules ]
; Compound        #mols
EOF

# Fill in the count of each molecule in the topology
for mol in "${keys[@]}"; do
    count="${molecules[$mol]}"
    printf "%-15s %d\n" "$mol" "$count" >> film.top
done

echo "" >> film.top

# ==========================================
# Step 3 - Solvate the system (SPC model)
# ==========================================
echo "=> Solvating with SPC water..."
gmx solvate -cp "$output_gro" -cs spc216.gro -o solucao_filmogenica.gro -p film.top

# ==========================================
# Step 4 - Energy Minimization
# ==========================================
echo "=> Executing Energy Minimization..."
gmx grompp -f ./mdp/em.mdp -c solucao_filmogenica.gro -p film.top -o em.tpr -maxwarn 2
gmx mdrun -v -deffnm em

# ==========================================
# Step 8 - Equilibration Simulations
# ==========================================

# 8.1 NVT Equilibration (Constant Volume and Temperature)
echo "=> Starting NVT Equilibration..."
gmx grompp -f ./mdp/nvt.mdp -c em.gro -r em.gro -p film.top -o nvt.tpr -maxwarn 3
gmx mdrun -v -deffnm nvt

# 8.2 NPT Equilibration (Constant Pressure and Temperature)
echo "=> Starting NPT Equilibration..."
# Note: Assumed npt.mdp is in the ./mdp/ folder like the other files.
gmx grompp -f ./mdp/npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p film.top -o npt.tpr -maxwarn 3
gmx mdrun -v -deffnm npt

echo "=> Preparation completed successfully!"
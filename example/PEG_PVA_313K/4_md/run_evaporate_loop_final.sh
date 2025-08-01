#!/usr/bin/env bash

# ------------ CONFIGURAÇÃO ------------
SEG_NS=1                      # Tempo de cada ciclo (ns)
CYCLES=20                     # Número total de ciclos
MDP=md.mdp                    # Arquivo .mdp de produção
TOP=film.top           # Topologia inicial (ajustada após o passo 14)
GRO=npt.gro           # Estrutura final do passo 14
TPR=md.tpr 
CPT=npt.cpt                   # Arquivo TPR temporário
PYTHON=python3               # Python com MDAnalysis instalado
plane=1.00


# ------------ LOOP PRINCIPAL ------------
for ((i=0; i<CYCLES; i++)); do
    echo "==== Iniciando ciclo $i ===="
    
    tag=$(printf "%d" ${i})
    
    if [[ -f "step${tag}_last.gro" ]]; then
        echo "Arquivo step${tag}_last.gro já existe. Pulando ciclo $i."
        continue
    fi
    
    # Se não for o primeiro ciclo, atualize os arquivos de entrada
    if [[ $i -ne 0 ]]; then
        prev_tag=$(printf "%d" $((i-1)))
        CPT="step${prev_tag}_npt.cpt"
        GRO="step${prev_tag}_npt.gro"
        TOP="step${prev_tag}_post.top"
        
        # --- novo nvt e mdp
        if [[ ! -f $GRO ]]; then
           GRO=step${prev_tag}_post.gro
           TOP=step${prev_tag}_post.top
           
             # 3. Determina altura Z do sistema
           echo ">>> Calculando nova altura da caixa..."
           zmax=$(tail -n 1 "$GRO" | awk '{print $3}')

           # 4. Proteção: garantir que Z >= 2.89 nm (para rlist padrão de 1.2 nm)
          minZ=2.75
          zsafe=$(echo "$zmax < $minZ" | bc -l)
        if [ "$zsafe" -eq 1 ]; then
            echo ">>> Altura Z calculada ($zmax nm) é muito pequena. Corrigindo para ${minZ} nm"
            zmax=$minZ
            cp ${GRO} old_${GRO}
            
            # 5. Redefine caixa com nova altura
            echo ">>> Aplicando nova altura da caixa Z = ${zmax} nm"
            gmx editconf -f step${prev_tag}_npt.gro -o newbox.gro -c -box 6.0 6.0 ${minZ}
            # 6. Prepara para próximo ciclo
            cp newbox.gro box_fixed_${prev_tag}.gro
            cp newbox.gro ${GRO}
            rm -rf newbox.gro
            
            echo "[INFO] Verificando consistência após ciclo $i"
            n_atoms_gro=$(tail -n +3 ${GRO} | wc -l)
            n_sol_top=$(grep "^SOL" ${TOP} | awk '{print $2}')
            echo "  Átomos no .gro: $n_atoms_gro"
            echo "  Moléculas de SOL no .top: $n_sol_top"           
            python3 redist_solv.py -i box_fixed_${tag}.gro -o ${GRO} -b 6.0 6.0 ${minZ} 
        fi 
           
           gmx grompp -f nvt_eq.mdp -c ${GRO} -r ${GRO} -p ${TOP} -o step${prev_tag}_nvt.tpr -maxwarn 10
           gmx mdrun -deffnm step${prev_tag}_nvt -v || { echo "Erro no NVT no ciclo $i"; exit 1; }

           gmx grompp -f npt_eq.mdp -c step${prev_tag}_nvt.gro -t step${prev_tag}_nvt.cpt -r step${prev_tag}_nvt.gro -p ${TOP} -o step${prev_tag}_npt.tpr -maxwarn 10
           gmx mdrun -deffnm step${prev_tag}_npt -v || { echo "Erro no NPT no ciclo $i"; exit 1; }
           
           CPT="step${prev_tag}_npt.cpt"
           GRO="step${prev_tag}_npt.gro"
           TOP="step${prev_tag}_post.top"
        fi   
                 
              
       
    fi


    # ------ Ajusta o plane_frac progressivamente ------
    if   [ $i -le 5 ]; then plane=0.95
    elif [ $i -le 10 ]; then plane=0.85
    elif [ $i -le 12 ]; then plane=0.80
    elif [ $i -le 14 ]; then plane=0.75
    elif [ $i -le 16 ]; then plane=0.70
    elif [ $i -le 18 ]; then plane=0.60
    else                     plane=0.50
    fi

# 1. Preparar TPR
    if [[ $i -eq 0 ]]; then
        gmx grompp -f ${MDP} -c ${GRO} -t ${CPT} -p ${TOP} -o step${tag}.tpr -maxwarn 10  #[1]
    else
        gmx grompp -f ${MDP} -c ${GRO} -t ${CPT} -p ${TOP} -o step${tag}.tpr -maxwarn 10
    fi
   
    # ------ Roda MD ------
    gmx mdrun -deffnm step${tag} -v || { echo "Erro no mdrun no ciclo $i"; exit 1; }


   # 3. Extrair o último quadro
    echo 0 | gmx trjconv -f step${tag}.trr -s step${tag}.tpr -dump $((SEG_NS*1000)) -o step${tag}_last.gro -pbc mol -ur compact
    
    rm -rf step${tag}.trr 
    # ------ Evaporação ------
    echo "Removendo solvente (plane_frac=$plane)..."
    python3 evaporate.py step${tag}_last.gro ${TOP} step${tag}_post.gro step${tag}_post.top ${plane} || { echo "Erro na evaporação no ciclo $i"; exit 1; }

    # ------ Atualiza para próximo ciclo ------
    GRO=step${tag}_post.gro
    TOP=step${tag}_post.top
    
    # 3. Determina altura Z do sistema
    echo ">>> Calculando nova altura da caixa..."
    zmax=$(tail -n 1 "$GRO" | awk '{print $3}')

    # 4. Proteção: garantir que Z >= 2.75 nm (para rlist padrão de 1.2 nm)
    minZ=2.75
    zsafe=$(echo "$zmax < $minZ" | bc -l)
    if [ "$zsafe" -eq 1 ]; then
        echo ">>> Altura Z calculada ($zmax nm) é muito pequena. Corrigindo para ${minZ} nm"
        zmax=$minZ
        cp ${GRO} old_${GRO}
        
                   
        # 5. Redefine caixa com nova altura
        echo ">>> Aplicando nova altura da caixa Z = ${zmax} nm"
        gmx editconf -f step${tag}_post.gro -o newbox.gro -c -box 6.0 6.0 ${minZ}
        
        # 6. Prepara para próximo ciclo
        cp newbox.gro box_fixed_${tag}.gro
        cp newbox.gro ${GRO}
        rm -rf newbox.gro
        
        echo "[INFO] Verificando consistência após ciclo $i"
        n_atoms_gro=$(tail -n +3 ${GRO} | wc -l)
        n_sol_top=$(grep "^SOL" ${TOP} | awk '{print $2}')
        echo "  Átomos no .gro: $n_atoms_gro"
        echo "  Moléculas de SOL no .top: $n_sol_top"
        python3 redist_solv.py -i box_fixed_${tag}.gro -o ${GRO} -b 6.0 6.0 ${minZ} 
        
    fi 
    
    # --- minimização
    gmx grompp -f em.mdp -c ${GRO} -p ${TOP} -o step${tag}_min.tpr -maxwarn 10
    gmx mdrun -deffnm step${tag}_min -v || { echo "Erro no EM no ciclo $i"; exit 1; }
    
    # --- novo nvt e mdp
    gmx grompp -f nvt_eq.mdp -c step${tag}_min.gro -r step${tag}_min.gro -p ${TOP} -o step${tag}_nvt.tpr -maxwarn 10
    gmx mdrun -deffnm step${tag}_nvt -v || { echo "Erro no NVT no ciclo $i"; exit 1; }

    gmx grompp -f npt_eq.mdp -c step${tag}_nvt.gro -t step${tag}_nvt.cpt -r step${tag}_nvt.gro -p ${TOP} -o step${tag}_npt.tpr -maxwarn 10
    gmx mdrun -deffnm step${tag}_npt -v || { echo "Erro no NPT no ciclo $i"; exit 1; }
    
    CPT=step${tag}_npt.cpt
    GRO=step${tag}_npt.gro
   
done

rm -rf *#*

# Faz anneling final

echo "==== Simulação completa até ciclo $CYCLES, resfriando sistema... ===="

CPT=step19_npt.cpt
GRO=step19_npt.gro
TOP=step19_post.top

 # 3. Determina altura Z do sistema
    echo ">>> Calculando nova altura da caixa..."
    zmax=$(tail -n 1 "$GRO" | awk '{print $3}')

    # 4. Proteção: garantir que Z >= 2.75 nm (para rlist padrão de 1.2 nm)
    minZ=2.75
    zsafe=$(echo "$zmax < $minZ" | bc -l)
    if [ "$zsafe" -eq 1 ]; then
        echo ">>> Altura Z calculada ($zmax nm) é muito pequena. Corrigindo para ${minZ} nm"
        zmax=$minZ
        cp ${GRO} old_${GRO}
        
                   
        # 5. Redefine caixa com nova altura
        echo ">>> Aplicando nova altura da caixa Z = ${zmax} nm"
        gmx editconf -f step19_post.gro -o newbox.gro -c -box 6.0 6.0 ${minZ}
        
        # 6. Prepara para próximo ciclo
        cp newbox.gro box_fixed_${tag}.gro
        cp newbox.gro ${GRO}
        rm -rf newbox.gro
        
        echo "[INFO] Verificando consistência após ciclo $i"
        n_atoms_gro=$(tail -n +3 ${GRO} | wc -l)
        n_sol_top=$(grep "^SOL" ${TOP} | awk '{print $2}')
        echo "  Átomos no .gro: $n_atoms_gro"
        echo "  Moléculas de SOL no .top: $n_sol_top"
    fi 

gmx grompp -f md_anneling.mdp -c step19_npt.gro -t step19_npt.cpt -p step19_post.top -o filme_final.tpr -maxwarn 10 || { echo "Erro no resfiamento"; exit 1; }
gmx mdrun -deffnm filme_final -v || { echo "Erro no resfiamento"; exit 1; }

echo "==== Simulação completa !! ===="

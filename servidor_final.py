import os
import subprocess
import io
import html
import signal
from flask import Flask, render_template, request, jsonify
from Bio import SeqIO, Phylo, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

app = Flask(__name__)

# =====================================================================
# NCBI RULE: Added cenimp email
# =====================================================================
Entrez.email = "cenimp13@gmail.com"

@app.route('/')
def home():
    return render_template('renimpblast.html')

@app.route('/rodar-script', methods=['POST'])
def rodar_script():
    dados = request.get_json(silent=True) or {}
    sequencia_bruta = dados.get('sequencia', '').strip()

    if not sequencia_bruta:
        return jsonify({"mensagem": "No sequence was submitted!"})

    try:
        # Remove any existing FASTA header if the user pasted one, keeping only the sequence bases
        linhas_seq = [
            linha.strip() for linha in sequencia_bruta.splitlines()
            if linha.strip() and not linha.strip().startswith(">")
        ]
        sequencia_limpa = "".join(linhas_seq)

        if not sequencia_limpa:
            return jsonify({"mensagem": "Invalid sequence format! Please paste the nucleotide bases."})

        with open("query_temp.fasta", "w") as f:
            f.write(">Your_Sequence\n")
            f.write(sequencia_limpa + "\n")

        hits_combinados = []

        # =========================================================
        # 1. LOCAL SEARCH (RENIMP DATABASE)
        # Added '-task blastn' and '-evalue 10' to guarantee at least 5 hits
        # =========================================================
        comando_blast = [
            "blastn", "-task", "blastn", "-query", "query_temp.fasta", "-db", "banco_peixes",
            "-evalue", "10", "-outfmt", "6 qseqid sseqid pident length evalue", "-max_target_seqs", "50"
        ]
        resultado_blast_local = subprocess.run(comando_blast, capture_output=True, text=True)
        
        if resultado_blast_local.returncode == 0 and resultado_blast_local.stdout.strip():
            linhas_local = resultado_blast_local.stdout.strip().split('\n')
            ids_locais_vistos = set()
            for linha in linhas_local:
                colunas = linha.split('\t')
                if len(colunas) < 5:
                    continue
                sseqid = colunas[1]
                if sseqid in ids_locais_vistos:
                    continue
                ids_locais_vistos.add(sseqid)

                pident = float(colunas[2])
                evalue = float(colunas[4])
                
                partes = sseqid.split('_')
                especie = f"{partes[0]}_{partes[1]}" if len(partes) >= 2 else sseqid
                
                hits_combinados.append({
                    'original_id': sseqid,
                    'especie': especie,
                    'pident': pident,
                    'evalue': evalue,
                    'source': 'LOCAL'
                })

        # =========================================================
        # 2. ONLINE SEARCH (NCBI 'nt') - With 18s Cloud Timeout
        # =========================================================
        def timeout_handler(signum, frame):
            raise TimeoutError("NCBI queue took too long")

        try:
            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(18)  # Max 18 seconds waiting for NCBI so Render never times out
            
            resultado_blast_online = NCBIWWW.qblast(
                "blastn", "nt", sequencia_limpa,
                entrez_query="txid7898[ORGN]", hitlist_size=30
            )
            blast_record = NCBIXML.read(resultado_blast_online)
            signal.alarm(0)  # Disable alarm if NCBI responded in time
            
            if blast_record.alignments:
                for alignment in blast_record.alignments:
                    accession = alignment.accession
                    titulo = alignment.hit_def
                    partes = titulo.split()
                    especie = f"{partes[0]}_{partes[1]}" if len(partes) >= 2 else partes[0]
                    
                    hsp = alignment.hsps[0]
                    evalue = float(hsp.expect)
                    pident = (hsp.identities / hsp.align_length) * 100
                    
                    hits_combinados.append({
                        'original_id': accession,
                        'especie': especie,
                        'pident': pident,
                        'evalue': evalue,
                        'source': 'NCBI'
                    })
        except Exception as e_ncbi:
            signal.alarm(0)
            print(f"Warning: NCBI online search skipped to prevent cloud timeout: {e_ncbi}")

        if not hits_combinados:
            return jsonify({"mensagem": "No results found in any database."})

        # =========================================================
        # 3. UNIFIED RANKING AND DECISION MAKING (MANDATORY 5 HITS)
        # =========================================================
        hits_ordenados = sorted(
            hits_combinados,
            key=lambda h: (h['evalue'], -h['pident'], 0 if h['source'] == 'LOCAL' else 1)
        )

        hits_unicos_final = []
        especies_vistas = set()
        ids_escolhidos = set()

        # Step 3.1: Prioritize unique species first
        for h in hits_ordenados:
            if h['especie'] not in especies_vistas:
                especies_vistas.add(h['especie'])
                ids_escolhidos.add(h['original_id'])
                hits_unicos_final.append(h)
                if len(hits_unicos_final) == 5:
                    break

        # Step 3.2: If fewer than 5 unique species exist, fill remaining slots with next best hits
        if len(hits_unicos_final) < 5:
            for h in hits_ordenados:
                if h['original_id'] not in ids_escolhidos:
                    ids_escolhidos.add(h['original_id'])
                    hits_unicos_final.append(h)
                    if len(hits_unicos_final) == 5:
                        break

        # Assign unique terminal names for IQ-TREE in case multiple hits share the same species name
        contagem_nomes = {}
        for h in hits_unicos_final:
            base_name = f"{h['source']}_{h['especie']}"
            contagem_nomes[base_name] = contagem_nomes.get(base_name, 0) + 1
            if contagem_nomes[base_name] > 1:
                h['tree_id'] = f"{base_name}_{contagem_nomes[base_name]}"
            else:
                h['tree_id'] = base_name

        tabela = ">>> UNIFIED BLAST RANKING (Priority: Local RENIMP - Top 5 Hits) <<<\n\n"
        tabela += f"{'Source':<8} | {'ID / Accession':<22} | {'Identified Species':<25} | {'Identity':<11} | {'E-value':<10}\n"
        tabela += "-" * 88 + "\n"

        for h in hits_unicos_final:
            orig = f"[{h['source']}]"
            id_str = h['original_id'][:20]
            esp = h['especie'].replace('_', ' ')[:23]
            pid = f"{h['pident']:.2f}%"
            ev = f"{h['evalue']:.2e}" if h['evalue'] < 0.001 else str(h['evalue'])
            tabela += f"{orig:<8} | {id_str:<22} | {esp:<25} | {pid:<11} | {ev:<10}\n"

        # =========================================================
        # 4. GATHER SEQUENCES FOR THE TREE
        # =========================================================
        sequencias_para_arvore = [SeqIO.read("query_temp.fasta", "fasta")]

        # 4.1 Get winners from Local database
        locais_ids = {h["original_id"] for h in hits_unicos_final if h['source'] == 'LOCAL'}
        if locais_ids:
            for registro in SeqIO.parse("renimp_database.fasta", "fasta"):
                if registro.id in locais_ids:
                    for hit in hits_unicos_final:
                        if hit['source'] == 'LOCAL' and hit['original_id'] == registro.id:
                            registro.id = hit['tree_id']
                            registro.description = ""
                            sequencias_para_arvore.append(registro)
                            break

        # 4.2 Get winners from NCBI
        ncbi_ids = [h["original_id"] for h in hits_unicos_final if h['source'] == 'NCBI']
        if ncbi_ids:
            try:
                handle = Entrez.efetch(db="nucleotide", id=ncbi_ids, rettype="fasta", retmode="text")
                registros_baixados = list(SeqIO.parse(handle, "fasta"))
                handle.close()
                
                for registro in registros_baixados:
                    acc_limpo = registro.id.split('.')[0]
                    for hit in hits_unicos_final:
                        if hit['source'] == 'NCBI' and hit['original_id'] == acc_limpo:
                            registro.id = hit['tree_id']
                            registro.description = ""
                            sequencias_para_arvore.append(registro)
                            break
            except Exception as e_fetch:
                print(f"Warning: Failed to fetch NCBI sequences for tree: {e_fetch}")

        SeqIO.write(sequencias_para_arvore, "tree_input.fasta", "fasta")

        # =========================================================
        # 5. ALIGNMENT AND ML TREE (MAFFT + IQ-TREE)
        # =========================================================
        with open("tree_aln.fasta", "w") as aln_out:
            subprocess.run(
                ["mafft", "--auto", "--adjustdirection", "--quiet", "tree_input.fasta"],
                stdout=aln_out,
                check=True
            )

        iqtree_bin = "iqtree"
        if subprocess.run(["which", "iqtree"], capture_output=True).returncode != 0:
            iqtree_bin = "iqtree2"

        subprocess.run(
            [iqtree_bin, "-s", "tree_aln.fasta", "-m", "GTR+G", "-fast", "-nt", "1", "-redo"],
            capture_output=True,
            check=True
        )

        tree = Phylo.read("tree_aln.fasta.treefile", "newick")

        # Root with the last sequence in the list (most distant from ranking)
        outgroup_hit = hits_unicos_final[-1]
        outgroup_prefix = outgroup_hit['tree_id']
        for clado in tree.get_terminals():
            if clado.name and clado.name.replace("_R_", "") == outgroup_prefix:
                try:
                    tree.root_with_outgroup({"name": clado.name})
                except Exception:
                    pass
                break

        # =========================================================
        # 6. DECISION (>= 98% THRESHOLD) AND TIE-BREAKER
        # =========================================================
        query_node = next(c for c in tree.get_terminals() if c.name and "Your_Sequence" in c.name)
        vizinhos_mais_proximos = []
        menor_distancia = float('inf')
        tolerancia = 1e-4  # Margem para considerar empate no comprimento do ramo (IQ-TREE)

        for clade in tree.get_terminals():
            if clade != query_node:
                dist = tree.distance(query_node, clade)
                if dist < menor_distancia - tolerancia:
                    menor_distancia = dist
                    vizinhos_mais_proximos = [clade.name]
                elif abs(dist - menor_distancia) <= tolerancia:
                    vizinhos_mais_proximos.append(clade.name)

        # Coleta dados de todos os empatados
        especies_empatadas = []
        origens_empatadas = []
        pident_vencedor = 0.0

        for vizinho in vizinhos_mais_proximos:
            nome_terminal_limpo = vizinho.replace("_R_", "")
            hit_vencedor = next(
                (h for h in hits_unicos_final if h['tree_id'] == nome_terminal_limpo),
                hits_unicos_final[0]
            )
            especie_nome = hit_vencedor['especie'].replace('_', ' ')
            if especie_nome not in especies_empatadas:
                especies_empatadas.append(especie_nome)
            
            origem = "LOCAL (RENIMP Database)" if hit_vencedor['source'] == "LOCAL" else "NCBI (Online)"
            if origem not in origens_empatadas:
                origens_empatadas.append(origem)
            
            if hit_vencedor['pident'] > pident_vencedor:
                pident_vencedor = hit_vencedor['pident']

        # Formata o texto final com suporte a empates
        if len(especies_empatadas) > 1:
            especie_final_str = " AND ".join(especies_empatadas)
            origem_final_str = " & ".join(origens_empatadas)
            alerta_empate = " [TIE IN PHYLOGENY]"
        else:
            especie_final_str = especies_empatadas[0]
            origem_final_str = origens_empatadas[0]
            alerta_empate = ""

        if pident_vencedor >= 98.0:
            decisao_final = (
                f"[ANALYSIS COMPLETED] - Final Verdict: The identified species is '{especie_final_str}'{alerta_empate} "
                f"({pident_vencedor:.2f}% identity).\n"
                f"The algorithm extracted this hit from source: {origem_final_str}."
            )
            relatorio_especie = f"Species: {especie_final_str} (Confirmed >= 98%)\n"
        else:
            decisao_final = (
                f"[INCONCLUSIVE VERDICT] - Could not reach a final verdict "
                f"(highest identity is {pident_vencedor:.2f}%, below the 98.0% threshold).\n"
                f"Closest evolutionary neighbor in the tree: '{especie_final_str}'{alerta_empate} (Source: {origem_final_str})."
            )
            relatorio_especie = f"Closest Neighbor (Unconfirmed - {pident_vencedor:.2f}%): {especie_final_str}\n"

        # Link do FishBase (usa a primeira espécie em caso de empate)
        try:
            partes_esp = especies_empatadas[0].split()
            genero, epiteto = partes_esp[0], partes_esp[1]
            link_fishbase = f"https://www.fishbase.se/summary/SpeciesSummary.php?genusname={genero}&speciesname={epiteto}"
        except Exception:
            link_fishbase = "Invalid name format for FishBase."

        # Omitir rótulos internos na árvore e desenhar ASCII
        for clado in tree.find_clades():
            if clado.name and not clado.is_terminal():
                clado.name = ""
        
        out_handle = io.StringIO()
        Phylo.draw_ascii(tree, file=out_handle)
        arvore_ascii = out_handle.getvalue()

        # =========================================================
        # 7. BUILD VISUAL HTML RESPONSE
        # =========================================================
        resposta_final = tabela + "\n"
        resposta_final += ">>> ALGORITHM DIAGNOSTIC (Decision Support) <<<\n"
        resposta_final += decisao_final + "\n\n"
        resposta_final += ">>> BIOLOGICAL REPORT OF THE SPECIES <<<\n"
        resposta_final += relatorio_especie
        resposta_final += f"More details on FishBase: {link_fishbase}\n\n"
        resposta_final += ">>> PHYLOGENETIC TREE (IQ-TREE) <<<\n"
        resposta_final += arvore_ascii

        resposta_html = html.escape(resposta_final)

        resposta_html = resposta_html.replace(
            "_R_Your_Sequence",
            "<span style='color: #fbbf24; font-weight: bold; text-decoration: underline; font-size: 1.1em;'>Your_Sequence (Reverse)</span>"
        )
        resposta_html = resposta_html.replace(
            "Your_Sequence",
            "<span style='color: #fbbf24; font-weight: bold; text-decoration: underline; font-size: 1.1em;'>Your_Sequence</span>"
        )
        resposta_html = resposta_html.replace(
            "[ANALYSIS COMPLETED]",
            "<span style='color: #4ade80; font-weight: bold;'>[ANALYSIS COMPLETED]</span>"
        )
        resposta_html = resposta_html.replace(
            "[INCONCLUSIVE VERDICT]",
            "<span style='color: #f87171; font-weight: bold;'>[INCONCLUSIVE VERDICT]</span>"
        )
        resposta_html = resposta_html.replace(
            "[LOCAL]",
            "<span style='color: #a78bfa; font-weight: bold;'>[LOCAL]</span>"
        )
        resposta_html = resposta_html.replace(
            "[NCBI]",
            "<span style='color: #38bdf8; font-weight: bold;'>[NCBI]</span>"
        )
        resposta_html = resposta_html.replace(
            "LOCAL_",
            "<span style='color: #a78bfa;'>LOCAL_</span>"
        )
        resposta_html = resposta_html.replace(
            "NCBI_",
            "<span style='color: #38bdf8;'>NCBI_</span>"
        )
        if link_fishbase.startswith("http"):
            link_escapado = html.escape(link_fishbase)
            resposta_html = resposta_html.replace(
                link_escapado,
                f"<a href='{link_escapado}' target='_blank' style='color: #60a5fa; text-decoration: underline;'>{link_escapado}</a>"
            )

        return jsonify({"mensagem": resposta_html})

    except Exception as e:
        return jsonify({"mensagem": "Error in decision pipeline: " + str(e)})

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port, debug=False)

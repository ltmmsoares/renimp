import os
import subprocess
import io
import html
import requests
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
    dados = request.get_json()
    sequencia = dados.get('sequencia', '').strip()

    if not sequencia:
        return jsonify({"mensagem": "No sequence was submitted!"})

    try:
        with open("query_temp.fasta", "w") as f:
            f.write(">Your_Sequence\n")
            f.write(sequencia + "\n")

        hits_combinados = []

        # =========================================================
        # 1. LOCAL SEARCH (RENIMP DATABASE)
        # =========================================================
        comando_blast = [
            "blastn", "-query", "query_temp.fasta", "-db", "banco_peixes",
            "-outfmt", "6 qseqid sseqid pident length evalue", "-max_target_seqs", "20"
        ]
        resultado_blast_local = subprocess.run(comando_blast, capture_output=True, text=True)
        
        if resultado_blast_local.returncode == 0 and resultado_blast_local.stdout.strip():
            linhas_local = resultado_blast_local.stdout.strip().split('\n')
            for linha in linhas_local:
                colunas = linha.split('\t')
                if len(colunas) < 5: continue
                sseqid = colunas[1]
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
        # 2. ONLINE SEARCH (NCBI 'nt')
        # =========================================================
        resultado_blast_online = NCBIWWW.qblast("blastn", "nt", sequencia, entrez_query="txid7898[ORGN]", hitlist_size=20)
        blast_record = NCBIXML.read(resultado_blast_online)
        
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

        # =========================================================
        # 3. UNIFIED RANKING AND DECISION MAKING (THE SECRET!)
        # =========================================================
        # We sort all hits by: 1st E-value (lower is better), 2nd Identity (higher is better)
        # "0 if h['source'] == 'LOCAL' else 1" ensures that in case of a tie, LOCAL wins!
        hits_ordenados = sorted(hits_combinados, key=lambda h: (h['evalue'], -h['pident'], 0 if h['source'] == 'LOCAL' else 1))

        hits_unicos_final = []
        especies_vistas = set()

        tabela = ">>> UNIFIED BLAST RANKING (Priority: Local RENIMP) <<<\n\n"
        tabela += f"{'Source':<8} | {'ID / Accession':<22} | {'Identified Species':<25} | {'Identity':<11} | {'E-value':<10}\n"
        tabela += "-" * 88 + "\n"

        for h in hits_ordenados:
            if h['especie'] not in especies_vistas:
                especies_vistas.add(h['especie'])
                hits_unicos_final.append(h)
                
                orig = f"[{h['source']}]"
                id_str = h['original_id'][:20]
                esp = h['especie'].replace('_', ' ')[:23]
                pid = f"{h['pident']:.2f}%"
                ev = f"{h['evalue']:.2e}" if h['evalue'] < 0.001 else str(h['evalue'])
                
                tabela += f"{orig:<8} | {id_str:<22} | {esp:<25} | {pid:<11} | {ev:<10}\n"
                
                if len(hits_unicos_final) == 5: # Limited to Top 5 best unique hits
                    break

        if not hits_unicos_final:
            return jsonify({"mensagem": "No results found in any database."})

        # =========================================================
        # 4. GATHER SEQUENCES FOR THE TREE
        # =========================================================
        sequencias_para_arvore = [SeqIO.read("query_temp.fasta", "fasta")]

        # 4.1 Get winners from Local database
        locais_ids = {h["original_id"] for h in hits_unicos_final if h['source'] == 'LOCAL'}
        if locais_ids:
            for registro in SeqIO.parse("renimp_database.fasta", "fasta"):
                if registro.id in locais_ids:
                    # Find correct species in dictionary to rename
                    hit = next(h for h in hits_unicos_final if h['original_id'] == registro.id)
                    registro.id = f"LOCAL_{hit['especie']}"
                    registro.description = ""
                    sequencias_para_arvore.append(registro)

        # 4.2 Get winners from NCBI
        ncbi_ids = [h["original_id"] for h in hits_unicos_final if h['source'] == 'NCBI']
        if ncbi_ids:
            handle = Entrez.efetch(db="nucleotide", id=ncbi_ids, rettype="fasta", retmode="text")
            registros_baixados = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            
            for registro in registros_baixados:
                acc_limpo = registro.id.split('.')[0]
                for hit in hits_unicos_final:
                    if hit['source'] == 'NCBI' and hit['original_id'] == acc_limpo:
                        registro.id = f"NCBI_{hit['especie']}"
                        registro.description = ""
                        sequencias_para_arvore.append(registro)
                        break

        SeqIO.write(sequencias_para_arvore, "tree_input.fasta", "fasta")

        # =========================================================
        # 5. ALIGNMENT AND ML TREE (MAFFT + IQ-TREE)
        # =========================================================
        subprocess.run(["mafft", "--auto", "--adjustdirection", "--quiet", "tree_input.fasta"], stdout=open("tree_aln.fasta", "w"))
        subprocess.run(["iqtree", "-s", "tree_aln.fasta", "-redo"], capture_output=True)

        tree = Phylo.read("tree_aln.fasta.treefile", "newick")

        # Root with the last sequence in the list (most distant from ranking)
        outgroup_hit = hits_unicos_final[-1]
        outgroup_prefix = f"{outgroup_hit['source']}_{outgroup_hit['especie']}"
        for clado in tree.get_terminals():
            if clado.name.startswith(outgroup_prefix):
                try: tree.root_with_outgroup({"name": clado.name})
                except: pass
                break

        # =========================================================
        # 6. DECISION, WIKIPEDIA AND FISHBASE
        # =========================================================
        query_node = next(c for c in tree.get_terminals() if "Your_Sequence" in c.name)
        vizinho_mais_proximo = None
        menor_distancia = float('inf')

        for clade in tree.get_terminals():
            if clade != query_node:
                dist = tree.distance(query_node, clade)
                if dist < menor_distancia:
                    menor_distancia = dist
                    vizinho_mais_proximo = clade.name

        origem_vizinho = "LOCAL (RENIMP Database)" if vizinho_mais_proximo.startswith("LOCAL") else "NCBI (Online)"
        
        nome_limpo = vizinho_mais_proximo.replace("LOCAL_", "").replace("NCBI_", "")
        especie_final = nome_limpo.replace('_', ' ')

        decisao_final = f"[ANALYSIS COMPLETED] - The closest genetic neighbor is '{especie_final}'.\n"
        decisao_final += f"The algorithm extracted this hit from source: {origem_vizinho}."

        try:
            genero, epiteto = especie_final.split()[0], especie_final.split()[1]
            link_fishbase = f"https://www.fishbase.se/summary/SpeciesSummary.php?genusname={genero}&speciesname={epiteto}"
        except:
            link_fishbase = "Invalid name format for FishBase."

        # Note: Using English Wikipedia endpoint for English translation
        url_wiki = f"https://en.wikipedia.org/api/rest_v1/page/summary/{especie_final.replace(' ', '_')}"
        try:
            req = requests.get(url_wiki)
            relato_biologico = req.json().get('extract', 'Summary not available in the open database.') if req.status_code == 200 else "Detailed report not found."
        except:
            relato_biologico = "Connection error while trying to fetch the report."

        for clado in tree.find_clades():
            if clado.name and not clado.is_terminal(): clado.name = ""

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
        resposta_final += f"Species: {especie_final}\n"
        resposta_final += f"About: {relato_biologico}\n"
        resposta_final += f"More details on FishBase: {link_fishbase}\n\n"
        resposta_final += ">>> PHYLOGENETIC TREE (IQ-TREE) <<<\n"
        resposta_final += arvore_ascii

        resposta_html = html.escape(resposta_final)

        resposta_html = resposta_html.replace("Your_Sequence", "<span style='color: #fbbf24; font-weight: bold; text-decoration: underline; font-size: 1.1em;'>Your_Sequence</span>")
        resposta_html = resposta_html.replace("_R_Your_Sequence", "<span style='color: #fbbf24; font-weight: bold; text-decoration: underline; font-size: 1.1em;'>Your_Sequence (Reverse)</span>")
        resposta_html = resposta_html.replace("[ANALYSIS COMPLETED]", "<span style='color: #4ade80; font-weight: bold;'>[ANALYSIS COMPLETED]</span>")
        resposta_html = resposta_html.replace("[LOCAL]", "<span style='color: #a78bfa; font-weight: bold;'>[LOCAL]</span>")
        resposta_html = resposta_html.replace("[NCBI]", "<span style='color: #38bdf8; font-weight: bold;'>[NCBI]</span>")
        resposta_html = resposta_html.replace("LOCAL_", "<span style='color: #a78bfa;'>LOCAL_</span>")
        resposta_html = resposta_html.replace("NCBI_", "<span style='color: #38bdf8;'>NCBI_</span>")
        resposta_html = resposta_html.replace(link_fishbase, f"<a href='{link_fishbase}' target='_blank' style='color: #60a5fa; text-decoration: underline;'>{link_fishbase}</a>")

        return jsonify({"mensagem": resposta_html})

    except Exception as e:
        return jsonify({"mensagem": "Error in decision pipeline: " + str(e)})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
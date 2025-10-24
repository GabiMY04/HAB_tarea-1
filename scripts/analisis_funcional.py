
####################################################################
############################## LIBRERIAS ###########################
####################################################################

import mygene 
import requests, os, json
import pandas as pd
from pathlib import Path

#####################################################################
############################### FUNCIONES ###########################
#####################################################################


#########################################
######## leer_txt #######################
#########################################

def leer_txt(ruta):
    """
    Leer la ruta de archivo de texto y devuelve una lista de símbolos.
    """
    with open(ruta, 'r', encoding='utf-8') as archivo: 
        text = archivo.read().strip()
    if not text:
        return []
    partes = [p.strip() for p in text.replace('\n', ',').split(',')]
    return [p for p in partes if p]


#########################################
######## anyadir_prefijo ################
#########################################

def anyadir_prefijo(g):
    """
    Añade prefijo 'MT-' a genes mitocondriales en caso de que esten en la lista
    """
    g = (g or '').strip()
    if not g:
        return g
    up = g.upper()
    mito = {'ATP6','ATP8','ND1','ND2','ND3','ND4','ND4L','ND5','ND6',
            'COX1','COX2','COX3','CYB'}  # genes mitocondriales comunes
    if up.startswith('MT-') or up.startswith('MT'):
        return up if up.startswith('MT-') else 'MT-' + up[2:]
    if up in mito:
        return 'MT-' + up
    return up




#########################################
######## mapear_genes ###################
#########################################  

def mapear_genes(symbols, species='human'):
    """
    Mapea los simbolos de genes usando mygene y devuelve los simbolos oficiales no-pseudo ordenada y sin duplicados
    """
    mg = mygene.MyGeneInfo()
    queries = [anyadir_prefijo(s) for s in symbols if s and s.strip()]
    if not queries:
        return []
    results = mg.querymany(queries, scopes='symbol', fields='symbol,type_of_gene', species=species)
    mapped = []
    for entry in results:
        if not entry.get('notfound'):
            gene_type = (entry.get('type_of_gene') or '').lower()
            sym = entry.get('symbol') or (entry.get('hits') and entry['hits'][0].get('symbol'))
            # anyadir solo si hay símbolo y no es pseudo
            if sym and ('pseudo' not in gene_type):
                mapped.append(sym)
    # quitar duplicados preservando orden
    return list(dict.fromkeys(mapped))


#########################################
######## guardar_resultados #############
#########################################  



def guardar_resultados(results, prefix='string_enrichment'):
    """
    Guarda results (lista de dicts) como CSV y JSON en ./results (cwd).
    Devuelve rutas (csv_path, json_path) o (None, None) en caso de error.
    """
    try:
        outdir = Path.cwd() / 'results'
        outdir.mkdir(parents=True, exist_ok=True)
        rows = [{
            "term": r.get("term"),
            "genes": ",".join(r.get("preferredNames", [])),
            "fdr": r.get("fdr"),
            "description": r.get("description", "")
        } for r in results]
        df = pd.DataFrame(rows)
        csv_path = outdir / f"{prefix}.csv"
        df.to_csv(csv_path, index=False)
        json_path = outdir / f"{prefix}_raw.json"
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump([r.get("raw", r) for r in results], fh, ensure_ascii=False, indent=2)
        print(f"Saved CSV -> {csv_path}")
        print(f"Saved RAW JSON -> {json_path}")
        return str(csv_path), str(json_path)
    except Exception as e:
        print("Error saving results:", e)
        return None, None

#########################################
######## consultar_genes ################
#########################################  

def consultar_genes(genes, species=9606, caller_identity='test_HAB'):
    """
    Consulta STRING enrichment para los genes y devuelve terminos GO
    """
    if not genes:
        print("No se han introducido los genes correctamente")
        return []

    url = "https://version-11-5.string-db.org/api/json/enrichment"
    params = {
        "identifiers": "\n".join(genes),
        "species": species,
        "caller_identity": caller_identity
    }

    try:
        resp = requests.post(url, data=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except requests.RequestException as e:
        print(f"Error en la petición a STRING: {e}")
        return []

    results = []
    for row in data:
        if row.get("category") == "Process":
            try:
                fdr = float(row.get("fdr", 1.0))
            except Exception:
                fdr = 1.0
            results.append({
                "term": row.get("term"),
                "preferredNames": row.get("preferredNames", []),
                "fdr": fdr,
                "description": row.get("description"),
                "category": row.get("category"),
                "raw": row
            })

    if not results:
        print("No se encontraron términos 'Process' para los genes solicitados")
        return []

    for r in results:
        names = ",".join(r["preferredNames"]) if r["preferredNames"] else ""
        print(f"{r['term']}\t{names}\t{r['fdr']:.3f}\t{r['category']} {r['description']}")

    guardar_resultados(results)
    return results

###################################################
######################## MAIN #####################
###################################################

def main():
    # Leemos el archivo con los nombres de los genes
    genes_symbols = leer_txt('data/genes_input.txt')

    #Añadimos los prefijos y normalizamos
    genes_symbols= [anyadir_prefijo(g) for g in genes_symbols if anyadir_prefijo(g)]
    print(f"Genes normalizados: {genes_symbols}")

    # Mapeo de los simbolos quitando los pseudogenes
    mapp_simbols = mapear_genes(genes_symbols, species='human')
    print("Símbolos mapeados (no-pseudo):", mapp_simbols)

    # Analisis STRING para obtener terminos GO 
    consultar_genes(mapp_simbols, species=9606, caller_identity='analisis_funcional_script')

if __name__ == "__main__":
    main()
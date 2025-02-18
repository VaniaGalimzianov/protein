from flask import Flask, render_template, request, jsonify, send_file
import requests
import pandas as pd
import os

app = Flask(name)

BASE_DIR = "generated_files"
os.makedirs(BASE_DIR, exist_ok=True)

def fetch_protein_data(protein_id):
    fasta_url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta"
    json_url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json?fields=ft_var_seq,ft_variant,ft_non_cons,ft_non_std,ft_non_ter,ft_conflict,ft_unsure,ft_act_site,ft_binding,ft_dna_bind,ft_site,ft_mutagen,ft_intramem,ft_topo_dom,ft_transmem,ft_chain,ft_crosslnk,ft_disulfid,ft_carbohyd,ft_init_met,ft_lipid,ft_mod_res,ft_peptide,ft_propep,ft_signal,ft_transit,ft_strand,ft_helix,ft_turn,ft_coiled,ft_compbias,ft_domain,ft_motif,ft_region,ft_repeat,ft_zn_fing"

    fasta_response = requests.get(fasta_url)
    json_response = requests.get(json_url)

    if fasta_response.status_code != 200 or json_response.status_code != 200:
        return None, None

    return fasta_response.text, json_response.json()

def analyze_protein(fasta, json_data):
    lines = fasta.strip().split("\n")
    sequence = "".join(lines[1:])  # Убираем заголовок FASTA

    structure_data = {"HELIX": [], "STRAND": [], "TURN": []}

    for feature in json_data.get("features", []):
        if "location" in feature and "start" in feature["location"] and "end" in feature["location"]:
            start = int(feature["location"]["start"])
            end = int(feature["location"]["end"])
            feature_type = feature["type"]

            if feature_type == "HELIX" and (end - start + 1) >= 6:
                structure_data["HELIX"].extend(sequence[start-1:end])
            elif feature_type == "STRAND" and (end - start + 1) >= 5:
                structure_data["STRAND"].extend(sequence[start-1:end])

    missing_start = 1
    first_structure_start = min(
        [int(feature["location"]["start"]) for feature in json_data.get("features", [])],
        default=len(sequence)
    )

    if first_structure_start > missing_start:
        structure_data["TURN"].extend(sequence[missing_start-1:first_structure_start-1])

    unique_amino_acids = sorted(set(sequence))
    df = pd.DataFrame(index=unique_amino_acids, columns=["HELIX", "STRAND", "TURN", "TOTAL"])

    for aa in unique_amino_acids:
        df.loc[aa, "HELIX"] = structure_data["HELIX"].count(aa)
        df.loc[aa, "STRAND"] = structure_data["STRAND"].count(aa)
        df.loc[aa, "TURN"] = structure_data["TURN"].count(aa)
        df.loc[aa, "TOTAL"] = df.loc[aa, ["HELIX", "STRAND", "TURN"]].sum()

    file_path = os.path.join(BASE_DIR, "protein_analysis.xlsx")
    df.to_excel(file_path)

    return file_path

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/analyze/<protein_id>")
def analyze(protein_id):
    fasta, json_data = fetch_protein_data(protein_id)
    
    if not fasta or not json_data:
        return jsonify({"success": False})

    file_path = analyze_protein(fasta, json_data)
    
    if file_path:
        return jsonify({"success": True})
    else:
        return jsonify({"success": False})

@app.route("/download")
def download():
    file_path = os.path.join(BASE_DIR, "protein_analysis.xlsx")
    return send_file(file_path, as_attachment=True)

if name == "main":
    app.run(debug=True)

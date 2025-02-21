from flask import Flask, render_template, request, jsonify, send_file
import requests
import pandas as pd
import os

app = Flask(__name__, template_folder='templates')
BASE_DIR = "generated_files"
os.makedirs(BASE_DIR, exist_ok=True)

def fetch_protein_data(protein_id):
    fasta_url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta"
    json_url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json?fields=ft_strand,ft_helix,ft_turn,ft_coiled"

    fasta_response = requests.get(fasta_url)
    json_response = requests.get(json_url)
    
    if fasta_response.status_code != 200 or json_response.status_code != 200:
        return None, None
    
    return fasta_response.text, json_response.json()

def analyze_protein(protein_id, fasta, json_data):
    lines = fasta.strip().split("\n")
    sequence = "".join(lines[1:])  # Убираем заголовок FASTA
    
    structure_data = {"HELIX": [], "STRAND": [], "TURN": [], "COILED": []}
    helix_positions, strand_positions, coiled_positions = set(), set(), set()

    for feature in json_data.get("features", []):
        feature_type = feature.get("type", "").upper()
        if feature_type == "COILED COIL":
            feature_type = "COILED"
        elif feature_type == "BETA STRAND":
            feature_type = "STRAND"
        
        if feature_type in structure_data:
            start = feature["location"]["start"]["value"]
            end = feature["location"]["end"]["value"]
            
            if 1 <= start <= len(sequence) and 1 <= end <= len(sequence):
                fragment = sequence[start - 1:end]
                
                if feature_type == "HELIX" and len(fragment) >= 6:
                    structure_data["HELIX"].extend(fragment)
                    helix_positions.update(range(start, end + 1))
                elif feature_type == "STRAND" and len(fragment) >= 5:
                    structure_data["STRAND"].extend(fragment)
                    strand_positions.update(range(start, end + 1))
                elif feature_type == "COILED":
                    structure_data["COILED"].extend(fragment)
                    coiled_positions.update(range(start, end + 1))
    
    all_positions = set(range(1, len(sequence) + 1))
    turn_positions = all_positions - helix_positions - strand_positions - coiled_positions
    structure_data["TURN"] = [sequence[pos - 1] for pos in turn_positions if 1 <= pos <= len(sequence)]
    
    unique_amino_acids = sorted(set(sequence))
    columns = [key for key in structure_data if structure_data[key]] + ["TOTAL"]
    df = pd.DataFrame(index=unique_amino_acids, columns=columns)

    for aa in unique_amino_acids:
        for key in structure_data:
            if key in df.columns:
                df.loc[aa, key] = structure_data[key].count(aa)
        df.loc[aa, "TOTAL"] = sum(df.loc[aa, key] for key in structure_data if key in df.columns)
    
    for key in ["HELIX", "STRAND"]:
        total_count = sum(df[key].dropna()) if key in df.columns else 0
        freq_col = f"{key}_FREQ"
        if key in df.columns and total_count > 0:
            df[freq_col] = df[key] / total_count
    
    if "HELIX_FREQ" in df.columns and "STRAND_FREQ" in df.columns:
        df["HELIX_STRAND_RATIO"] = df["HELIX_FREQ"].fillna(0) / df["STRAND_FREQ"].replace(0, float('nan'))
    
    df.dropna(how='all', axis=1, inplace=True)
    
    file_path = os.path.join(BASE_DIR, f"{protein_id}_analysis.xlsx")
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

    file_path = analyze_protein(protein_id, fasta, json_data)
    
    if file_path:
        return jsonify({"success": True, "file": file_path})
    else:
        return jsonify({"success": False})

@app.route("/download/<protein_id>")
def download(protein_id):
    file_path = os.path.join(BASE_DIR, f"{protein_id}_analysis.xlsx")
    if os.path.exists(file_path):
        return send_file(file_path, as_attachment=True)
    else:
        return "File not found", 404

if __name__ == "__main__":
    app.run(debug=True)

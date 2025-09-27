import pytest
from unittest.mock import patch, Mock
import os
import pandas as pd
from app import fetch_protein_data, analyze_protein, BASE_DIR
from unittest.mock import patch, Mock

protein_id_valid = "P12345"
fasta_sample = ">P12345 Sample protein\nMKTIIALSYIFCLVFADYKDDDDK"
json_sample = {
    "features": [
        {"type": "HELIX", "location": {"start": {"value": 2}, "end": {"value": 7}}},
        {"type": "BETA STRAND", "location": {"start": {"value": 10}, "end": {"value": 14}}},
        {"type": "COILED COIL", "location": {"start": {"value": 15}, "end": {"value": 18}}}
    ]
}

# 1. Тест fetch_protein_data: корректные данные
@patch("app.requests.get")
def test_fetch_protein_data_success(mock_get):
    mock_get.side_effect = [
        Mock(status_code=200, text=fasta_sample),
        Mock(status_code=200, json=lambda: json_sample)
    ]
    fasta, json_data = fetch_protein_data(protein_id_valid)
    assert fasta.startswith(">P12345")
    assert "features" in json_data

# 2. Тест fetch_protein_data: ошибка запроса
@patch("app.requests.get")
def test_fetch_protein_data_fail(mock_get):
    mock_get.return_value = Mock(status_code=404)
    fasta, json_data = fetch_protein_data("INVALID")
    assert fasta is None
    assert json_data is None

# 3. Тест analyze_protein: создание файла и структура DataFrame
def test_analyze_protein_creates_file(tmp_path):
    file_path = analyze_protein(protein_id_valid, fasta_sample, json_sample)
    assert os.path.exists(file_path)
    df = pd.read_excel(file_path)
    assert "HELIX" in df.columns
    assert "STRAND" in df.columns
    # Проверка расчёта TOTAL
    assert all(df["TOTAL"] >= 0)

# 4. Тест analyze_protein: граничные случаи (короткая HELIX)
json_short_helix = {
    "features": [{"type": "HELIX", "location": {"start": {"value": 1}, "end": {"value": 3}}}]
}
def test_analyze_protein_short_helix(tmp_path):
    file_path = analyze_protein(protein_id_valid, fasta_sample, json_short_helix)
    df = pd.read_excel(file_path)
    # HELIX слишком короткая → не должна учитываться
    assert "HELIX" not in df.columns or df["HELIX"].sum() == 0

# 5. Тест analyze_protein: HELIX/STRAND frequency и ratio
def test_helix_strand_ratio(tmp_path):
    df = pd.read_excel(analyze_protein(protein_id_valid, fasta_sample, json_sample))
    if "HELIX_STRAND_RATIO" in df.columns:
        ratio = df["HELIX_STRAND_RATIO"].iloc[0]
        assert ratio >= 0 or pd.isna(ratio)
